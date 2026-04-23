import time
import threading
import multiprocessing as mp
import traceback
import sys
import queue as pyqueue
import itertools
import os
import tempfile
import shutil
from multiprocessing.connection import wait

try:
    import cplex
except Exception:
    cplex = None

try:
    import gurobipy as gp
except Exception:
    gp = None

from pysat.solvers import Solver
from solvers.backend import (
    EXTERNAL_SOLVERS,
    run_external_solver,
    cpu_self_seconds,
    solve_radius_cplex,
    solve_radius_gurobi,
)


def _safe_put(q, item):
    try:
        q.put(item)
    except Exception:
        pass


def _safe_close_queue_like(q):
    if q is None:
        return
    try:
        q.close()
    except Exception:
        pass
    try:
        q.join_thread()
    except Exception:
        pass


def _safe_close_conn(conn):
    if conn is None:
        return
    try:
        conn.close()
    except Exception:
        pass


def _terminate_process_fast(proc, *, grace=0.02):
    if proc is None:
        return
    try:
        if not proc.is_alive():
            return
    except Exception:
        return

    try:
        proc.terminate()
    except Exception:
        pass

    try:
        proc.join(timeout=grace)
    except Exception:
        pass

    try:
        if proc.is_alive():
            proc.kill()
            proc.join(timeout=grace)
    except Exception:
        pass


def _recv_if_matching(conn, step):
    """
    Drain exactly one message if available and return it only if step matches.
    Non-matching stale messages are ignored.
    """
    try:
        if conn.poll():
            item = conn.recv()
        else:
            return None
    except EOFError:
        return None
    except Exception:
        return None

    if not isinstance(item, tuple) or len(item) != 3:
        return None

    tag, step_got, payload = item
    if step_got != step:
        return None
    return (tag, payload)


def _check_mip_backend_available(mip_backend: str):
    if mip_backend == "cplex_mip":
        if cplex is None:
            raise ImportError("hybrid_race with mip_backend='cplex_mip' needs CPLEX installed.")
    elif mip_backend == "gurobi_mip":
        if gp is None:
            raise ImportError("hybrid_race with mip_backend='gurobi_mip' needs Gurobi installed.")
    else:
        raise ValueError(f"Unsupported mip_backend: {mip_backend}")


def _solve_radius_external_once(inst, idx, encoding, solver_name, radius, time_limit, cancel_ev=None):
    cnf, varmap = inst._encode_cnf(radius, encoding)

    if cnf is None:
        return idx, radius, "unsat", 0.0, None, None, None

    base_tmp = "/dev/shm" if os.path.isdir("/dev/shm") else None
    local_tmpdir = tempfile.mkdtemp(prefix="pcsat_race_", dir=base_tmp)

    try:
        status, cpu_sec, model = run_external_solver(
            solver_name=solver_name,
            cnf=cnf,
            time_limit=time_limit,
            cancel_ev=cancel_ev,
            tmpdir=local_tmpdir,
        )

        nvars = cnf.nv
        nclauses = len(cnf.clauses)

        if status in ("timeout", "error"):
            return idx, radius, status, cpu_sec, nvars, nclauses, None

        if status == "unsat":
            return idx, radius, "unsat", cpu_sec, nvars, nclauses, None

        model_set = set(model or [])
        y_vars = varmap.get("y", [])
        candidates = set(varmap.get("candidates", []))
        chosen = {
            j for j, v in enumerate(y_vars)
            if (v in model_set) and (j in candidates)
        }
        centers = sorted(chosen)

        return idx, radius, "sat", cpu_sec, nvars, nclauses, centers

    finally:
        try:
            shutil.rmtree(local_tmpdir, ignore_errors=True)
        except Exception:
            pass


# -------------------------------
# Persistent external SAT worker
# -------------------------------

def _external_sat_control_loop(control_q, state):
    while True:
        msg = control_q.get()
        if msg is None:
            return

        cmd = msg.get("cmd")
        if cmd == "stop":
            with state["lock"]:
                state["stopped"] = True
                ev = state["current_cancel_ev"]
                if ev is not None:
                    ev.set()
            return

        if cmd == "cancel":
            step = msg.get("step")
            with state["lock"]:
                state["pending_cancel_steps"].add(step)
                if state["current_step"] == step and state["current_cancel_ev"] is not None:
                    state["current_cancel_ev"].set()


def _external_sat_worker_loop(inst, encoding, solver_name, job_q, control_q, result_send, time_limit):
    state = {
        "lock": threading.Lock(),
        "pending_cancel_steps": set(),
        "current_step": None,
        "current_cancel_ev": None,
        "stopped": False,
    }

    ctl_th = threading.Thread(
        target=_external_sat_control_loop,
        args=(control_q, state),
        daemon=True,
    )
    ctl_th.start()

    while True:
        with state["lock"]:
            if state["stopped"]:
                break

        msg = job_q.get()
        if msg is None:
            break
        if msg.get("cmd") != "solve":
            continue

        step = msg["step"]
        idx = msg["idx"]
        radius = msg["radius"]

        with state["lock"]:
            if step in state["pending_cancel_steps"]:
                state["pending_cancel_steps"].discard(step)
                result_send.send(("cancelled", step, None))
                continue

            cancel_ev = threading.Event()
            state["current_step"] = step
            state["current_cancel_ev"] = cancel_ev

        try:
            res = _solve_radius_external_once(
                inst=inst,
                idx=idx,
                encoding=encoding,
                solver_name=solver_name,
                radius=radius,
                time_limit=time_limit,
                cancel_ev=cancel_ev,
            )

            cancelled = cancel_ev.is_set()
            if cancelled and res[2] not in ("sat", "unsat"):
                result_send.send(("cancelled", step, None))
            else:
                result_send.send(("result", step, res))

        except Exception:
            result_send.send(("error", step, traceback.format_exc()))

        finally:
            with state["lock"]:
                state["current_step"] = None
                state["current_cancel_ev"] = None

    _safe_close_conn(result_send)


# -------------------------------
# Persistent incremental SAT worker
# -------------------------------

class _PersistentIncrementalSATWorker:
    def __init__(self, inst, encoding, solver_name, time_limit):
        self.inst = inst
        self.encoding = encoding
        self.solver_name = solver_name
        self.time_limit = time_limit

        base_cnf, info = inst._build_incremental_base_cnf(encoding=encoding)
        self.base_cnf = base_cnf
        self.radii = info["radii"]
        self.y_vars = info["y"]

        self.next_var = base_cnf.nv + 1
        self.tested = {}
        self.added_clause_count = 0

        self.solver = Solver(name=solver_name, bootstrap_with=base_cnf.clauses)

        self.job_q = pyqueue.Queue()

        self.result_recv, self.result_send = mp.Pipe(duplex=False)

        self._lock = threading.Lock()
        self._current_step = None
        self._cancelled_steps = set()
        self._stopped = False

        self._thread = threading.Thread(target=self._worker_loop, daemon=True)
        self._thread.start()

    def close(self):
        with self._lock:
            self._stopped = True
            try:
                self.solver.interrupt()
            except Exception:
                pass
        try:
            self.job_q.put(None)
        except Exception:
            pass
        try:
            self._thread.join(timeout=0.05)
        except Exception:
            pass
        try:
            self.solver.delete()
        except Exception:
            pass
        _safe_close_conn(self.result_send)
        _safe_close_conn(self.result_recv)

    def submit(self, step, idx, radius):
        self.job_q.put((step, idx, radius))

    def cancel(self, step):
        with self._lock:
            self._cancelled_steps.add(step)
            if self._current_step == step and hasattr(self.solver, "interrupt"):
                try:
                    self.solver.interrupt()
                except Exception:
                    pass

    def _worker_loop(self):
        while True:
            msg = self.job_q.get()
            if msg is None:
                break

            step, idx, radius = msg

            with self._lock:
                if self._stopped:
                    break
                self._current_step = step

            try:
                if idx not in self.tested:
                    alpha = self.next_var
                    self.next_var += 1
                    self.tested[idx] = alpha

                    step_clauses = self.inst._build_incremental_guarded_clauses(
                        radius=radius,
                        selector_lit=alpha,
                    )
                    for cl in step_clauses:
                        self.solver.add_clause(cl)
                    self.added_clause_count += len(step_clauses)
                else:
                    alpha = self.tested[idx]

                timer = None
                cpu0 = None
                cpu1 = None
                try:
                    if self.time_limit and self.time_limit > 0 and hasattr(self.solver, "interrupt"):
                        timer = threading.Timer(self.time_limit, self.solver.interrupt)
                        timer.start()
                        cpu0 = cpu_self_seconds()
                        try:
                            sat_res = self.solver.solve_limited(
                                assumptions=[alpha],
                                expect_interrupt=True,
                            )
                        except NotImplementedError:
                            sat_res = self.solver.solve(assumptions=[alpha])
                        cpu1 = cpu_self_seconds()
                    else:
                        cpu0 = cpu_self_seconds()
                        sat_res = self.solver.solve(assumptions=[alpha])
                        cpu1 = cpu_self_seconds()
                finally:
                    if timer:
                        timer.cancel()

                cpu_sec = (cpu1 - cpu0) if (cpu0 is not None and cpu1 is not None) else 0.0

                with self._lock:
                    was_cancelled = step in self._cancelled_steps
                    if was_cancelled:
                        self._cancelled_steps.discard(step)

                if sat_res is None:
                    if was_cancelled:
                        self.result_send.send(("cancelled", step, None))
                    else:
                        self.result_send.send(("timeout", step, cpu_sec))
                elif sat_res:
                    model = self.solver.get_model() or []
                    model_set = set(model)
                    centers = sorted([j for j, v in enumerate(self.y_vars) if v in model_set])
                    self.result_send.send(("result", step, (idx, radius, "sat", cpu_sec, None, None, centers)))
                else:
                    self.result_send.send(("result", step, (idx, radius, "unsat", cpu_sec, None, None, None)))

            except Exception:
                self.result_send.send(("error", step, traceback.format_exc()))

            finally:
                if hasattr(self.solver, "clear_interrupt"):
                    try:
                        self.solver.clear_interrupt()
                    except Exception:
                        pass
                with self._lock:
                    self._current_step = None


# -------------------------------
# Persistent MIP worker
# Cancel by terminating worker proc
# -------------------------------

def _mip_worker_loop(inst, mip_backend, job_q, result_send, time_limit):
    while True:
        msg = job_q.get()
        if msg is None:
            break
        if msg.get("cmd") != "solve":
            continue

        step = msg["step"]
        idx = msg["idx"]
        radius = msg["radius"]

        try:
            if mip_backend == "cplex_mip":
                res = solve_radius_cplex(
                    inst=inst,
                    idx=idx,
                    radius=radius,
                    time_limit=time_limit,
                    data=None,
                )
            elif mip_backend == "gurobi_mip":
                res = solve_radius_gurobi(
                    inst=inst,
                    idx=idx,
                    radius=radius,
                    time_limit=time_limit,
                    data=None,
                )
            else:
                raise ValueError(f"Unsupported mip_backend: {mip_backend}")

            result_send.send(("result", step, res))

        except Exception:
            result_send.send(("error", step, traceback.format_exc()))

    _safe_close_conn(result_send)


def _spawn_mip_worker(inst, mip_backend, time_limit):
    job_q = mp.Queue()
    result_recv, result_send = mp.Pipe(duplex=False)
    proc = mp.Process(
        target=_mip_worker_loop,
        args=(inst, mip_backend, job_q, result_send, time_limit),
        daemon=True,
    )
    proc.start()
    _safe_close_conn(result_send)
    return proc, job_q, result_recv


def _restart_mip_worker(inst, mip_backend, time_limit, proc, job_q, result_recv):
    _terminate_process_fast(proc)
    _safe_close_queue_like(job_q)
    _safe_close_conn(result_recv)
    return _spawn_mip_worker(inst, mip_backend, time_limit)


# -------------------------------
# Main
# -------------------------------

def search_min_radius_hybrid_race(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    mip_backend="cplex_mip",
    radii_workers,
    seed_idx=None,
    mgr=None,
    cancel_dict=None,
):
    search_t0 = time.perf_counter()

    if solver_name == mip_backend:
        raise ValueError(
            f"hybrid_race expects a SAT solver in --solvers, not '{mip_backend}'."
        )

    _check_mip_backend_available(mip_backend)

    use_incremental_sat = solver_name not in EXTERNAL_SOLVERS
    radii = inst.radii
    if not radii:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    lo = 0
    hi = len(radii) - 1

    best_sat_idx = None
    best_centers = None
    best_cpu = None
    best_nvars = None
    best_nclauses = None

    step_counter = itertools.count(1)

    mip_proc, mip_job_q, mip_result_recv = _spawn_mip_worker(inst, mip_backend, time_limit)

    sat_proc = None
    sat_job_q = None
    sat_ctl_q = None
    sat_result_recv = None
    sat_worker = None

    try:
        if use_incremental_sat:
            sat_worker = _PersistentIncrementalSATWorker(
                inst=inst,
                encoding=encoding,
                solver_name=solver_name,
                time_limit=time_limit,
            )
            sat_conn = sat_worker.result_recv
        else:
            sat_job_q = mp.Queue()
            sat_ctl_q = mp.Queue()
            sat_result_recv, sat_result_send = mp.Pipe(duplex=False)
            sat_proc = mp.Process(
                target=_external_sat_worker_loop,
                args=(inst, encoding, solver_name, sat_job_q, sat_ctl_q, sat_result_send, time_limit),
                daemon=True,
            )
            sat_proc.start()
            _safe_close_conn(sat_result_send)
            sat_conn = sat_result_recv

        print(
            f"[HYBRID-RACE-INIT] encoding={encoding} solver={solver_name} "
            f"mip_backend={mip_backend} "
            f"sat_mode={'incremental' if use_incremental_sat else 'external_per_radius'} "
            f"p={inst.p} lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={len(radii)}",
            flush=True,
        )

        while lo <= hi:
            if mip_proc is None or (not mip_proc.is_alive()):
                mip_proc, mip_job_q, mip_result_recv = _spawn_mip_worker(inst, mip_backend, time_limit)

            mid = (lo + hi) // 2
            R = radii[mid]
            step = next(step_counter)

            print(
                f"[HYBRID-RACE-STEP] step={step} lo={lo} hi={hi} mid={mid} R={R} "
                f"mip_backend={mip_backend}",
                flush=True,
            )

            if use_incremental_sat:
                sat_worker.submit(step=step, idx=mid, radius=R)
            else:
                sat_job_q.put({"cmd": "solve", "step": step, "idx": mid, "radius": R})

            mip_job_q.put({"cmd": "solve", "step": step, "idx": mid, "radius": R})

            winner_backend = None
            winner_result = None

            while True:
                conns = [mip_result_recv, sat_conn]
                ready = wait(conns, timeout=0.01)
                if not ready:
                    continue

                for conn in ready:
                    msg = _recv_if_matching(conn, step)
                    if msg is None:
                        continue

                    tag, payload = msg

                    if conn is sat_conn:
                        if tag == "result":
                            idx1, R1, s_status, s_cpu, s_nvars, s_nclauses, s_centers = payload
                            if s_status in ("sat", "unsat"):
                                winner_backend = "sat"
                                winner_result = payload

                                mip_proc, mip_job_q, mip_result_recv = _restart_mip_worker(
                                    inst, mip_backend, time_limit, mip_proc, mip_job_q, mip_result_recv
                                )

                                print(
                                    f"[HYBRID-RACE-WIN] step={step} idx={idx1} R={R1} "
                                    f"winner=SAT status={s_status} cpu={s_cpu:.6f}s",
                                    flush=True,
                                )
                                break
                        elif tag == "error":
                            print(f"[HYBRID-RACE-SAT-ERROR]\n{payload}", file=sys.stderr, flush=True)

                    else:
                        if tag == "result":
                            idx2, R2, m_status, m_cpu, m_nvars, m_nclauses, m_centers = payload
                            if m_status in ("sat", "unsat"):
                                winner_backend = mip_backend
                                winner_result = payload

                                if use_incremental_sat:
                                    sat_worker.cancel(step)
                                else:
                                    _safe_put(sat_ctl_q, {"cmd": "cancel", "step": step})

                                print(
                                    f"[HYBRID-RACE-WIN] step={step} idx={idx2} R={R2} "
                                    f"winner={mip_backend} status={m_status} cpu={m_cpu:.6f}s",
                                    flush=True,
                                )
                                break
                        elif tag == "error":
                            print(f"[HYBRID-RACE-MIP-ERROR]\n{payload}", file=sys.stderr, flush=True)

                if winner_backend is not None:
                    break

            if winner_backend == "sat":
                idx1, R1, s_status, s_cpu, s_nvars, s_nclauses, s_centers = winner_result
                if s_status == "sat":
                    best_sat_idx = idx1
                    best_centers = s_centers
                    best_cpu = s_cpu
                    best_nvars = s_nvars if s_nvars is not None else (
                        sat_worker.next_var - 1 if use_incremental_sat else None
                    )
                    best_nclauses = s_nclauses if s_nclauses is not None else (
                        len(sat_worker.base_cnf.clauses) + sat_worker.added_clause_count
                        if use_incremental_sat else None
                    )
                    lo = idx1 + 1
                else:
                    hi = idx1 - 1

            elif winner_backend == mip_backend:
                idx2, R2, m_status, m_cpu, m_nvars, m_nclauses, m_centers = winner_result
                if m_status == "sat":
                    best_sat_idx = idx2
                    best_centers = m_centers
                    best_cpu = m_cpu
                    best_nvars = m_nvars
                    best_nclauses = m_nclauses
                    lo = idx2 + 1
                else:
                    hi = idx2 - 1

            else:
                if best_sat_idx is not None:
                    break
                search_elapsed = time.perf_counter() - search_t0
                return "error", None, None, None, None, None, search_elapsed

    finally:
        if sat_worker is not None:
            sat_worker.close()

        if sat_job_q is not None:
            _safe_put(sat_ctl_q, {"cmd": "stop"})
            _safe_put(sat_ctl_q, None)
            _safe_put(sat_job_q, None)

        if sat_proc is not None:
            sat_proc.join(timeout=0.05)
            _terminate_process_fast(sat_proc)

        _safe_put(mip_job_q, None)
        if mip_proc is not None:
            mip_proc.join(timeout=0.05)
            _terminate_process_fast(mip_proc)

        _safe_close_queue_like(mip_job_q)
        _safe_close_conn(mip_result_recv)
        _safe_close_queue_like(sat_job_q)
        _safe_close_queue_like(sat_ctl_q)
        _safe_close_conn(sat_result_recv)

    if best_sat_idx is None:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    best_radius = radii[best_sat_idx]
    search_elapsed = time.perf_counter() - search_t0
    return (
        "OK",
        best_radius,
        best_nvars,
        best_nclauses,
        best_centers,
        best_cpu,
        search_elapsed,
    )