import time
import threading
import multiprocessing as mp
import traceback
import sys
import queue as pyqueue
import itertools

try:
    import cplex
except Exception:
    cplex = None

try:
    import gurobipy as gp
except Exception:
    gp = None

from pysat.solvers import Solver


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


def _drain_control_queue(control_q, pending_cancel_steps):
    while True:
        try:
            msg = control_q.get_nowait()
        except pyqueue.Empty:
            break
        except Exception:
            break
        if msg is None:
            continue
        cmd = msg.get("cmd")
        if cmd == "cancel":
            pending_cancel_steps.add(msg.get("step"))
        elif cmd == "stop":
            pending_cancel_steps.add("__STOP__")


def _terminate_process_fast(proc, *, grace=0.05):
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


def _wait_for_cancel_ack(result_q, step, timeout=1.0):
    deadline = time.perf_counter() + timeout
    while time.perf_counter() < deadline:
        remaining = max(0.0, deadline - time.perf_counter())
        try:
            tag, step_got, payload = result_q.get(timeout=min(0.02, remaining))
        except pyqueue.Empty:
            continue
        except Exception:
            return None

        if step_got != step:
            continue

        if tag in ("cancelled", "result", "error"):
            return (tag, payload)

    return None


def _check_mip_backend_available(mip_backend: str):
    if mip_backend == "cplex_mip":
        if cplex is None:
            raise ImportError("hybrid_race with mip_backend='cplex_mip' needs CPLEX installed.")
    elif mip_backend == "gurobi_mip":
        if gp is None:
            raise ImportError("hybrid_race with mip_backend='gurobi_mip' needs Gurobi installed.")
    else:
        raise ValueError(f"Unsupported mip_backend: {mip_backend}")


def _mip_solve_child(inst, mip_backend, step, idx, radius, time_limit, out_q):
    try:
        if mip_backend == "cplex_mip":
            from solvers_cplex import _solve_radius_cplex
            res = _solve_radius_cplex(
                inst=inst,
                idx=idx,
                radius=radius,
                time_limit=time_limit,
                data=None,
            )
        elif mip_backend == "gurobi_mip":
            from solvers_gurobi import _solve_radius_gurobi
            res = _solve_radius_gurobi(
                inst=inst,
                idx=idx,
                radius=radius,
                time_limit=time_limit,
                data=None,
            )
        else:
            raise ValueError(f"Unsupported mip_backend: {mip_backend}")

        out_q.put(("result", step, res))
    except Exception:
        out_q.put(("error", step, traceback.format_exc()))


def _mip_worker_loop(inst, mip_backend, job_q, control_q, result_q, time_limit):
    pending_cancel_steps = set()

    while True:
        _drain_control_queue(control_q, pending_cancel_steps)
        if "__STOP__" in pending_cancel_steps:
            break

        msg = job_q.get()
        if msg is None:
            break
        if msg.get("cmd") != "solve":
            continue

        step = msg["step"]
        idx = msg["idx"]
        radius = msg["radius"]

        if step in pending_cancel_steps:
            pending_cancel_steps.discard(step)
            result_q.put(("cancelled", step, None))
            continue

        child_q = mp.Queue()
        child = mp.Process(
            target=_mip_solve_child,
            args=(inst, mip_backend, step, idx, radius, time_limit, child_q),
            daemon=True,
        )
        child.start()

        sent_terminal_msg = False

        while True:
            _drain_control_queue(control_q, pending_cancel_steps)

            if step in pending_cancel_steps:
                pending_cancel_steps.discard(step)
                _terminate_process_fast(child)
                result_q.put(("cancelled", step, None))
                sent_terminal_msg = True
                break

            try:
                item = child_q.get_nowait()
                result_q.put(item)
                sent_terminal_msg = True
                break
            except pyqueue.Empty:
                pass
            except Exception:
                pass

            if not child.is_alive():
                break

            time.sleep(0.0005)

        try:
            child.join(timeout=0.1)
        except Exception:
            pass

        if not sent_terminal_msg:
            try:
                item = child_q.get_nowait()
                result_q.put(item)
                sent_terminal_msg = True
            except pyqueue.Empty:
                pass
            except Exception:
                pass

        if (not sent_terminal_msg) and (not child.is_alive()):
            result_q.put(("cancelled", step, None))

        try:
            child_q.close()
            child_q.join_thread()
        except Exception:
            pass


def _solve_radius_external_once(inst, idx, encoding, solver_name, radius, time_limit, cancel_ev=None):
    cnf, varmap = inst._encode_cnf(radius, encoding)

    if cnf is None:
        return idx, radius, "unsat", 0.0, None, None, None

    from experiments import _run_external_solver
    status, cpu_sec, model = _run_external_solver(
        solver_name=solver_name,
        cnf=cnf,
        time_limit=time_limit,
        cancel_ev=cancel_ev,
    )

    nvars = cnf.nv
    nclauses = len(cnf.clauses)

    if status in ("timeout", "error"):
        return idx, radius, status, cpu_sec, nvars, nclauses, None

    if status == "unsat":
        return idx, radius, "unsat", cpu_sec, nvars, nclauses, None

    model_set = set(model or [])
    y_vars = varmap.get("y", [])
    Nc = set(varmap.get("Nc", []))
    candidates = set(varmap.get("candidates", []))
    chosen = {
        j for j, v in enumerate(y_vars)
        if (v in model_set) and (j in candidates)
    }
    centers = sorted(chosen | Nc)

    return idx, radius, "sat", cpu_sec, nvars, nclauses, centers


def _external_sat_solve_child(inst, encoding, solver_name, step, idx, radius, time_limit, out_q):
    try:
        res = _solve_radius_external_once(
            inst, idx, encoding, solver_name, radius, time_limit, cancel_ev=None
        )
        out_q.put(("result", step, res))
    except Exception:
        out_q.put(("error", step, traceback.format_exc()))


def _external_sat_worker_loop(inst, encoding, solver_name, job_q, control_q, result_q, time_limit):
    pending_cancel_steps = set()

    while True:
        _drain_control_queue(control_q, pending_cancel_steps)
        if "__STOP__" in pending_cancel_steps:
            break

        msg = job_q.get()
        if msg is None:
            break
        if msg.get("cmd") != "solve":
            continue

        step = msg["step"]
        idx = msg["idx"]
        radius = msg["radius"]

        if step in pending_cancel_steps:
            pending_cancel_steps.discard(step)
            result_q.put(("cancelled", step, None))
            continue

        child_q = mp.Queue()
        child = mp.Process(
            target=_external_sat_solve_child,
            args=(inst, encoding, solver_name, step, idx, radius, time_limit, child_q),
            daemon=True,
        )
        child.start()

        sent_terminal_msg = False

        while True:
            _drain_control_queue(control_q, pending_cancel_steps)

            if step in pending_cancel_steps:
                pending_cancel_steps.discard(step)
                _terminate_process_fast(child)
                result_q.put(("cancelled", step, None))
                sent_terminal_msg = True
                break

            try:
                item = child_q.get_nowait()
                result_q.put(item)
                sent_terminal_msg = True
                break
            except pyqueue.Empty:
                pass
            except Exception:
                pass

            if not child.is_alive():
                break

            time.sleep(0.0005)

        try:
            child.join(timeout=0.1)
        except Exception:
            pass

        if not sent_terminal_msg:
            try:
                item = child_q.get_nowait()
                result_q.put(item)
                sent_terminal_msg = True
            except pyqueue.Empty:
                pass
            except Exception:
                pass

        if (not sent_terminal_msg) and (not child.is_alive()):
            result_q.put(("cancelled", step, None))

        try:
            child_q.close()
            child_q.join_thread()
        except Exception:
            pass


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

    def close(self):
        try:
            self.solver.delete()
        except Exception:
            pass

    def solve_radius(self, idx, radius, cancel_ev):
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

        done_ev = threading.Event()
        out_q = pyqueue.Queue()

        def _worker():
            timer = None
            cpu0 = None
            cpu1 = None
            try:
                if cancel_ev is not None and hasattr(self.solver, "interrupt"):
                    def _watch_cancel():
                        cancel_ev.wait()
                        if not done_ev.is_set():
                            try:
                                self.solver.interrupt()
                            except Exception:
                                pass
                    threading.Thread(target=_watch_cancel, daemon=True).start()

                if self.time_limit and self.time_limit > 0 and hasattr(self.solver, "interrupt"):
                    timer = threading.Timer(self.time_limit, self.solver.interrupt)
                    timer.start()
                    from experiments import _cpu_self_seconds as get_cpu_seconds
                    cpu0 = get_cpu_seconds()
                    try:
                        sat_res = self.solver.solve_limited(
                            assumptions=[alpha],
                            expect_interrupt=True
                        )
                    except NotImplementedError:
                        sat_res = self.solver.solve(assumptions=[alpha])
                    cpu1 = get_cpu_seconds()
                else:
                    from experiments import _cpu_self_seconds as get_cpu_seconds
                    cpu0 = get_cpu_seconds()
                    sat_res = self.solver.solve(assumptions=[alpha])
                    cpu1 = get_cpu_seconds()

                cpu_sec = (cpu1 - cpu0) if (cpu0 is not None and cpu1 is not None) else 0.0

                if sat_res is None:
                    if cancel_ev is not None and cancel_ev.is_set():
                        out_q.put(("cancelled", cpu_sec, None))
                    else:
                        out_q.put(("timeout", cpu_sec, None))
                    return

                if sat_res:
                    model = self.solver.get_model() or []
                    model_set = set(model)
                    centers = sorted([j for j, v in enumerate(self.y_vars) if v in model_set])
                    out_q.put(("sat", cpu_sec, centers))
                else:
                    out_q.put(("unsat", cpu_sec, None))

            except Exception:
                out_q.put(("error", 0.0, None))
            finally:
                if timer:
                    timer.cancel()
                done_ev.set()
                if hasattr(self.solver, "clear_interrupt"):
                    try:
                        self.solver.clear_interrupt()
                    except Exception:
                        pass

        th = threading.Thread(target=_worker, daemon=True)
        th.start()
        return th, out_q


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

    from experiments import EXTERNAL_SOLVERS
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

    mip_job_q = mp.Queue()
    mip_ctl_q = mp.Queue()
    mip_result_q = mp.Queue()
    mip_proc = mp.Process(
        target=_mip_worker_loop,
        args=(inst, mip_backend, mip_job_q, mip_ctl_q, mip_result_q, time_limit),
    )
    mip_proc.start()

    sat_proc = None
    sat_job_q = None
    sat_ctl_q = None
    sat_result_q = None
    sat_worker = None

    try:
        if use_incremental_sat:
            sat_worker = _PersistentIncrementalSATWorker(
                inst=inst,
                encoding=encoding,
                solver_name=solver_name,
                time_limit=time_limit,
            )
        else:
            sat_job_q = mp.Queue()
            sat_ctl_q = mp.Queue()
            sat_result_q = mp.Queue()
            sat_proc = mp.Process(
                target=_external_sat_worker_loop,
                args=(inst, encoding, solver_name, sat_job_q, sat_ctl_q, sat_result_q, time_limit),
            )
            sat_proc.start()

        print(
            f"[HYBRID-RACE-INIT] encoding={encoding} solver={solver_name} "
            f"mip_backend={mip_backend} "
            f"sat_mode={'incremental' if use_incremental_sat else 'external_per_radius'} "
            f"p={inst.p} lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={len(radii)}",
            flush=True,
        )

        while lo <= hi:
            mid = (lo + hi) // 2
            R = radii[mid]
            step = next(step_counter)

            print(
                f"[HYBRID-RACE-STEP] step={step} lo={lo} hi={hi} mid={mid} R={R} "
                f"mip_backend={mip_backend}",
                flush=True,
            )

            sat_thread = None
            sat_local_q = None
            sat_cancel_ev = threading.Event() if use_incremental_sat else None

            if use_incremental_sat:
                sat_thread, sat_local_q = sat_worker.solve_radius(
                    idx=mid, radius=R, cancel_ev=sat_cancel_ev
                )
            else:
                sat_job_q.put({"cmd": "solve", "step": step, "idx": mid, "radius": R})

            mip_job_q.put({"cmd": "solve", "step": step, "idx": mid, "radius": R})

            winner_backend = None
            winner_result = None

            while True:
                if use_incremental_sat and sat_local_q is not None:
                    try:
                        sat_status, sat_cpu, sat_centers = sat_local_q.get_nowait()
                        if sat_status in ("sat", "unsat"):
                            winner_backend = "sat"
                            winner_result = (mid, R, sat_status, sat_cpu, None, None, sat_centers)

                            _safe_put(mip_ctl_q, {"cmd": "cancel", "step": step})
                            _wait_for_cancel_ack(mip_result_q, step, timeout=1.0)

                            print(
                                f"[HYBRID-RACE-WIN] step={step} idx={mid} R={R} "
                                f"winner=SAT status={sat_status} cpu={sat_cpu:.6f}s",
                                flush=True,
                            )
                            break
                    except pyqueue.Empty:
                        pass

                if (not use_incremental_sat) and sat_result_q is not None:
                    try:
                        tag, step_got, payload = sat_result_q.get_nowait()
                        if step_got == step:
                            if tag == "result":
                                idx1, R1, s_status, s_cpu, s_nvars, s_nclauses, s_centers = payload
                                if s_status in ("sat", "unsat"):
                                    winner_backend = "sat"
                                    winner_result = payload

                                    _safe_put(mip_ctl_q, {"cmd": "cancel", "step": step})
                                    _wait_for_cancel_ack(mip_result_q, step, timeout=1.0)

                                    print(
                                        f"[HYBRID-RACE-WIN] step={step} idx={idx1} R={R1} "
                                        f"winner=SAT status={s_status} cpu={s_cpu:.6f}s",
                                        flush=True,
                                    )
                                    break
                            elif tag == "error":
                                print(f"[HYBRID-RACE-SAT-ERROR]\n{payload}", file=sys.stderr, flush=True)
                    except pyqueue.Empty:
                        pass
                    except Exception:
                        pass

                try:
                    tag, step_got, payload = mip_result_q.get_nowait()
                    if step_got == step:
                        if tag == "result":
                            idx2, R2, m_status, m_cpu, m_nvars, m_nclauses, m_centers = payload
                            if m_status in ("sat", "unsat"):
                                winner_backend = mip_backend
                                winner_result = payload

                                if use_incremental_sat:
                                    sat_cancel_ev.set()

                                    deadline = time.perf_counter() + 1.0
                                    while time.perf_counter() < deadline:
                                        try:
                                            sat_local_q.get(timeout=0.02)
                                            break
                                        except pyqueue.Empty:
                                            pass
                                        except Exception:
                                            break

                                    try:
                                        sat_thread.join(timeout=0.1)
                                    except Exception:
                                        pass

                                    if sat_thread.is_alive():
                                        raise RuntimeError(
                                            f"Incremental SAT thread did not stop cleanly at step={step}"
                                        )
                                else:
                                    _safe_put(sat_ctl_q, {"cmd": "cancel", "step": step})
                                    _wait_for_cancel_ack(sat_result_q, step, timeout=1.0)

                                print(
                                    f"[HYBRID-RACE-WIN] step={step} idx={idx2} R={R2} "
                                    f"winner={mip_backend} status={m_status} cpu={m_cpu:.6f}s",
                                    flush=True,
                                )
                                break
                        elif tag == "error":
                            print(f"[HYBRID-RACE-MIP-ERROR]\n{payload}", file=sys.stderr, flush=True)
                except pyqueue.Empty:
                    pass
                except Exception:
                    pass

                time.sleep(0.0005)

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
            _safe_put(sat_job_q, None)
        if sat_proc is not None:
            sat_proc.join(timeout=1)

        _safe_put(mip_ctl_q, {"cmd": "stop"})
        _safe_put(mip_job_q, None)
        mip_proc.join(timeout=1)

        _safe_close_queue_like(mip_job_q)
        _safe_close_queue_like(mip_ctl_q)
        _safe_close_queue_like(mip_result_q)
        _safe_close_queue_like(sat_job_q)
        _safe_close_queue_like(sat_ctl_q)
        _safe_close_queue_like(sat_result_q)

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