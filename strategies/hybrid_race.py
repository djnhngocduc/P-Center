import time
import multiprocessing as mp
import threading
import traceback
import sys
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
    solve_radius_cplex,
    solve_radius_gurobi,
)


def _safe_close_conn(conn):
    if conn is None:
        return
    try:
        conn.close()
    except Exception:
        pass


def _send_msg(conn, obj):
    try:
        conn.send(obj)
        return True
    except Exception:
        return False


def _terminate_process_fast(proc, *, grace=0.2):
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


def _solve_radius_sat_once(
    inst,
    idx,
    encoding,
    solver_name,
    radius,
    time_limit,
    cancel_ev=None,
    current_cancel_ref=None,
):
    call_t0 = time.perf_counter()

    def _remaining_time():
        if time_limit and time_limit > 0:
            return max(0.0, float(time_limit) - (time.perf_counter() - call_t0))
        return time_limit

    cnf, varmap = inst._encode_cnf(radius, encoding)

    if cnf is None:
        return idx, radius, "unsat", None, None, None

    nvars = cnf.nv
    nclauses = len(cnf.clauses)
    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        return idx, radius, "timeout", nvars, nclauses, None

    if solver_name in EXTERNAL_SOLVERS:
        base_tmp = "/dev/shm" if os.path.isdir("/dev/shm") else None
        local_tmpdir = tempfile.mkdtemp(prefix="pcsat_race_", dir=base_tmp)

        try:
            if current_cancel_ref is not None and cancel_ev is not None:
                current_cancel_ref["fn"] = cancel_ev.set

            status, model = run_external_solver(
                solver_name=solver_name,
                cnf=cnf,
                time_limit=step_time_limit,
                cancel_ev=cancel_ev,
                tmpdir=local_tmpdir,
            )

            if status in ("timeout", "error", "cancelled"):
                return idx, radius, status, nvars, nclauses, None

            if status == "unsat":
                return idx, radius, "unsat", nvars, nclauses, None

            model_set = set(model or [])
            y_vars = varmap.get("y", [])
            candidates = set(varmap.get("candidates", []))
            chosen = {
                j for j, v in enumerate(y_vars)
                if (v in model_set) and (j in candidates)
            }
            centers = sorted(chosen)
            return idx, radius, "sat", nvars, nclauses, centers

        finally:
            if current_cancel_ref is not None:
                current_cancel_ref["fn"] = None
            try:
                shutil.rmtree(local_tmpdir, ignore_errors=True)
            except Exception:
                pass

    solver = Solver(name=solver_name, bootstrap_with=cnf.clauses)
    try:
        if cancel_ev is not None and hasattr(solver, "interrupt"):
            if current_cancel_ref is not None:
                current_cancel_ref["fn"] = solver.interrupt

            def _watch_cancel():
                cancel_ev.wait()
                try:
                    solver.interrupt()
                except Exception:
                    pass

            threading.Thread(target=_watch_cancel, daemon=True).start()

        timer = None
        try:
            if step_time_limit and step_time_limit > 0 and hasattr(solver, "interrupt"):
                timer = threading.Timer(step_time_limit, solver.interrupt)
                timer.start()
                try:
                    sat = solver.solve_limited(expect_interrupt=True)
                except NotImplementedError:
                    sat = solver.solve()
            else:
                sat = solver.solve()
        finally:
            if timer:
                timer.cancel()

        if sat is None:
            if cancel_ev is not None and cancel_ev.is_set():
                return idx, radius, "cancelled", nvars, nclauses, None
            return idx, radius, "timeout", nvars, nclauses, None

        if not sat:
            return idx, radius, "unsat", nvars, nclauses, None

        model = solver.get_model() or []
        model_set = set(model)
        y_vars = varmap.get("y", [])
        candidates = set(varmap.get("candidates", []))
        chosen = {
            j for j, v in enumerate(y_vars)
            if (v in model_set) and (j in candidates)
        }
        centers = sorted(chosen)
        return idx, radius, "sat", nvars, nclauses, centers

    finally:
        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = None
        try:
            if hasattr(solver, "clear_interrupt"):
                solver.clear_interrupt()
        except Exception:
            pass
        try:
            solver.delete()
        except Exception:
            pass


def _sat_worker_loop(inst, encoding, solver_name, job_recv, ctl_recv, result_send, time_limit):
    current_cancel_ref = {"fn": None}
    stop_flag = {"stop": False}

    def _ctl_loop():
        while True:
            try:
                msg = ctl_recv.recv()
            except Exception:
                return
            if msg is None:
                return

            cmd = msg.get("cmd")
            if cmd == "stop":
                stop_flag["stop"] = True
                fn = current_cancel_ref["fn"]
                if fn is not None:
                    try:
                        fn()
                    except Exception:
                        pass
                return

            if cmd == "cancel":
                fn = current_cancel_ref["fn"]
                if fn is not None:
                    try:
                        fn()
                    except Exception:
                        pass

    threading.Thread(target=_ctl_loop, daemon=True).start()

    while not stop_flag["stop"]:
        try:
            msg = job_recv.recv()
        except EOFError:
            break
        except Exception:
            break

        if msg is None:
            break
        if not isinstance(msg, dict) or msg.get("cmd") != "solve":
            continue

        step = msg["step"]
        idx = msg["idx"]
        radius = msg["radius"]
        step_time_limit = msg.get("time_limit", time_limit)

        cancel_ev = threading.Event()

        try:
            res = _solve_radius_sat_once(
                inst=inst,
                idx=idx,
                encoding=encoding,
                solver_name=solver_name,
                radius=radius,
                time_limit=step_time_limit,
                cancel_ev=cancel_ev,
                current_cancel_ref=current_cancel_ref,
            )
            result_send.send(("result", step, res))
        except Exception:
            result_send.send(("error", step, traceback.format_exc()))
        finally:
            current_cancel_ref["fn"] = None

    _safe_close_conn(job_recv)
    _safe_close_conn(ctl_recv)
    _safe_close_conn(result_send)


def _mip_worker_loop(inst, mip_backend, job_recv, ctl_recv, result_send, time_limit):
    current_cancel_ref = {"fn": None}
    stop_flag = {"stop": False}
    env = None

    def _ctl_loop():
        while True:
            try:
                msg = ctl_recv.recv()
            except Exception:
                return
            if msg is None:
                return

            cmd = msg.get("cmd")
            if cmd == "stop":
                stop_flag["stop"] = True
                fn = current_cancel_ref["fn"]
                if fn is not None:
                    try:
                        fn()
                    except Exception:
                        pass
                return

            if cmd == "cancel":
                fn = current_cancel_ref["fn"]
                if fn is not None:
                    try:
                        fn()
                    except Exception:
                        pass

    threading.Thread(target=_ctl_loop, daemon=True).start()

    try:
        if mip_backend == "gurobi_mip":
            env = gp.Env(empty=True)
            env.start()

        while not stop_flag["stop"]:
            try:
                msg = job_recv.recv()
            except EOFError:
                break
            except Exception:
                break

            if msg is None:
                break
            if not isinstance(msg, dict) or msg.get("cmd") != "solve":
                continue

            step = msg["step"]
            idx = msg["idx"]
            radius = msg["radius"]
            step_time_limit = msg.get("time_limit", time_limit)

            try:
                if mip_backend == "cplex_mip":
                    res = solve_radius_cplex(
                        inst=inst,
                        idx=idx,
                        radius=radius,
                        time_limit=step_time_limit,
                        data=None,
                        current_cancel_ref=current_cancel_ref,
                    )
                elif mip_backend == "gurobi_mip":
                    res = solve_radius_gurobi(
                        inst=inst,
                        idx=idx,
                        radius=radius,
                        time_limit=step_time_limit,
                        data=None,
                        env=env,
                        current_cancel_ref=current_cancel_ref,
                    )
                else:
                    raise ValueError(f"Unsupported mip_backend: {mip_backend}")

                result_send.send(("result", step, res))
            except Exception:
                result_send.send(("error", step, traceback.format_exc()))
            finally:
                current_cancel_ref["fn"] = None

    finally:
        if env is not None:
            try:
                env.dispose()
            except Exception:
                pass
        _safe_close_conn(job_recv)
        _safe_close_conn(ctl_recv)
        _safe_close_conn(result_send)


def _spawn_sat_worker(inst, encoding, solver_name, time_limit):
    job_recv, job_send = mp.Pipe(duplex=False)
    ctl_recv, ctl_send = mp.Pipe(duplex=False)
    result_recv, result_send = mp.Pipe(duplex=False)
    proc = mp.Process(
        target=_sat_worker_loop,
        args=(inst, encoding, solver_name, job_recv, ctl_recv, result_send, time_limit),
        daemon=True,
    )
    proc.start()
    _safe_close_conn(job_recv)
    _safe_close_conn(ctl_recv)
    _safe_close_conn(result_send)
    return proc, job_send, ctl_send, result_recv


def _spawn_mip_worker(inst, mip_backend, time_limit):
    job_recv, job_send = mp.Pipe(duplex=False)
    ctl_recv, ctl_send = mp.Pipe(duplex=False)
    result_recv, result_send = mp.Pipe(duplex=False)
    proc = mp.Process(
        target=_mip_worker_loop,
        args=(inst, mip_backend, job_recv, ctl_recv, result_send, time_limit),
        daemon=True,
    )
    proc.start()
    _safe_close_conn(job_recv)
    _safe_close_conn(ctl_recv)
    _safe_close_conn(result_send)
    return proc, job_send, ctl_send, result_recv


def search_min_radius_hybrid_race(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    mip_backend="cplex_mip",
    seed_idx=None,
    verbose=True,
):
    search_t0 = time.perf_counter()
    deadline = (
        search_t0 + float(time_limit)
        if time_limit and time_limit > 0
        else None
    )

    def _remaining_time():
        if deadline is None:
            return time_limit
        return max(0.0, deadline - time.perf_counter())

    def _log(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    if solver_name == mip_backend:
        raise ValueError(
            f"hybrid_race expects a SAT solver in --solvers, not '{mip_backend}'."
        )

    _check_mip_backend_available(mip_backend)

    radii = inst.radii
    if not radii:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, search_elapsed

    lo = 0
    hi = len(radii) - 1

    best_sat_idx = None
    best_centers = None
    best_nvars = None
    best_nclauses = None
    decided = {}

    def _timeout_result():
        search_elapsed = time.perf_counter() - search_t0
        if best_sat_idx is not None:
            return (
                "timeout_with_incumbent",
                radii[best_sat_idx],
                best_nvars,
                best_nclauses,
                best_centers,
                search_elapsed,
            )
        return "timeout", None, None, None, None, search_elapsed

    step_counter = itertools.count(1)

    sat_proc, sat_job_send, sat_ctl_send, sat_result_recv = _spawn_sat_worker(
        inst, encoding, solver_name, time_limit
    )
    mip_proc, mip_job_send, mip_ctl_send, mip_result_recv = _spawn_mip_worker(
        inst, mip_backend, time_limit
    )

    try:
        _log(
            f"[HYBRID-RACE-INIT] encoding={encoding} solver={solver_name} "
            f"mip_backend={mip_backend} "
            f"sat_mode=worker+control mip_mode=worker+control "
            f"p={inst.p} lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={len(radii)}",
            flush=True,
        )

        while lo <= hi:
            step_time_limit = _remaining_time()
            if deadline is not None and step_time_limit <= 0:
                _log("[HYBRID-RACE-TIMEOUT] instance-level time budget exhausted", flush=True)
                return _timeout_result()

            mid = (lo + hi) // 2
            R = radii[mid]
            step = next(step_counter)

            _log(
                f"[HYBRID-RACE-STEP] step={step} lo={lo} hi={hi} mid={mid} R={R} "
                f"mip_backend={mip_backend}",
                flush=True,
            )

            job = {
                "cmd": "solve",
                "step": step,
                "idx": mid,
                "radius": R,
                "time_limit": step_time_limit,
            }
            _send_msg(sat_job_send, job)
            _send_msg(mip_job_send, job)

            winner_backend = None
            winner_result = None

            sat_seen = None
            mip_seen = None

            while True:
                wait_timeout = 0.01
                if deadline is not None:
                    remaining = _remaining_time()
                    if remaining <= 0:
                        _send_msg(sat_ctl_send, {"cmd": "cancel", "step": step})
                        _send_msg(mip_ctl_send, {"cmd": "cancel", "step": step})
                        _log(
                            "[HYBRID-RACE-TIMEOUT] instance-level time budget exhausted",
                            flush=True,
                        )
                        return _timeout_result()
                    wait_timeout = min(wait_timeout, remaining)

                ready = wait([sat_result_recv, mip_result_recv], timeout=wait_timeout)
                if not ready:
                    continue

                for conn in ready:
                    msg = _recv_if_matching(conn, step)
                    if msg is None:
                        continue

                    tag, payload = msg

                    if conn is sat_result_recv:
                        if tag == "result":
                            sat_seen = payload
                            idx1, R1, s_status, s_nvars, s_nclauses, s_centers = payload
                            if s_status in ("sat", "unsat"):
                                winner_backend = "sat"
                                winner_result = payload
                                _send_msg(mip_ctl_send, {"cmd": "cancel", "step": step})
                                _log(
                                    f"[HYBRID-RACE-WIN] step={step} idx={idx1} R={R1} "
                                    f"winner=SAT status={s_status}",
                                    flush=True,
                                )
                                break
                        elif tag == "error":
                            sat_seen = ("error", payload)
                            print(f"[HYBRID-RACE-SAT-ERROR]\n{payload}", file=sys.stderr, flush=True)

                    else:
                        if tag == "result":
                            mip_seen = payload
                            idx2, R2, m_status, m_nvars, m_nclauses, m_centers = payload
                            if m_status in ("sat", "unsat"):
                                winner_backend = mip_backend
                                winner_result = payload
                                _send_msg(sat_ctl_send, {"cmd": "cancel", "step": step})
                                _log(
                                    f"[HYBRID-RACE-WIN] step={step} idx={idx2} R={R2} "
                                    f"winner={mip_backend} status={m_status}",
                                    flush=True,
                                )
                                break
                        elif tag == "error":
                            mip_seen = ("error", payload)
                            print(f"[HYBRID-RACE-MIP-ERROR]\n{payload}", file=sys.stderr, flush=True)

                if winner_backend is not None:
                    break

                sat_done = sat_seen is not None
                mip_done = mip_seen is not None
                if sat_done and mip_done:
                    break

            if winner_backend == "sat":
                idx1, R1, s_status, s_nvars, s_nclauses, s_centers = winner_result
                decided[idx1] = s_status

                if s_status == "sat":
                    best_sat_idx = idx1
                    best_centers = s_centers
                    best_nvars = s_nvars
                    best_nclauses = s_nclauses
                    lo = idx1 + 1
                else:
                    hi = idx1 - 1

            elif winner_backend == mip_backend:
                idx2, R2, m_status, m_nvars, m_nclauses, m_centers = winner_result
                decided[idx2] = m_status

                if m_status == "sat":
                    best_sat_idx = idx2
                    best_centers = m_centers
                    best_nvars = m_nvars
                    best_nclauses = m_nclauses
                    lo = idx2 + 1
                else:
                    hi = idx2 - 1

            else:
                return _timeout_result()

    finally:
        _send_msg(sat_ctl_send, {"cmd": "stop"})
        _send_msg(sat_job_send, None)

        _send_msg(mip_ctl_send, {"cmd": "stop"})
        _send_msg(mip_job_send, None)

        if sat_proc is not None:
            sat_proc.join(timeout=0.2)
            _terminate_process_fast(sat_proc, grace=0.2)

        if mip_proc is not None:
            mip_proc.join(timeout=0.2)
            _terminate_process_fast(mip_proc, grace=0.2)

        _safe_close_conn(sat_job_send)
        _safe_close_conn(sat_ctl_send)
        _safe_close_conn(sat_result_recv)

        _safe_close_conn(mip_job_send)
        _safe_close_conn(mip_ctl_send)
        _safe_close_conn(mip_result_recv)

    if best_sat_idx is None:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, search_elapsed

    best_radius = radii[best_sat_idx]
    cert_unsat_idx = best_sat_idx + 1
    nR = len(radii)
    is_smallest_candidate = best_sat_idx == nR - 1
    certified = is_smallest_candidate or (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        final_status = "OK"
        if is_smallest_candidate:
            _log(
                f"[HYBRID-RACE-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={best_radius} "
                f"is the smallest candidate radius",
                flush=True,
            )
        else:
            _log(
                f"[HYBRID-RACE-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={best_radius}, "
                f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={radii[cert_unsat_idx]}",
                flush=True,
            )
    else:
        final_status = "uncertified"
        _log(
            f"[HYBRID-RACE-FALLBACK] best SAT found but UNSAT boundary not certified; "
            f"best_sat_idx={best_sat_idx}, best_R={best_radius}",
            flush=True,
        )

    search_elapsed = time.perf_counter() - search_t0
    return (
        final_status,
        best_radius,
        best_nvars,
        best_nclauses,
        best_centers,
        search_elapsed,
    )
