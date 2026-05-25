import os
import sys
import time
import tempfile
import subprocess
import traceback
import multiprocessing as mp
import threading
import signal
from typing import Optional, Tuple, List

from pysat.formula import WCNF
from pysat.examples.rc2 import RC2

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PYSAT_MAXSAT_SOLVERS = {
    "rc2",
}

EXTERNAL_MAXSAT_SOLVERS = {
    "openwbo": (
        os.path.join(BASE_DIR, "solvers", "Open-WBO", "open-wbo_release"),
        []
    ),
}


def _kill_subprocess_tree(proc):
    if proc is None:
        return

    try:
        if proc.poll() is not None:
            return
    except Exception:
        return

    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except Exception:
        try:
            proc.kill()
        except Exception:
            pass


def _run_external_maxsat_solver(
    solver_name,
    wcnf,
    time_limit,
    cancel_ev=None,
    tmpdir=None,
):
    call_t0 = time.perf_counter()

    def _remaining_time():
        if time_limit and time_limit > 0:
            return max(0.0, float(time_limit) - (time.perf_counter() - call_t0))
        return time_limit

    bin_path, extra_args = EXTERNAL_MAXSAT_SOLVERS[solver_name]

    if (not os.path.exists(bin_path)) or (not os.access(bin_path, os.X_OK)):
        print(
            f"[MAXSAT-SOLVER-ERROR] solver={solver_name}: "
            f"binary not found or not executable: {bin_path}",
            file=sys.stderr,
            flush=True,
        )
        return "error", None, []

    base_tmp = tmpdir if tmpdir is not None else (
        "/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir()
    )

    wcnf_path = os.path.join(
        base_tmp,
        f"pcmaxsat_{os.getpid()}_{threading.get_ident()}.wcnf",
    )

    _write_wcnf(wcnf, wcnf_path)
    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        try:
            if os.path.exists(wcnf_path):
                os.remove(wcnf_path)
        except Exception:
            pass
        return "timeout", None, []

    solver_dir = os.path.dirname(bin_path)

    env = os.environ.copy()
    if tmpdir is not None:
        env["TMPDIR"] = tmpdir

    solver_args = extra_args + [wcnf_path]
    cmd = [bin_path] + solver_args

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="ignore",
            cwd=solver_dir,
            env=env,
            start_new_session=True,
        )
    except OSError as exc:
        try:
            if os.path.exists(wcnf_path):
                os.remove(wcnf_path)
        except Exception:
            pass
        print(
            f"[MAXSAT-SOLVER-ERROR] solver={solver_name}: could not start "
            f"{bin_path}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        return "error", None, []

    if cancel_ev is not None:
        def _watch_cancel():
            cancel_ev.wait()
            _kill_subprocess_tree(proc)

        threading.Thread(target=_watch_cancel, daemon=True).start()

    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        _kill_subprocess_tree(proc)
        try:
            if os.path.exists(wcnf_path):
                os.remove(wcnf_path)
        except Exception:
            pass
        return "timeout", None, []

    try:
        stdout, stderr = proc.communicate(
            timeout=step_time_limit if (step_time_limit and step_time_limit > 0) else None
        )
    except subprocess.TimeoutExpired:
        _kill_subprocess_tree(proc)

        try:
            stdout, stderr = proc.communicate(timeout=1)
        except Exception:
            stdout, stderr = "", ""

        print(f"[MAXSAT-SOLVER-TIMEOUT] solver={solver_name} cmd={cmd}", flush=True)
        return "timeout", None, []

    finally:
        try:
            if os.path.exists(wcnf_path):
                os.remove(wcnf_path)
        except Exception:
            pass

    status, cost, model = _parse_maxsat_model(stdout or "")

    if status == "error" and proc.returncode not in (0, 10, 20, 30):
        print(
            f"[MAXSAT-SOLVER-ERROR] solver={solver_name} returncode={proc.returncode}",
            file=sys.stderr,
            flush=True,
        )
        if stderr:
            print(stderr, file=sys.stderr, flush=True)

    if status in ("optimum", "sat") and not model:
        print(
            f"[MAXSAT-SOLVER-WARN] {solver_name} reported {status} but no model.",
            flush=True,
        )
        return "error", cost, []

    if status == "sat":
        print(
            f"[MAXSAT-SOLVER-WARN] {solver_name} reported SATISFIABLE "
            "without an optimality proof; treating result as incomplete.",
            flush=True,
        )
        return "incomplete", cost, model

    return status, cost, model


def _write_wcnf(wcnf, path: str) -> None:
    hard = list(getattr(wcnf, "hard", []))
    soft = list(getattr(wcnf, "soft", []))
    weights = list(getattr(wcnf, "wght", []))

    top = int(sum(weights) + 1)
    nclauses = len(hard) + len(soft)

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"p wcnf {wcnf.nv} {nclauses} {top}\n")

        for cl in hard:
            f.write(f"{top} " + " ".join(map(str, cl)) + " 0\n")

        for cl, w in zip(soft, weights):
            f.write(f"{int(w)} " + " ".join(map(str, cl)) + " 0\n")


def _parse_maxsat_model(stdout: str) -> Tuple[str, Optional[int], List[int]]:
    status = "error"
    cost = None
    model = []

    for raw in stdout.splitlines():
        line = raw.strip()
        if not line:
            continue

        low = line.lower()

        if line.startswith("s "):
            if "optimum" in low:
                status = "optimum"
            elif "unsatisfiable" in low:
                status = "unsat"
            elif "satisfiable" in low:
                status = "sat"

        elif line.startswith("o "):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    cost = int(float(parts[1]))
                except Exception:
                    pass

        elif line.startswith("v ") or line.startswith("V "):
            for tok in line.split()[1:]:
                if tok == "0":
                    continue
                try:
                    model.append(int(tok))
                except ValueError:
                    pass

    return status, cost, model


def _wcnf_to_data(wcnf):
    return {
        "nv": int(getattr(wcnf, "nv", 0)),
        "hard": [list(map(int, cl)) for cl in getattr(wcnf, "hard", [])],
        "soft": [list(map(int, cl)) for cl in getattr(wcnf, "soft", [])],
        "weights": [int(w) for w in getattr(wcnf, "wght", [])],
    }


def _wcnf_from_data(data):
    wcnf = WCNF()

    for cl in data["hard"]:
        wcnf.append(cl)

    for cl, weight in zip(data["soft"], data["weights"]):
        wcnf.append(cl, weight=weight)

    wcnf.nv = int(data["nv"])
    return wcnf


def _pysat_maxsat_worker(wcnf_data, solver_name, send_conn):
    try:
        wcnf = _wcnf_from_data(wcnf_data)

        if solver_name in PYSAT_MAXSAT_SOLVERS:
            with RC2(wcnf) as solver:
                model = solver.compute()
                cost = solver.cost
        else:
            raise ValueError(f"Unsupported PySAT MaxSAT solver: {solver_name}")

        if model is None:
            send_conn.send(("error", cost, []))
        else:
            send_conn.send(("optimum", cost, model))

    except Exception:
        send_conn.send(("error", None, [], traceback.format_exc()))

    finally:
        try:
            send_conn.close()
        except Exception:
            pass


def _solve_wcnf_pysat_direct(wcnf, solver_name):
    try:
        if solver_name in PYSAT_MAXSAT_SOLVERS:
            with RC2(wcnf) as solver:
                model = solver.compute()
                cost = solver.cost
        else:
            raise ValueError(f"Unsupported PySAT MaxSAT solver: {solver_name}")

        if model is None:
            return "error", cost, []
        return "optimum", cost, model

    except Exception:
        print(traceback.format_exc(), file=sys.stderr, flush=True)
        return "error", None, []


def _terminate_process_fast(proc, grace=0.2):
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


def _solve_wcnf_pysat(wcnf, solver_name, time_limit):
    call_t0 = time.perf_counter()

    def _remaining_time():
        if time_limit and time_limit > 0:
            return max(0.0, float(time_limit) - (time.perf_counter() - call_t0))
        return time_limit

    try:
        recv_conn, send_conn = mp.Pipe(duplex=False)
    except OSError as exc:
        print(
            f"[MAXSAT-WARN] multiprocessing Pipe unavailable ({exc}); "
            "running RC2 in-process without hard timeout enforcement.",
            file=sys.stderr,
            flush=True,
        )
        return _solve_wcnf_pysat_direct(wcnf, solver_name)

    proc = mp.Process(
        target=_pysat_maxsat_worker,
        args=(_wcnf_to_data(wcnf), solver_name, send_conn),
        daemon=True,
    )

    try:
        proc.start()
    except OSError as exc:
        print(
            f"[MAXSAT-WARN] multiprocessing process unavailable ({exc}); "
            "running RC2 in-process without hard timeout enforcement.",
            file=sys.stderr,
            flush=True,
        )
        try:
            recv_conn.close()
        except Exception:
            pass
        try:
            send_conn.close()
        except Exception:
            pass
        return _solve_wcnf_pysat_direct(wcnf, solver_name)

    try:
        send_conn.close()
    except Exception:
        pass

    timeout = _remaining_time() if (time_limit and time_limit > 0) else None
    if timeout is not None and timeout <= 0:
        _terminate_process_fast(proc)
        return "timeout", None, []

    try:
        if timeout is None:
            msg = recv_conn.recv()
        else:
            if not recv_conn.poll(timeout):
                _terminate_process_fast(proc)
                return "timeout", None, []

            msg = recv_conn.recv()

        proc.join(timeout=0.2)

        if not msg:
            return "error", None, []

        status, cost, model, *rest = msg

        if status == "error" and rest:
            print(rest[0], file=sys.stderr, flush=True)

        return status, cost, model

    except EOFError:
        return "error", None, []

    finally:
        _terminate_process_fast(proc)
        try:
            recv_conn.close()
        except Exception:
            pass


def _solve_wcnf(wcnf, solver_name, time_limit):
    if solver_name in PYSAT_MAXSAT_SOLVERS:
        return _solve_wcnf_pysat(wcnf, solver_name, time_limit)

    if solver_name in EXTERNAL_MAXSAT_SOLVERS:
        return _run_external_maxsat_solver(
            solver_name=solver_name,
            wcnf=wcnf,
            time_limit=time_limit,
        )

    supported = sorted(PYSAT_MAXSAT_SOLVERS) + sorted(EXTERNAL_MAXSAT_SOLVERS)
    raise ValueError(
        f"MaxSAT mode supports only solver_name in {supported}."
    )


def _decode_centers(model, y_vars):
    model_set = set(model or [])
    return sorted([j for j, v in enumerate(y_vars) if v in model_set])


def search_min_radius_maxsat(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    seed_idx=None,
):
    if encoding is not None:
        raise ValueError(
            "MaxSAT mode uses encoding=None and does not need a SAT encoding."
        )

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

    radii = inst.radii
    nR = len(radii)

    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, search_elapsed

    lo = 0
    hi = nR - 1

    best_sat_idx = None
    best_radius = None
    best_centers = None
    best_nvars = None
    best_nclauses = None

    last_nvars = None
    last_nclauses = None

    decided = {}

    while lo <= hi:
        step_time_limit = _remaining_time()
        if deadline is not None and step_time_limit <= 0:
            search_elapsed = time.perf_counter() - search_t0
            if best_sat_idx is not None:
                return (
                    "timeout_with_incumbent",
                    best_radius,
                    best_nvars,
                    best_nclauses,
                    best_centers,
                    search_elapsed,
                )
            return "timeout", None, last_nvars, last_nclauses, None, search_elapsed

        mid = (lo + hi) // 2
        radius = radii[mid]

        wcnf, info = inst._encode_wcnf_maxsat(radius)

        if wcnf is None:
            decided[mid] = "unsat"
            hi = mid - 1
            print(
                f"[MAXSAT-DONE] solver={solver_name} idx={mid} R={radius} "
                f"status=unsat_by_no_coverage cost=None wall=0.000000s",
                flush=True,
            )
            continue

        last_nvars = wcnf.nv
        last_nclauses = len(wcnf.hard) + len(wcnf.soft)
        step_time_limit = _remaining_time()
        if deadline is not None and step_time_limit <= 0:
            search_elapsed = time.perf_counter() - search_t0
            if best_sat_idx is not None:
                return (
                    "timeout_with_incumbent",
                    best_radius,
                    best_nvars,
                    best_nclauses,
                    best_centers,
                    search_elapsed,
                )
            return "timeout", None, last_nvars, last_nclauses, None, search_elapsed

        print(
            f"[MAXSAT-STEP] solver={solver_name} idx={mid} R={radius} "
            f"lo={lo} hi={hi} nvars={last_nvars} "
            f"nclauses={last_nclauses} time_limit={step_time_limit}",
            flush=True,
        )

        step_t0 = time.perf_counter()

        status, cost, model = _solve_wcnf(
            wcnf,
            solver_name,
            step_time_limit,
        )

        step_wall = time.perf_counter() - step_t0

        print(
            f"[MAXSAT-DONE] solver={solver_name} idx={mid} R={radius} "
            f"status={status} cost={cost} wall={step_wall:.6f}s",
            flush=True,
        )

        if status == "unsat":
            decided[mid] = "unsat"
            hi = mid - 1
            continue

        if status != "optimum":
            decided[mid] = status

            if status in ("timeout", "error", "incomplete") and best_sat_idx is not None:
                search_elapsed = time.perf_counter() - search_t0
                return (
                    f"{status}_with_incumbent",
                    best_radius,
                    best_nvars,
                    best_nclauses,
                    best_centers,
                    search_elapsed,
                )

            search_elapsed = time.perf_counter() - search_t0
            return (
                status,
                None if best_sat_idx is None else best_radius,
                last_nvars,
                last_nclauses,
                None if best_sat_idx is None else best_centers,
                search_elapsed,
            )

        centers = _decode_centers(model, info["y"])

        if cost is None:
            cost = len(centers)

        if cost <= inst.p:
            decided[mid] = "sat"

            best_sat_idx = mid
            best_radius = radius
            best_centers = centers
            best_nvars = last_nvars
            best_nclauses = last_nclauses

            lo = mid + 1
        else:
            decided[mid] = "unsat"
            hi = mid - 1

    if best_sat_idx is None:
        search_elapsed = time.perf_counter() - search_t0
        return (
            "infeasible",
            None,
            last_nvars,
            last_nclauses,
            None,
            search_elapsed,
        )

    cert_unsat_idx = best_sat_idx + 1
    is_smallest_candidate = best_sat_idx == nR - 1
    certified = is_smallest_candidate or (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        final_status = "OK"
        if is_smallest_candidate:
            print(
                f"[MAXSAT-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={best_radius} "
                f"is the smallest candidate radius",
                flush=True,
            )
        else:
            print(
                f"[MAXSAT-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={best_radius}, "
                f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={radii[cert_unsat_idx]}",
                flush=True,
            )
    else:
        final_status = "uncertified"
        print(
            f"[MAXSAT-FALLBACK] best SAT found but UNSAT boundary not certified; "
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
