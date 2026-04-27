import os
import sys
import time
import tempfile
import subprocess
import traceback
import multiprocessing as mp
import threading
from typing import Optional, Tuple, List
import resource

from pysat.formula import WCNF
from pysat.examples.rc2 import RC2

from solvers.backend import cpu_self_seconds

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PYSAT_MAXSAT_SOLVERS = {
    "rc2",
}

EXTERNAL_MAXSAT_SOLVERS = {
    "openwbo": (
        os.path.join(BASE_DIR, "solvers", "open-wbo", "open-wbo_release"),
        []
    ),
}


def _cpu_children_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return float(usage.ru_utime + usage.ru_stime)


def _extract_time_wrapper_cpu(stderr_text: str):
    if not stderr_text:
        return None

    user_t = None
    sys_t = None

    for line in stderr_text.splitlines():
        s = line.strip()
        if s.startswith("user "):
            try:
                user_t = float(s.split()[1])
            except Exception:
                pass
        elif s.startswith("sys "):
            try:
                sys_t = float(s.split()[1])
            except Exception:
                pass

    if user_t is not None and sys_t is not None:
        return user_t + sys_t

    return None

def _run_external_maxsat_solver(
    solver_name,
    wcnf,
    time_limit,
    cancel_ev=None,
    tmpdir=None,
):
    bin_path, extra_args = EXTERNAL_MAXSAT_SOLVERS[solver_name]

    if (not os.path.exists(bin_path)) or (not os.access(bin_path, os.X_OK)):
        print(
            f"[MAXSAT-SOLVER-ERROR] solver={solver_name}: "
            f"binary not found or not executable: {bin_path}",
            file=sys.stderr,
            flush=True,
        )
        return "error", None, [], 0.0

    base_tmp = tmpdir if tmpdir is not None else (
        "/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir()
    )

    wcnf_path = os.path.join(
        base_tmp,
        f"pcmaxsat_{os.getpid()}_{threading.get_ident()}.wcnf",
    )

    _write_wcnf(wcnf, wcnf_path)

    solver_dir = os.path.dirname(bin_path)

    env = os.environ.copy()
    if tmpdir is not None:
        env["TMPDIR"] = tmpdir

    time_bin = "/usr/bin/time"
    use_time_wrapper = os.path.exists(time_bin) and os.access(time_bin, os.X_OK)

    solver_args = extra_args + [wcnf_path]

    if use_time_wrapper:
        cmd = [time_bin, "-p", bin_path] + solver_args
    else:
        cmd = [bin_path] + solver_args

    child_cpu0 = _cpu_children_seconds()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="ignore",
        cwd=solver_dir,
        env=env,
    )

    if cancel_ev is not None:
        def _watch_cancel():
            cancel_ev.wait()
            if proc.poll() is None:
                try:
                    proc.kill()
                except Exception:
                    pass

        threading.Thread(target=_watch_cancel, daemon=True).start()

    try:
        stdout, stderr = proc.communicate(
            timeout=time_limit if (time_limit and time_limit > 0) else None
        )
    except subprocess.TimeoutExpired:
        try:
            proc.kill()
        except Exception:
            pass

        try:
            stdout, stderr = proc.communicate(timeout=1)
        except Exception:
            stdout, stderr = "", ""

        print(f"[MAXSAT-SOLVER-TIMEOUT] solver={solver_name} cmd={cmd}", flush=True)
        cpu_time = _extract_time_wrapper_cpu(stderr or "")
        if cpu_time is None:
            cpu_time = max(0.0, _cpu_children_seconds() - child_cpu0)
        return "timeout", None, [], cpu_time

    finally:
        try:
            if os.path.exists(wcnf_path):
                os.remove(wcnf_path)
        except Exception:
            pass

    status, cost, model = _parse_maxsat_model(stdout or "")

    cpu_time = _extract_time_wrapper_cpu(stderr or "") if use_time_wrapper else None

    if cpu_time is None:
        cpu_time = max(0.0, _cpu_children_seconds() - child_cpu0)

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
        return "error", cost, [], cpu_time

    return status, cost, model, cpu_time

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
            elif "satisfiable" in low:
                status = "sat"
            elif "unsatisfiable" in low:
                status = "unsat"

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

    if status == "sat":
        status = "optimum"

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
                cpu0 = cpu_self_seconds()
                model = solver.compute()
                cpu1 = cpu_self_seconds()
                cost = solver.cost
        else:
            raise ValueError(f"Unsupported PySAT MaxSAT solver: {solver_name}")

        cpu = (cpu1 - cpu0) if (cpu0 is not None and cpu1 is not None) else 0.0

        if model is None:
            send_conn.send(("error", cost, [], cpu))
        else:
            send_conn.send(("optimum", cost, model, cpu))

    except Exception:
        send_conn.send(("error", None, [], 0.0, traceback.format_exc()))

    finally:
        try:
            send_conn.close()
        except Exception:
            pass


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
    recv_conn, send_conn = mp.Pipe(duplex=False)
    proc = mp.Process(
        target=_pysat_maxsat_worker,
        args=(_wcnf_to_data(wcnf), solver_name, send_conn),
        daemon=True,
    )

    child_cpu0 = _cpu_children_seconds()
    proc.start()

    try:
        send_conn.close()
    except Exception:
        pass

    timeout = time_limit if (time_limit and time_limit > 0) else None

    try:
        if timeout is None:
            msg = recv_conn.recv()
        else:
            if not recv_conn.poll(timeout):
                _terminate_process_fast(proc)
                cpu = max(0.0, _cpu_children_seconds() - child_cpu0)
                return "timeout", None, [], cpu

            msg = recv_conn.recv()

        proc.join(timeout=0.2)

        if not msg:
            cpu = max(0.0, _cpu_children_seconds() - child_cpu0)
            return "error", None, [], cpu

        status, cost, model, cpu, *rest = msg

        if status == "error" and rest:
            print(rest[0], file=sys.stderr, flush=True)

        return status, cost, model, cpu

    except EOFError:
        cpu = max(0.0, _cpu_children_seconds() - child_cpu0)
        return "error", None, [], cpu

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
    radii_workers,
    seed_idx=None,
    mgr=None,
    cancel_dict=None,
):
    if encoding not in ("maxsat_setcover", "maxsat_cover"):
        raise ValueError(
            "MaxSAT set-cover mode expects --encodings maxsat_setcover "
            "(or alias maxsat_cover)."
        )

    search_t0 = time.perf_counter()

    radii = inst.radii
    nR = len(radii)

    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    lo = 0
    hi = nR - 1

    best_sat_idx = None
    best_radius = None
    best_centers = None
    best_nvars = None
    best_nclauses = None
    best_cpu = None

    last_nvars = None
    last_nclauses = None
    last_cpu = None

    decided = {}

    while lo <= hi:
        mid = (lo + hi) // 2
        radius = radii[mid]

        if cancel_dict is not None and cancel_dict.get("cancel", False):
            search_elapsed = time.perf_counter() - search_t0
            if best_sat_idx is not None:
                return (
                    "cancelled_with_incumbent",
                    best_radius,
                    best_nvars,
                    best_nclauses,
                    best_centers,
                    best_cpu,
                    search_elapsed,
                )
            return (
                "cancelled",
                None,
                last_nvars,
                last_nclauses,
                None,
                last_cpu,
                search_elapsed,
            )

        wcnf, info = inst._encode_wcnf_maxsat_setcover(radius)

        if wcnf is None:
            decided[mid] = "unsat"
            last_cpu = 0.0
            hi = mid - 1
            continue

        status, cost, model, cpu_sec = _solve_wcnf(
            wcnf,
            solver_name,
            time_limit,
        )

        last_nvars = wcnf.nv
        last_nclauses = len(wcnf.hard) + len(wcnf.soft)
        last_cpu = cpu_sec

        if status not in ("optimum", "sat"):
            decided[mid] = status

            if status in ("timeout", "error") and best_sat_idx is not None:
                search_elapsed = time.perf_counter() - search_t0
                return (
                    f"{status}_with_incumbent",
                    best_radius,
                    best_nvars,
                    best_nclauses,
                    best_centers,
                    best_cpu,
                    search_elapsed,
                )

            search_elapsed = time.perf_counter() - search_t0
            return (
                status,
                None if best_sat_idx is None else best_radius,
                last_nvars,
                last_nclauses,
                None if best_sat_idx is None else best_centers,
                last_cpu,
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
            best_cpu = last_cpu

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
            last_cpu,
            search_elapsed,
        )

    cert_unsat_idx = best_sat_idx + 1
    certified = (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        final_status = "OK"
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
        best_cpu,
        search_elapsed,
    )