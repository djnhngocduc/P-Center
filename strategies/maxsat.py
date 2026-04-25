import os
import sys
import time
import tempfile
import subprocess
import traceback
import multiprocessing as mp
import threading
from typing import Optional, Tuple, List

from pysat.formula import WCNF
from pysat.examples.rc2 import RC2

from solvers.backend import cpu_self_seconds

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PYSAT_MAXSAT_SOLVERS = {
    "rc2",
}

EXTERNAL_MAXSAT_SOLVERS = {
    "maxcdcl": (
        os.path.join(BASE_DIR, "solvers", "MaxCDCL", "core", "maxcdcl_release"),
        []
    ),
    "openwbo": (
        os.path.join(BASE_DIR, "solvers", "Open-WBO", "open-wbo_release"),
        []
    ),
}

def _run_external_maxsat_solver(
    solver_name,
    wcnf,
    time_limit,
    p_bound=None,
    cancel_ev=None,
    tmpdir=None,
):
    """
    Run an external MaxSAT solver on a temporary WCNF file.

    Supported solvers are defined in EXTERNAL_MAXSAT_SOLVERS.
    This follows the same style as run_external_solver() for external SAT solvers.
    """
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

    if solver_name == "maxcdcl":
        if p_bound is None:
            raise ValueError("maxcdcl requires p_bound for -p-centers.")

        solver_args = [
            f"-filename={wcnf_path}",
            f"-p-centers={int(p_bound)}",
            "-verb=0",
            "-mem-lim=10000",
        ]
    else:
        solver_args = extra_args + [wcnf_path]

    if use_time_wrapper:
        cmd = [time_bin, "-p", bin_path] + solver_args
    else:
        cmd = [bin_path] + solver_args

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
            proc.communicate(timeout=1)
        except Exception:
            pass

        print(f"[MAXSAT-SOLVER-TIMEOUT] solver={solver_name} cmd={cmd}", flush=True)
        return "timeout", None, [], float(time_limit if time_limit else 0.0)

    finally:
        try:
            if os.path.exists(wcnf_path):
                os.remove(wcnf_path)
        except Exception:
            pass

    if solver_name == "maxcdcl":
        status, cost, model = _parse_pcenter_maxcdcl_output(
            stdout or "",
            p_bound=int(p_bound),
        )
    else:
        status, cost, model = _parse_maxsat_model(stdout or "")

    cpu_time = None
    if use_time_wrapper:
        user_t = None
        sys_t = None

        for line in (stderr or "").splitlines():
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
            cpu_time = user_t + sys_t

    if cpu_time is None:
        cpu_time = 0.0

    if status == "error" and proc.returncode not in (0, 10, 20, 30):
        print(
            f"[MAXSAT-SOLVER-ERROR] solver={solver_name} returncode={proc.returncode}",
            file=sys.stderr,
            flush=True,
        )
        if stderr:
            print(stderr, file=sys.stderr, flush=True)

    if status in ("optimum", "sat") and not model:
        if not (
            solver_name == "maxcdcl"
            and cost is not None
            and p_bound is not None
            and cost > int(p_bound)
        ):
            print(
                f"[MAXSAT-SOLVER-WARN] {solver_name} reported {status} but no model.",
                flush=True,
            )
            return "error", cost, [], cpu_time

    return status, cost, model, cpu_time

def _write_wcnf(wcnf, path: str) -> None:
    """
    Write a PySAT WCNF object to weighted DIMACS format.
    This is used by the external MaxCDCL backend.
    """
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
    """
    Parse common MaxSAT output:
      s OPTIMUM FOUND
      o <cost>
      v <lits...>
    """
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

def _parse_pcenter_maxcdcl_output(stdout: str, p_bound: int) -> Tuple[str, Optional[int], List[int]]:
    """
    Parse the custom p-center MaxCDCL output.

    This MaxCDCL build is not a standard MaxSAT CLI solver. It is called with:
      -filename=<wcnf_file> -p-centers=<p> -verb=<level> -mem-lim=<MB>

    It prints:
      v <selected-positive-lits> 0
      Optimal

    or:
      Infeasible

    For the outer binary search we only need:
      cost <= p  => feasible
      cost >  p  => infeasible
    """
    model = []
    final_status = None

    for raw in stdout.splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith("v ") or line.startswith("V "):
            for tok in line.split()[1:]:
                if tok == "0":
                    continue
                try:
                    model.append(int(tok))
                except ValueError:
                    pass

        elif line == "Optimal":
            final_status = "Optimal"

        elif line == "Infeasible":
            final_status = "Infeasible"

        elif "Out of Memory" in line or "Memory limit" in line:
            return "error", None, []

    if final_status == "Optimal":
        # Feasible with <= p centers. Returning p_bound is enough for feasibility.
        return "optimum", int(p_bound), model

    if final_status == "Infeasible":
        # Proven infeasible with <= p centers.
        return "optimum", int(p_bound) + 1, model

    return "error", None, model

def _wcnf_to_data(wcnf):
    """
    Convert WCNF into plain Python lists so the RC2 worker can be spawned safely.
    """
    return {
        "nv": int(getattr(wcnf, "nv", 0)),
        "hard": [list(map(int, cl)) for cl in getattr(wcnf, "hard", [])],
        "soft": [list(map(int, cl)) for cl in getattr(wcnf, "soft", [])],
        "weights": [int(w) for w in getattr(wcnf, "wght", [])],
    }


def _wcnf_from_data(data):
    """
    Rebuild a PySAT WCNF object from plain Python lists.
    """
    wcnf = WCNF()

    for cl in data["hard"]:
        wcnf.append(cl)

    for cl, weight in zip(data["soft"], data["weights"]):
        wcnf.append(cl, weight=weight)

    wcnf.nv = int(data["nv"])
    return wcnf


def _pysat_maxsat_worker(wcnf_data, solver_name, send_conn):
    """
    Run a PySAT MaxSAT solver in a child process so the parent can enforce
    a hard timeout.
    """
    try:
        wcnf = _wcnf_from_data(wcnf_data)

        cpu0 = cpu_self_seconds()
        wall0 = time.perf_counter()

        if solver_name in PYSAT_MAXSAT_SOLVERS:
            with RC2(wcnf) as solver:
                model = solver.compute()
                cost = solver.cost
        else:
            raise ValueError(f"Unsupported PySAT MaxSAT solver: {solver_name}")

        cpu1 = cpu_self_seconds()
        wall = time.perf_counter() - wall0
        cpu = (cpu1 - cpu0) if (cpu0 is not None and cpu1 is not None) else wall

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
    """
    Solve a weighted partial MaxSAT instance using a PySAT MaxSAT solver.

    The solver is executed in a child process so the parent process can enforce
    the experiment time limit.
    """
    recv_conn, send_conn = mp.Pipe(duplex=False)
    proc = mp.Process(
        target=_pysat_maxsat_worker,
        args=(_wcnf_to_data(wcnf), solver_name, send_conn),
        daemon=True,
    )

    wall0 = time.perf_counter()
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
                wall = time.perf_counter() - wall0
                _terminate_process_fast(proc)
                return "timeout", None, [], wall

            msg = recv_conn.recv()

        proc.join(timeout=0.2)

        if not msg:
            return "error", None, [], time.perf_counter() - wall0

        status, cost, model, cpu, *rest = msg

        if status == "error" and rest:
            print(rest[0], file=sys.stderr, flush=True)

        return status, cost, model, cpu

    except EOFError:
        return "error", None, [], time.perf_counter() - wall0

    finally:
        _terminate_process_fast(proc)
        try:
            recv_conn.close()
        except Exception:
            pass

def _solve_wcnf(wcnf, solver_name, time_limit, p_bound=None):
    if solver_name in PYSAT_MAXSAT_SOLVERS:
        return _solve_wcnf_pysat(wcnf, solver_name, time_limit)

    if solver_name in EXTERNAL_MAXSAT_SOLVERS:
        return _run_external_maxsat_solver(
            solver_name=solver_name,
            wcnf=wcnf,
            time_limit=time_limit,
            p_bound=p_bound,
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
    """
    MaxSAT-based feasibility solving for p-Center.

    For one tested radius R:
      - hard clauses = coverage(R)
      - soft clauses = (-y_j), weight 1, for all candidate centers j

    The MaxSAT optimum cost is the minimum number of centers required to cover all
    nodes under radius R. Therefore:
      cost <= p  => R is feasible
      cost >  p  => R is infeasible

    The outer optimization over R is still performed by binary search over the
    descending candidate radius list.
    """
    if encoding not in ("maxsat_setcover", "maxsat_cover"):
        raise ValueError(
            "MaxSAT set-cover mode expects --encodings maxsat_setcover "
            "(or alias maxsat_cover)."
        )

    search_t0 = time.perf_counter()

    lo = 0
    hi = len(inst.radii) - 1

    best_radius = None
    best_centers = None
    last_nvars = None
    last_nclauses = None
    last_cpu = None
    final_status = "infeasible"

    while lo <= hi:
        mid = (lo + hi) // 2
        radius = inst.radii[mid]

        if cancel_dict is not None and cancel_dict.get("cancel", False):
            return (
                "cancelled",
                best_radius,
                last_nvars,
                last_nclauses,
                best_centers,
                last_cpu,
                time.perf_counter() - search_t0,
            )

        wcnf, info = inst._encode_wcnf_maxsat_setcover(radius)
        if wcnf is None:
            hi = mid - 1
            continue

        status, cost, model, cpu_sec = _solve_wcnf(
            wcnf, 
            solver_name, 
            time_limit,
            p_bound=inst.p
        )

        last_nvars = wcnf.nv
        last_nclauses = len(wcnf.hard) + len(wcnf.soft)
        last_cpu = cpu_sec

        if status not in ("optimum", "sat"):
            return (
                status,
                best_radius,
                last_nvars,
                last_nclauses,
                best_centers,
                last_cpu,
                time.perf_counter() - search_t0,
            )

        centers = _decode_centers(model, info["y"])

        if cost is None:
            cost = len(centers)

        if cost <= inst.p:
            # Feasible radius: try smaller radii.
            best_radius = radius
            best_centers = centers
            final_status = "OK"
            lo = mid + 1
        else:
            # Infeasible radius: need larger radius.
            hi = mid - 1

    return (
        final_status,
        best_radius,
        last_nvars,
        last_nclauses,
        best_centers,
        last_cpu,
        time.perf_counter() - search_t0,
    )
