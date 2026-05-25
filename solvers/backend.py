import os
import sys
import tempfile
import threading
import time
import subprocess

try:
    import cplex
except Exception:
    cplex = None

try:
    import gurobipy as gp
    from gurobipy import GRB
except Exception:
    gp = None
    GRB = None

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

EXTERNAL_SOLVERS = {
    "kissat": (
        os.path.join(BASE_DIR, "kissat", "build", "kissat"),
        []
    )
}

def write_dimacs(cnf, path):
    nvars = cnf.nv
    nclauses = len(cnf.clauses)

    lines = [f"p cnf {nvars} {nclauses}\n"]
    for cl in cnf.clauses:
        if not cl:
            lines.append("0\n")
        else:
            lines.append(" ".join(str(l) for l in cl) + " 0\n")

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)


def run_external_solver(solver_name, cnf, time_limit, cancel_ev=None, tmpdir=None):
    call_t0 = time.perf_counter()

    def _remaining_time():
        if time_limit and time_limit > 0:
            return max(0.0, float(time_limit) - (time.perf_counter() - call_t0))
        return time_limit

    bin_path, extra_args = EXTERNAL_SOLVERS[solver_name]

    if (not os.path.exists(bin_path)) or (not os.access(bin_path, os.X_OK)):
        print(
            f"[SOLVER-ERROR] solver={solver_name}: binary not found or not executable: {bin_path}",
            file=sys.stderr,
            flush=True
        )
        return "error", None

    base_tmp = tmpdir if tmpdir is not None else ("/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir())

    cnf_path = os.path.join(
        base_tmp,
        f"pcsat_{os.getpid()}_{threading.get_ident()}.cnf",
    )

    write_dimacs(cnf, cnf_path)
    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        try:
            if os.path.exists(cnf_path):
                os.remove(cnf_path)
        except Exception:
            pass
        return "timeout", None

    solver_dir = os.path.dirname(bin_path)

    env = os.environ.copy()
    if tmpdir is not None:
        env["TMPDIR"] = tmpdir

    cmd = [bin_path] + extra_args + [cnf_path]

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

    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        try:
            proc.kill()
        except Exception:
            pass
        try:
            if os.path.exists(cnf_path):
                os.remove(cnf_path)
        except Exception:
            pass
        return "timeout", None

    try:
        stdout, stderr = proc.communicate(
            timeout=step_time_limit if (step_time_limit and step_time_limit > 0) else None
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

        print(f"[SOLVER-TIMEOUT] solver={solver_name} cmd={cmd}", flush=True)
        return "timeout", None
    finally:
        try:
            if os.path.exists(cnf_path):
                os.remove(cnf_path)
        except Exception:
            pass

    status = None
    model_lits = []

    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("s "):
            low = line.lower()
            if "unsat" in low:
                status = "unsat"
            elif "sat" in low:
                status = "sat"
        elif line.startswith("v ") or line.startswith("V "):
            parts = line.split()[1:]
            for tok in parts:
                if tok == "0":
                    continue
                try:
                    lit = int(tok)
                    model_lits.append(lit)
                except ValueError:
                    continue

    if status is None:
        if model_lits:
            status = "sat"
        else:
            status = "error" if proc.returncode not in (0, 10, 20) else "unsat"

    if status == "sat" and not model_lits:
        print(f"[SOLVER-WARN] {solver_name} reported SAT but no model.", flush=True)
        return "error", None

    return status, (model_lits if status == "sat" else None)


def solve_radius_cplex(inst, idx, radius, time_limit, data=None, current_cancel_ref=None):
    call_t0 = time.perf_counter()

    def _remaining_time():
        if time_limit and time_limit > 0:
            return max(0.0, float(time_limit) - (time.perf_counter() - call_t0))
        return time_limit

    pid = os.getpid()

    if cplex is None:
        raise ImportError(
            "Could not import the IBM CPLEX Python API. "
            "Make sure cplex_studio2211 is installed and PYTHONPATH includes "
            ".../cplex/python/<python-version>/x86-64_linux"
        )

    if data is None:
        data = inst._build_setcover_data(radius)

    if data is None:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} -> UNSAT(no coverage)",
            flush=True
        )
        return idx, radius, "unsat", inst.n, 0, None

    n = data["n"]
    p = data["p"]
    cover_rows = data["cover_rows"]
    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        return idx, radius, "timeout", n, len(cover_rows) + 1, None

    model = cplex.Cplex()
    aborter = None

    model.set_results_stream(None)
    model.set_warning_stream(None)
    model.set_error_stream(None)
    model.set_log_stream(None)
    model.parameters.threads.set(1)

    if step_time_limit and step_time_limit > 0:
        model.parameters.timelimit.set(float(step_time_limit))

    try:
        aborter = cplex.Aborter()
        model.use_aborter(aborter)
        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = aborter.abort
    except Exception:
        aborter = None
        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = None

    names = [f"x_{j}" for j in range(n)]
    model.variables.add(
        obj=[0.0] * n,
        lb=[0.0] * n,
        ub=[1.0] * n,
        types=[model.variables.type.binary] * n,
        names=names,
    )

    for u, allowed in cover_rows:
        model.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=allowed, val=[1.0] * len(allowed))],
            senses=["G"],
            rhs=[1.0],
            names=[f"cover_{u}"],
        )

    model.linear_constraints.add(
        lin_expr=[cplex.SparsePair(ind=list(range(n)), val=[1.0] * n)],
        senses=["L"],
        rhs=[float(p)],
        names=["atmost_p"],
    )

    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        nvars = model.variables.get_num()
        nclauses = model.linear_constraints.get_num()
        return idx, radius, "timeout", nvars, nclauses, None
    if step_time_limit and step_time_limit > 0:
        model.parameters.timelimit.set(float(step_time_limit))

    model.solve()

    nvars = model.variables.get_num()
    nclauses = model.linear_constraints.get_num()

    status_code = model.solution.get_status()
    try:
        status_name = model.solution.get_status_string(status_code)
    except Exception:
        try:
            status_name = model.solution.status[status_code]
        except Exception:
            status_name = str(status_code)
    status_name = str(status_name).lower()

    try:
        primal_feasible = model.solution.is_primal_feasible()
    except Exception:
        primal_feasible = False

    try:
        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = None
    except Exception:
        pass

    if primal_feasible or ("optimal" in status_name) or ("feasible" in status_name and "infeasible" not in status_name):
        vals = model.solution.get_values()
        centers = [j for j, v in enumerate(vals) if v > 0.5]
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> SAT(CPLEX) status={status_name} centers={centers}",
            flush=True
        )
        return idx, radius, "sat", nvars, nclauses, centers

    if "infeasible" in status_name:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> UNSAT(CPLEX) status={status_name}",
            flush=True
        )
        return idx, radius, "unsat", nvars, nclauses, None

    if ("abort" in status_name) or ("interrupted" in status_name):
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> CANCELLED(CPLEX) status={status_name}",
            flush=True
        )
        return idx, radius, "cancelled", nvars, nclauses, None

    if ("time limit" in status_name) or ("limit" in status_name):
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> TIMEOUT(CPLEX) status={status_name}",
            flush=True
        )
        return idx, radius, "timeout", nvars, nclauses, None

    print(
        f"[WORKER-END] pid={pid} idx={idx} R={radius} "
        f"-> ERROR(CPLEX) status={status_name}",
        flush=True
    )
    return idx, radius, "error", nvars, nclauses, None


def solve_radius_gurobi(inst, idx, radius, time_limit, data=None, env=None, current_cancel_ref=None):
    call_t0 = time.perf_counter()

    def _remaining_time():
        if time_limit and time_limit > 0:
            return max(0.0, float(time_limit) - (time.perf_counter() - call_t0))
        return time_limit

    pid = os.getpid()

    if gp is None or GRB is None:
        raise ImportError(
            "Could not import gurobipy. Make sure Gurobi is installed and the license is available."
        )

    if data is None:
        data = inst._build_setcover_data(radius)

    if data is None:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} -> UNSAT(GUROBI:no coverage)",
            flush=True
        )
        return idx, radius, "unsat", inst.n, 0, None

    n = data["n"]
    p = data["p"]
    cover_rows = data["cover_rows"]
    step_time_limit = _remaining_time()
    if time_limit and time_limit > 0 and step_time_limit <= 0:
        return idx, radius, "timeout", n, len(cover_rows) + 1, None

    model = None
    try:
        if env is None:
            model = gp.Model(f"pcenter_feas_R_{radius}")
        else:
            model = gp.Model(f"pcenter_feas_R_{radius}", env=env)

        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = model.terminate

        model.Params.OutputFlag = 0
        model.Params.Threads = 1

        if step_time_limit and step_time_limit > 0:
            model.Params.TimeLimit = float(step_time_limit)

        x = model.addVars(n, vtype=GRB.BINARY, name="x")

        for u, allowed in cover_rows:
            model.addConstr(
                gp.quicksum(x[j] for j in allowed) >= 1,
                name=f"cover_{u}"
            )

        model.addConstr(
            gp.quicksum(x[j] for j in range(n)) <= p,
            name="atmost_p"
        )

        step_time_limit = _remaining_time()
        if time_limit and time_limit > 0 and step_time_limit <= 0:
            model.update()
            return idx, radius, "timeout", model.NumVars, model.NumConstrs, None
        if step_time_limit and step_time_limit > 0:
            model.Params.TimeLimit = float(step_time_limit)

        model.optimize()

        nvars = model.NumVars
        nclauses = model.NumConstrs
        status = model.Status
        solcount = int(getattr(model, "SolCount", 0) or 0)

        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = None

        if solcount > 0:
            centers = [j for j in range(n) if x[j].X > 0.5]
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> SAT(GUROBI) status={status} centers={centers}",
                flush=True
            )
            return idx, radius, "sat", nvars, nclauses, centers

        if status == GRB.INFEASIBLE:
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> UNSAT(GUROBI) status={status}",
                flush=True
            )
            return idx, radius, "unsat", nvars, nclauses, None

        if status in (
            GRB.TIME_LIMIT,
            GRB.NODE_LIMIT,
            GRB.ITERATION_LIMIT,
            GRB.SOLUTION_LIMIT,
        ):
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> TIMEOUT(GUROBI) status={status}",
                flush=True
            )
            return idx, radius, "timeout", nvars, nclauses, None

        if status == GRB.INTERRUPTED:
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> CANCELLED(GUROBI) status={status}",
                flush=True
            )
            return idx, radius, "cancelled", nvars, nclauses, None

        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> ERROR(GUROBI) status={status}",
            flush=True
        )
        return idx, radius, "error", nvars, nclauses, None

    finally:
        if current_cancel_ref is not None:
            current_cancel_ref["fn"] = None
        if model is not None:
            try:
                model.dispose()
            except Exception:
                pass

