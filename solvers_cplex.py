try:
    import cplex
except Exception:
    cplex = None


def _solve_radius_cplex(inst, idx, radius, time_limit, data=None):
    import os
    from experiments import _cpu_self_seconds
    
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
        return idx, radius, "unsat", 0.0, inst.n, 0, None

    n = data["n"]
    p = data["p"]
    cover_rows = data["cover_rows"]

    model = cplex.Cplex()
    model.set_results_stream(None)
    model.set_warning_stream(None)
    model.set_error_stream(None)
    model.set_log_stream(None)
    model.parameters.threads.set(1)

    if time_limit and time_limit > 0:
        model.parameters.timelimit.set(float(time_limit))

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

    cpu0 = _cpu_self_seconds()
    model.solve()
    cpu1 = _cpu_self_seconds()
    cpu_sec = cpu1 - cpu0

    status_code = model.solution.get_status()
    nvars = model.variables.get_num()
    nclauses = model.linear_constraints.get_num()

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

    if primal_feasible or ("optimal" in status_name) or ("feasible" in status_name and "infeasible" not in status_name):
        vals = model.solution.get_values()
        centers = [j for j, v in enumerate(vals) if v > 0.5]
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> SAT(CPLEX) cpu={cpu_sec:.6f}s status={status_name} centers={centers}",
            flush=True
        )
        return idx, radius, "sat", cpu_sec, nvars, nclauses, centers

    if "infeasible" in status_name:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> UNSAT(CPLEX) cpu={cpu_sec:.6f}s status={status_name}",
            flush=True
        )
        return idx, radius, "unsat", cpu_sec, nvars, nclauses, None

    if ("time limit" in status_name) or ("abort" in status_name) or ("limit" in status_name):
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> TIMEOUT(CPLEX) cpu={cpu_sec:.6f}s status={status_name}",
            flush=True
        )
        return idx, radius, "timeout", cpu_sec, nvars, nclauses, None

    print(
        f"[WORKER-END] pid={pid} idx={idx} R={radius} "
        f"-> ERROR(CPLEX) cpu={cpu_sec:.6f}s status={status_name}",
        flush=True
    )
    return idx, radius, "error", cpu_sec, nvars, nclauses, None
