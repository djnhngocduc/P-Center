try:
    from docplex.cp.model import CpoModel
    from docplex.cp.parameters import CpoParameters
except Exception:
    CpoModel = None
    CpoParameters = None


def _solve_radius_cpo(inst, idx, radius, time_limit, data=None):
    import os
    from experiments import _cpu_self_seconds
    
    pid = os.getpid()

    if CpoModel is None or CpoParameters is None:
        raise ImportError(
            "Could not import docplex.cp.model.CpoModel. "
            "Install docplex and make sure CP Optimizer is available."
        )

    if data is None:
        data = inst._build_setcover_data(radius)

    if data is None:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} -> UNSAT(CPO:no coverage)",
            flush=True
        )
        return idx, radius, "unsat", 0.0, inst.n, 0, None

    n = data["n"]
    p = data["p"]
    cover_rows = data["cover_rows"]

    mdl = CpoModel(name=f"pcenter_feas_R_{radius}")
    x = mdl.binary_var_list(n, name="x")

    rows_sorted = sorted(cover_rows, key=lambda t: len(t[1]))
    for u, allowed in rows_sorted:
        mdl.add(mdl.sum(x[j] for j in allowed) >= 1)

    total_open = mdl.sum(x)
    mdl.add(total_open <= p)

    mdl.set_search_phases([mdl.search_phase(x)])

    params = CpoParameters()
    params.LogVerbosity = "Quiet"
    params.TimeLimit = float(time_limit) if (time_limit and time_limit > 0) else None
    params.TimeMode = "ElapsedTime"
    params.Presolve = "On"
    params.DefaultInferenceLevel = "Extended"
    params.DynamicProbing = "On"
    params.SearchType = "Auto"
    params.RandomSeed = 0

    cpu0 = _cpu_self_seconds()
    try:
        sol = mdl.solve(params=params)
    except Exception as e:
        cpu1 = _cpu_self_seconds()
        cpu_sec = cpu1 - cpu0
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} -> ERROR(CPO) cpu={cpu_sec:.6f}s err={e}",
            flush=True
        )
        return idx, radius, "error", cpu_sec, n, len(cover_rows) + 1, None

    cpu1 = _cpu_self_seconds()
    cpu_sec = cpu1 - cpu0

    nvars = n
    nclauses = len(cover_rows) + 1

    if sol is None:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} -> TIMEOUT(CPO) cpu={cpu_sec:.6f}s",
            flush=True
        )
        return idx, radius, "timeout", cpu_sec, nvars, nclauses, None

    try:
        solve_status = str(sol.get_solve_status()).lower()
    except Exception:
        solve_status = "unknown"

    try:
        fail_status = str(sol.get_fail_status()).lower()
    except Exception:
        fail_status = ""

    centers = []
    has_solution = False

    try:
        if sol.is_solution():
            vals = [sol[x[j]] for j in range(n)]
            centers = [j for j, v in enumerate(vals) if v is not None and float(v) > 0.5]
            has_solution = True
    except Exception:
        try:
            vals = [sol[x[j]] for j in range(n)]
            centers = [j for j, v in enumerate(vals) if v is not None and float(v) > 0.5]
            has_solution = True
        except Exception:
            has_solution = False

    if has_solution or ("feasible" in solve_status) or ("optimal" in solve_status):
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> SAT(CPO) cpu={cpu_sec:.6f}s status={solve_status} centers={centers}",
            flush=True
        )
        return idx, radius, "sat", cpu_sec, nvars, nclauses, centers

    if ("infeasible" in solve_status) or ("searchcompleted" in fail_status and not has_solution):
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> UNSAT(CPO) cpu={cpu_sec:.6f}s status={solve_status} fail={fail_status}",
            flush=True
        )
        return idx, radius, "unsat", cpu_sec, nvars, nclauses, None

    print(
        f"[WORKER-END] pid={pid} idx={idx} R={radius} "
        f"-> TIMEOUT(CPO) cpu={cpu_sec:.6f}s status={solve_status} fail={fail_status}",
        flush=True
    )
    return idx, radius, "timeout", cpu_sec, nvars, nclauses, None
