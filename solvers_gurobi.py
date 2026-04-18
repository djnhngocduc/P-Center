try:
    import gurobipy as gp
    from gurobipy import GRB
except Exception:
    gp = None
    GRB = None


def _solve_radius_gurobi(inst, idx, radius, time_limit, data=None):
    import os
    from experiments import _cpu_self_seconds
    
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
        return idx, radius, "unsat", 0.0, inst.n, 0, None

    n = data["n"]
    p = data["p"]
    cover_rows = data["cover_rows"]

    model = None
    try:
        model = gp.Model(f"pcenter_feas_R_{radius}")
        model.Params.OutputFlag = 0
        model.Params.Threads = 1

        if time_limit and time_limit > 0:
            model.Params.TimeLimit = float(time_limit)

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

        cpu0 = _cpu_self_seconds()
        model.optimize()
        cpu1 = _cpu_self_seconds()
        cpu_sec = cpu1 - cpu0

        nvars = model.NumVars
        nclauses = model.NumConstrs
        status = model.Status
        solcount = int(getattr(model, "SolCount", 0) or 0)

        if solcount > 0:
            centers = [j for j in range(n) if x[j].X > 0.5]
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> SAT(GUROBI) cpu={cpu_sec:.6f}s status={status} centers={centers}",
                flush=True
            )
            return idx, radius, "sat", cpu_sec, nvars, nclauses, centers

        if status == GRB.INFEASIBLE:
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> UNSAT(GUROBI) cpu={cpu_sec:.6f}s status={status}",
                flush=True
            )
            return idx, radius, "unsat", cpu_sec, nvars, nclauses, None

        if status in (
            GRB.TIME_LIMIT,
            GRB.NODE_LIMIT,
            GRB.ITERATION_LIMIT,
            GRB.SOLUTION_LIMIT,
        ):
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> TIMEOUT(GUROBI) cpu={cpu_sec:.6f}s status={status}",
                flush=True
            )
            return idx, radius, "timeout", cpu_sec, nvars, nclauses, None

        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} "
            f"-> ERROR(GUROBI) cpu={cpu_sec:.6f}s status={status}",
            flush=True
        )
        return idx, radius, "error", cpu_sec, nvars, nclauses, None

    finally:
        if model is not None:
            try:
                model.dispose()
            except Exception:
                pass
