import time
import threading

from pysat.solvers import Solver

from solvers.backend import EXTERNAL_SOLVERS

def search_min_radius_incremental(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    seed_idx=None,
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

    if solver_name in EXTERNAL_SOLVERS:
        raise ValueError(
            f"Incremental mode does not support external solver '{solver_name}'. "
            f"Use an internal PySAT solver such as glucose4/maplecm/maplechrono."
        )

    base_cnf, info = inst._build_incremental_base_cnf(encoding=encoding)
    radii = info["radii"]
    y_vars = info["y"]

    nR = len(radii)
    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, base_cnf.nv, len(base_cnf.clauses), None, search_elapsed

    lo = 0
    hi = nR - 1

    best_sat_idx = None
    best_centers = None

    next_var = base_cnf.nv + 1
    tested = {}
    decided = {}

    def _timeout_result():
        search_elapsed = time.perf_counter() - search_t0
        nvars = next_var - 1
        nclauses = len(base_cnf.clauses) + (len(tested) * inst.n)
        if best_sat_idx is not None:
            return (
                "timeout_with_incumbent",
                radii[best_sat_idx],
                nvars,
                nclauses,
                best_centers,
                search_elapsed,
            )
        return "timeout", None, nvars, nclauses, None, search_elapsed

    print(
        f"[INC-INIT] encoding={encoding} solver={solver_name} p={inst.p} "
        f"lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={nR} "
        f"base_clauses={len(base_cnf.clauses)} base_vars={base_cnf.nv}",
        flush=True
    )

    with Solver(name=solver_name, bootstrap_with=base_cnf.clauses) as solver:
        while lo <= hi:
            step_time_limit = _remaining_time()
            if deadline is not None and step_time_limit <= 0:
                print("[INC-TIMEOUT] instance-level time budget exhausted", flush=True)
                return _timeout_result()

            mid = (lo + hi) // 2
            R = radii[mid]

            if mid not in tested:
                alpha = next_var
                next_var += 1
                tested[mid] = alpha

                step_clauses = inst._build_incremental_guarded_clauses(
                    radius=R,
                    selector_lit=alpha,
                )
                for cl in step_clauses:
                    solver.add_clause(cl)

                print(
                    f"[INC-ADD] idx={mid} R={R} selector={alpha} "
                    f"added_clauses={len(step_clauses)}",
                    flush=True
                )
            else:
                alpha = tested[mid]

            step_time_limit = _remaining_time()
            if deadline is not None and step_time_limit <= 0:
                print("[INC-TIMEOUT] instance-level time budget exhausted", flush=True)
                return _timeout_result()

            print(
                f"[INC-STEP] lo={lo} hi={hi} mid={mid} R={R}",
                flush=True
            )

            timer = None
            try:
                if step_time_limit and step_time_limit > 0 and hasattr(solver, "interrupt"):
                    timer = threading.Timer(step_time_limit, solver.interrupt)
                    timer.start()
                    try:
                        sat = solver.solve_limited(
                            assumptions=[alpha],
                            expect_interrupt=True
                        )
                    except NotImplementedError:
                        sat = solver.solve(assumptions=[alpha])
                else:
                    sat = solver.solve(assumptions=[alpha])
            finally:
                if timer:
                    timer.cancel()

            if sat is None:
                decided[mid] = "timeout"
                print(
                    f"[INC-DONE] idx={mid} R={R} status=timeout",
                    flush=True
                )
                return _timeout_result()

            if sat:
                decided[mid] = "sat"
                model = solver.get_model() or []
                model_set = set(model)
                centers = sorted([j for j, v in enumerate(y_vars) if v in model_set])

                best_sat_idx = mid
                best_centers = centers

                print(
                    f"[INC-DONE] idx={mid} R={R} status=sat centers={centers}",
                    flush=True
                )

                lo = mid + 1
            else:
                decided[mid] = "unsat"
                print(
                    f"[INC-DONE] idx={mid} R={R} status=unsat",
                    flush=True
                )

                hi = mid - 1

    if best_sat_idx is None:
        print("[INC-RESULT] infeasible (no SAT found)", flush=True)
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, next_var - 1, None, None, search_elapsed

    best_radius = radii[best_sat_idx]
    nvars = next_var - 1
    nclauses = len(base_cnf.clauses) + (len(tested) * inst.n)

    cert_unsat_idx = best_sat_idx + 1
    is_smallest_candidate = best_sat_idx == nR - 1
    certified = is_smallest_candidate or (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        final_status = "OK"
        if is_smallest_candidate:
            print(
                f"[INC-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={best_radius} "
                f"is the smallest candidate radius",
                flush=True,
            )
        else:
            print(
                f"[INC-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={best_radius}, "
                f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={radii[cert_unsat_idx]}",
                flush=True
            )
    else:
        final_status = "uncertified"
        print(
            f"[INC-FALLBACK] best SAT found but UNSAT boundary not certified; "
            f"best_sat_idx={best_sat_idx}, best_R={best_radius}",
            flush=True
        )

    print(
        f"[INC-RESULT] status={final_status} best_idx={best_sat_idx} "
        f"best_R={best_radius}",
        flush=True
    )

    search_elapsed = time.perf_counter() - search_t0
    return (
        final_status,
        best_radius,
        nvars,
        nclauses,
        best_centers,
        search_elapsed,
    )
