import time
import threading
from concurrent.futures import ProcessPoolExecutor
import traceback
import sys

from pysat.solvers import Solver


def search_min_radius_hybrid_threshold(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    hybrid_threshold,
    radii_workers,
    seed_idx=None,
    mgr=None,
    cancel_dict=None
):
    search_t0 = time.perf_counter()

    if solver_name == "cplex_mip":
        raise ValueError("hybrid_threshold expects a SAT solver in --solvers, not 'cplex_mip'.")

    from experiments import EXTERNAL_SOLVERS
    try:
        import cplex
    except Exception:
        cplex = None
    
    if cplex is None:
        raise ImportError("hybrid_threshold needs CPLEX installed.")

    use_incremental_sat = solver_name not in EXTERNAL_SOLVERS

    if use_incremental_sat:
        base_cnf, info = inst._build_incremental_base_cnf(encoding=encoding)
        radii = info["radii"]
        y_vars = info["y"]

        if not radii:
            search_elapsed = time.perf_counter() - search_t0
            return "infeasible", None, base_cnf.nv, len(base_cnf.clauses), None, None, search_elapsed

        lo = 0
        hi = len(radii) - 1
        next_var = base_cnf.nv + 1
        tested = {}
        added_clause_count = 0

        print(
            f"[HYBRID-INIT] encoding={encoding} solver={solver_name} sat_mode=incremental "
            f"threshold={hybrid_threshold} p={inst.p} lo={lo} hi={hi} "
            f"R_lo={radii[lo]} R_hi={radii[hi]} nR={len(radii)} "
            f"base_clauses={len(base_cnf.clauses)} base_vars={base_cnf.nv}",
            flush=True
        )

        best_sat_idx = None
        best_centers = None
        best_cpu = None
        best_nvars = None
        best_nclauses = None

        with Solver(name=solver_name, bootstrap_with=base_cnf.clauses) as sat_solver:
            while lo <= hi:
                mid = (lo + hi) // 2
                R = radii[mid]
                from experiments import _remaining_iters_from_bounds
                remaining_iters = _remaining_iters_from_bounds(lo, hi)
                use_sat = remaining_iters > hybrid_threshold

                print(
                    f"[HYBRID-STEP] lo={lo} hi={hi} mid={mid} R={R} "
                    f"remaining_iters={remaining_iters} backend={'SAT' if use_sat else 'CPLEX'}",
                    flush=True
                )

                if use_sat:
                    if mid not in tested:
                        alpha = next_var
                        next_var += 1
                        tested[mid] = alpha

                        step_clauses = inst._build_incremental_guarded_clauses(
                            radius=R,
                            selector_lit=alpha,
                        )
                        for cl in step_clauses:
                            sat_solver.add_clause(cl)
                        added_clause_count += len(step_clauses)

                        print(
                            f"[HYBRID-SAT-ADD] idx={mid} R={R} selector={alpha} "
                            f"added_clauses={len(step_clauses)}",
                            flush=True
                        )
                    else:
                        alpha = tested[mid]

                    timer = None
                    cpu0 = None
                    cpu1 = None
                    try:
                        from experiments import _cpu_self_seconds
                        if time_limit and time_limit > 0 and hasattr(sat_solver, "interrupt"):
                            timer = threading.Timer(time_limit, sat_solver.interrupt)
                            timer.start()
                            cpu0 = _cpu_self_seconds()
                            try:
                                sat = sat_solver.solve_limited(
                                    assumptions=[alpha],
                                    expect_interrupt=True
                                )
                            except NotImplementedError:
                                sat = sat_solver.solve(assumptions=[alpha])
                            cpu1 = _cpu_self_seconds()
                        else:
                            cpu0 = _cpu_self_seconds()
                            sat = sat_solver.solve(assumptions=[alpha])
                            cpu1 = _cpu_self_seconds()
                    finally:
                        if timer:
                            timer.cancel()

                    cpu_sec = (cpu1 - cpu0) if (cpu0 is not None and cpu1 is not None) else 0.0

                    if sat is None:
                        print(
                            f"[HYBRID-SAT-DONE] idx={mid} R={R} status=timeout cpu={cpu_sec:.6f}s",
                            flush=True
                        )
                        search_elapsed = time.perf_counter() - search_t0
                        return "timeout", None, next_var - 1, None, None, cpu_sec, search_elapsed

                    if sat:
                        model = sat_solver.get_model() or []
                        model_set = set(model)
                        centers = sorted([j for j, v in enumerate(y_vars) if v in model_set])

                        best_sat_idx = mid
                        best_centers = centers
                        best_cpu = cpu_sec
                        best_nvars = next_var - 1
                        best_nclauses = len(base_cnf.clauses) + added_clause_count

                        print(
                            f"[HYBRID-SAT-DONE] idx={mid} R={R} status=sat "
                            f"cpu={cpu_sec:.6f}s centers={centers}",
                            flush=True
                        )
                        lo = mid + 1
                    else:
                        print(
                            f"[HYBRID-SAT-DONE] idx={mid} R={R} status=unsat cpu={cpu_sec:.6f}s",
                            flush=True
                        )
                        hi = mid - 1

                else:
                    from experiments import _solve_radius_cplex
                    idx, Rret, status, cpu_sec, nvars_cpx, nclauses_cpx, centers = _solve_radius_cplex(
                        inst=inst,
                        idx=mid,
                        radius=R,
                        time_limit=time_limit,
                        data=None,
                    )

                    print(
                        f"[HYBRID-CPLEX-DONE] idx={idx} R={Rret} status={status} "
                        f"cpu={cpu_sec:.6f}s nvars={nvars_cpx} nclauses={nclauses_cpx}",
                        flush=True
                    )

                    if status == "sat":
                        best_sat_idx = idx
                        best_centers = centers
                        best_cpu = cpu_sec
                        best_nvars = nvars_cpx
                        best_nclauses = nclauses_cpx
                        lo = idx + 1

                    elif status == "unsat":
                        hi = idx - 1

                    elif status in ("timeout", "error"):
                        if best_sat_idx is not None:
                            print(
                                f"[HYBRID-FALLBACK] stop on {status}; "
                                f"use best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                                flush=True
                            )
                            break

                        search_elapsed = time.perf_counter() - search_t0
                        return status, None, nvars_cpx, nclauses_cpx, None, cpu_sec, search_elapsed

        if best_sat_idx is None:
            print("[HYBRID-RESULT] infeasible (no SAT found)", flush=True)
            search_elapsed = time.perf_counter() - search_t0
            return "infeasible", None, None, None, None, None, search_elapsed

        best_radius = radii[best_sat_idx]
        print(
            f"[HYBRID-RESULT] status=OK best_idx={best_sat_idx} "
            f"best_R={best_radius} cpu={best_cpu} threshold={hybrid_threshold}",
            flush=True
        )

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

    else:
        radii = inst.radii

        if not radii:
            search_elapsed = time.perf_counter() - search_t0
            return "infeasible", None, None, None, None, None, search_elapsed

        lo = 0
        hi = len(radii) - 1

        print(
            f"[HYBRID-INIT] encoding={encoding} solver={solver_name} sat_mode=external_per_radius "
            f"threshold={hybrid_threshold} p={inst.p} lo={lo} hi={hi} "
            f"R_lo={radii[lo]} R_hi={radii[hi]} nR={len(radii)}",
            flush=True
        )

        best_sat_idx = None
        best_centers = None
        best_cpu = None
        best_nvars = None
        best_nclauses = None

        from experiments import _pool_initializer
        with ProcessPoolExecutor(
            max_workers=1,
            initializer=_pool_initializer,
            initargs=(cancel_dict, inst),
        ) as ex:
            while lo <= hi:
                mid = (lo + hi) // 2
                R = radii[mid]
                from experiments import _remaining_iters_from_bounds
                remaining_iters = _remaining_iters_from_bounds(lo, hi)
                use_sat = remaining_iters > hybrid_threshold

                print(
                    f"[HYBRID-STEP] lo={lo} hi={hi} mid={mid} R={R} "
                    f"remaining_iters={remaining_iters} backend={'SAT' if use_sat else 'CPLEX'}",
                    flush=True
                )

                if use_sat:
                    from experiments import _solve_radius_worker_proc
                    fut = ex.submit(
                        _solve_radius_worker_proc,
                        mid,
                        encoding,
                        solver_name,
                        R,
                        time_limit,
                    )
                    try:
                        idx, Rret, status, cpu_sec, nvars_sat, nclauses_sat, centers = fut.result()
                    except Exception:
                        tb = traceback.format_exc()
                        print(
                            f"[HYBRID-EXTERNAL-EXCEPTION] idx={mid} R={R}\n{tb}",
                            file=sys.stderr,
                            flush=True
                        )
                        idx, Rret, status, cpu_sec, nvars_sat, nclauses_sat, centers = (
                            mid, R, "error", 0.0, None, None, None
                        )

                    print(
                        f"[HYBRID-SAT-DONE] idx={idx} R={Rret} status={status} "
                        f"cpu={cpu_sec:.6f}s nvars={nvars_sat} nclauses={nclauses_sat}",
                        flush=True
                    )

                    if status == "sat":
                        best_sat_idx = idx
                        best_centers = centers
                        best_cpu = cpu_sec
                        best_nvars = nvars_sat
                        best_nclauses = nclauses_sat
                        lo = idx + 1

                    elif status == "unsat":
                        hi = idx - 1

                    elif status in ("timeout", "error"):
                        if best_sat_idx is not None:
                            print(
                                f"[HYBRID-FALLBACK] stop on {status}; "
                                f"use best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                                flush=True
                            )
                            break

                        search_elapsed = time.perf_counter() - search_t0
                        return status, None, nvars_sat, nclauses_sat, None, cpu_sec, search_elapsed

                else:
                    from experiments import _solve_radius_cplex
                    idx, Rret, status, cpu_sec, nvars_cpx, nclauses_cpx, centers = _solve_radius_cplex(
                        inst=inst,
                        idx=mid,
                        radius=R,
                        time_limit=time_limit,
                        data=None,
                    )

                    print(
                        f"[HYBRID-CPLEX-DONE] idx={idx} R={Rret} status={status} "
                        f"cpu={cpu_sec:.6f}s nvars={nvars_cpx} nclauses={nclauses_cpx}",
                        flush=True
                    )

                    if status == "sat":
                        best_sat_idx = idx
                        best_centers = centers
                        best_cpu = cpu_sec
                        best_nvars = nvars_cpx
                        best_nclauses = nclauses_cpx
                        lo = idx + 1

                    elif status == "unsat":
                        hi = idx - 1

                    elif status in ("timeout", "error"):
                        if best_sat_idx is not None:
                            print(
                                f"[HYBRID-FALLBACK] stop on {status}; "
                                f"use best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                                flush=True
                            )
                            break

                        search_elapsed = time.perf_counter() - search_t0
                        return status, None, nvars_cpx, nclauses_cpx, None, cpu_sec, search_elapsed

        if best_sat_idx is None:
            print("[HYBRID-RESULT] infeasible (no SAT found)", flush=True)
            search_elapsed = time.perf_counter() - search_t0
            return "infeasible", None, None, None, None, None, search_elapsed

        best_radius = radii[best_sat_idx]
        print(
            f"[HYBRID-RESULT] status=OK best_idx={best_sat_idx} "
            f"best_R={best_radius} cpu={best_cpu} threshold={hybrid_threshold}",
            flush=True
        )

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
