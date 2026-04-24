import time
import threading
import json
from collections import defaultdict
import statistics as _stats
from solvers.backend import EXTERNAL_SOLVERS, cpu_self_seconds, solve_radius_gurobi, solve_radius_cplex
from utils.threshold import remaining_iters_from_bounds, solver_wall_cost
from pysat.solvers import Solver
from utils.io import write_csv_rows

try:
    import cplex
except Exception:
    cplex = None

try:
    import gurobipy as gp
except Exception:
    gp = None


def _check_mip_backend_available(mip_backend: str):
    if mip_backend == "cplex_mip":
        if cplex is None:
            raise ImportError("threshold_profile with mip_backend='cplex_mip' needs CPLEX installed.")
    elif mip_backend == "gurobi_mip":
        if gp is None:
            raise ImportError("threshold_profile with mip_backend='gurobi_mip' needs Gurobi installed.")
    else:
        raise ValueError(f"Unsupported mip_backend: {mip_backend}")


def _solve_radius_mip(inst, idx, radius, time_limit, mip_backend):
    if mip_backend == "cplex_mip":
        return solve_radius_cplex(
            inst=inst,
            idx=idx,
            radius=radius,
            time_limit=time_limit,
            data=None,
        )
    elif mip_backend == "gurobi_mip":
        return solve_radius_gurobi(
            inst=inst,
            idx=idx,
            radius=radius,
            time_limit=time_limit,
            data=None,
        )
    else:
        raise ValueError(f"Unsupported mip_backend: {mip_backend}")


def profile_threshold_sat_vs_mip(
    inst,
    encoding,
    sat_solver_name,
    time_limit,
    *,
    mip_backend="cplex_mip",
    instance_name,
    detail_out=None,
    summary_out=None,
):
    if sat_solver_name in EXTERNAL_SOLVERS:
        raise ValueError(
            f"threshold_profile requires an INTERNAL PySAT solver for incremental SAT, "
            f"but got external solver '{sat_solver_name}'."
        )

    _check_mip_backend_available(mip_backend)

    search_t0 = time.perf_counter()

    base_cnf, info = inst._build_incremental_base_cnf(encoding=encoding)
    radii = info["radii"]
    y_vars = info["y"]

    nR = len(radii)
    if nR == 0:
        return {
            "status": "infeasible",
            "best_radius": None,
            "detail_rows": [],
            "summary_rows": [],
            "best_threshold": None,
            "best_threshold_total_wall": None,
            "search_time": time.perf_counter() - search_t0,
        }

    lo = 0
    hi = nR - 1

    best_sat_idx = None
    best_centers = None

    next_var = base_cnf.nv + 1
    tested = {}

    detail_rows = []

    print(
        f"[THRESHOLD-INIT] encoding={encoding} sat_solver={sat_solver_name} mip_backend={mip_backend} "
        f"p={inst.p} n={inst.n} lo={lo} hi={hi} "
        f"R_lo={radii[lo]} R_hi={radii[hi]} nR={nR} "
        f"base_clauses={len(base_cnf.clauses)} base_vars={base_cnf.nv}",
        flush=True
    )

    with Solver(name=sat_solver_name, bootstrap_with=base_cnf.clauses) as solver:
        step_no = 0

        while lo <= hi:
            step_no += 1
            mid = (lo + hi) // 2
            R = radii[mid]
            remaining_iters = remaining_iters_from_bounds(lo, hi)

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
                    f"[THRESHOLD-ADD] step={step_no} idx={mid} R={R} "
                    f"selector={alpha} added_clauses={len(step_clauses)}",
                    flush=True
                )
            else:
                alpha = tested[mid]

            print(
                f"[THRESHOLD-STEP] step={step_no} lo={lo} hi={hi} mid={mid} "
                f"R={R} remaining_iters={remaining_iters}",
                flush=True
            )

            sat_timer = None
            sat_cpu0 = None
            sat_cpu1 = None
            sat_wall0 = time.perf_counter()
            sat_status = "error"
            sat_cpu = 0.0
            sat_centers = None

            try:
                if time_limit and time_limit > 0 and hasattr(solver, "interrupt"):
                    sat_timer = threading.Timer(time_limit, solver.interrupt)
                    sat_timer.start()
                    sat_cpu0 = cpu_self_seconds()
                    try:
                        sat_res = solver.solve_limited(
                            assumptions=[alpha],
                            expect_interrupt=True
                        )
                    except NotImplementedError:
                        sat_res = solver.solve(assumptions=[alpha])
                    sat_cpu1 = cpu_self_seconds()
                else:
                    sat_cpu0 = cpu_self_seconds()
                    sat_res = solver.solve(assumptions=[alpha])
                    sat_cpu1 = cpu_self_seconds()
            finally:
                if sat_timer:
                    sat_timer.cancel()

            sat_wall = time.perf_counter() - sat_wall0
            sat_cpu = (sat_cpu1 - sat_cpu0) if (sat_cpu0 is not None and sat_cpu1 is not None) else 0.0

            if sat_res is None:
                sat_status = "timeout"
            elif sat_res:
                model = solver.get_model() or []
                model_set = set(model)
                sat_centers = sorted([j for j, v in enumerate(y_vars) if v in model_set])
                sat_status = "sat"
            else:
                sat_status = "unsat"

            mip_wall0 = time.perf_counter()
            m_idx, m_R, mip_status, mip_cpu, mip_nvars, mip_nclauses, mip_centers = _solve_radius_mip(
                inst=inst,
                idx=mid,
                radius=R,
                time_limit=time_limit,
                mip_backend=mip_backend,
            )
            mip_wall = time.perf_counter() - mip_wall0

            status_mismatch = False
            if sat_status in ("sat", "unsat") and mip_status in ("sat", "unsat"):
                if sat_status != mip_status:
                    status_mismatch = True
                    print(
                        f"[THRESHOLD-WARN] mismatch at idx={mid} R={R}: "
                        f"SAT={sat_status}, {mip_backend}={mip_status}",
                        flush=True
                    )

            winner = "tie"
            if sat_status in ("sat", "unsat") and mip_status in ("sat", "unsat"):
                if sat_wall < mip_wall - 1e-12:
                    winner = "sat"
                elif mip_wall < sat_wall - 1e-12:
                    winner = mip_backend
            elif sat_status in ("sat", "unsat"):
                winner = "sat"
            elif mip_status in ("sat", "unsat"):
                winner = mip_backend

            row = {
                "instance": instance_name,
                "mip_backend": mip_backend,
                "step_no": step_no,
                "lo_before": lo,
                "hi_before": hi,
                "mid_idx": mid,
                "radius": R,
                "remaining_iters": remaining_iters,
                "sat_status": sat_status,
                "sat_wall": sat_wall,
                "sat_cpu": sat_cpu,
                "sat_centers": json.dumps(sat_centers if sat_centers is not None else []),
                "mip_status": mip_status,
                "mip_wall": mip_wall,
                "mip_cpu": mip_cpu,
                "mip_centers": json.dumps(mip_centers if mip_centers is not None else []),
                "winner": winner,
                "status_mismatch": int(status_mismatch),
            }
            detail_rows.append(row)

            print(
                f"[THRESHOLD-DONE] step={step_no} idx={mid} R={R} "
                f"remaining_iters={remaining_iters} "
                f"SAT=({sat_status}, wall={sat_wall:.6f}s, cpu={sat_cpu:.6f}s) "
                f"{mip_backend}=({mip_status}, wall={mip_wall:.6f}s, cpu={mip_cpu:.6f}s) "
                f"winner={winner}",
                flush=True
            )

            decision_status = sat_status
            if decision_status not in ("sat", "unsat"):
                decision_status = mip_status

            if decision_status == "sat":
                best_sat_idx = mid
                best_centers = sat_centers if sat_centers is not None else mip_centers
                lo = mid + 1
            elif decision_status == "unsat":
                hi = mid - 1
            elif decision_status == "timeout":
                return {
                    "status": "timeout",
                    "best_radius": None,
                    "detail_rows": detail_rows,
                    "summary_rows": [],
                    "best_threshold": None,
                    "best_threshold_total_wall": None,
                    "search_time": time.perf_counter() - search_t0,
                }
            else:
                return {
                    "status": "error",
                    "best_radius": None,
                    "detail_rows": detail_rows,
                    "summary_rows": [],
                    "best_threshold": None,
                    "best_threshold_total_wall": None,
                    "search_time": time.perf_counter() - search_t0,
                }

    grouped = defaultdict(list)
    for r in detail_rows:
        grouped[int(r["remaining_iters"])].append(r)

    summary_rows = []

    max_k = max(grouped.keys()) if grouped else 0
    for k in sorted(grouped.keys(), reverse=True):
        rows_k = grouped[k]
        sat_wins = sum(1 for r in rows_k if r["winner"] == "sat")
        mip_wins = sum(1 for r in rows_k if r["winner"] == mip_backend)
        ties = sum(1 for r in rows_k if r["winner"] == "tie")
        mismatches = sum(int(r["status_mismatch"]) for r in rows_k)

        sat_wall_vals = [float(r["sat_wall"]) for r in rows_k]
        mip_wall_vals = [float(r["mip_wall"]) for r in rows_k]
        sat_cpu_vals = [float(r["sat_cpu"]) for r in rows_k]
        mip_cpu_vals = [float(r["mip_cpu"]) for r in rows_k]

        summary_rows.append({
            "instance": instance_name,
            "mip_backend": mip_backend,
            "remaining_iters": k,
            "count": len(rows_k),
            "sat_wins": sat_wins,
            "mip_wins": mip_wins,
            "ties": ties,
            "sat_win_rate": (sat_wins / len(rows_k)) if rows_k else 0.0,
            "mip_win_rate": (mip_wins / len(rows_k)) if rows_k else 0.0,
            "sat_wall_mean": _stats.mean(sat_wall_vals) if sat_wall_vals else None,
            "mip_wall_mean": _stats.mean(mip_wall_vals) if mip_wall_vals else None,
            "sat_wall_median": _stats.median(sat_wall_vals) if sat_wall_vals else None,
            "mip_wall_median": _stats.median(mip_wall_vals) if mip_wall_vals else None,
            "sat_cpu_mean": _stats.mean(sat_cpu_vals) if sat_cpu_vals else None,
            "mip_cpu_mean": _stats.mean(mip_cpu_vals) if mip_cpu_vals else None,
            "mismatches": mismatches,
        })

    threshold_rows = []
    best_threshold = None
    best_threshold_total_wall = None

    for T in range(0, max_k + 1):
        total_cost = 0.0
        sat_steps = 0
        mip_steps = 0

        for r in detail_rows:
            k = int(r["remaining_iters"])
            if k > T:
                total_cost += solver_wall_cost(
                    r["sat_status"], r["sat_wall"], time_limit
                )
                sat_steps += 1
            else:
                total_cost += solver_wall_cost(
                    r["mip_status"], r["mip_wall"], time_limit
                )
                mip_steps += 1

        threshold_rows.append({
            "instance": instance_name,
            "mip_backend": mip_backend,
            "threshold": T,
            "policy": f"use SAT if remaining_iters > T else {mip_backend}",
            "total_hybrid_wall": total_cost,
            "sat_steps": sat_steps,
            "mip_steps": mip_steps,
        })

        if best_threshold_total_wall is None or total_cost < best_threshold_total_wall:
            best_threshold_total_wall = total_cost
            best_threshold = T

    summary_export_rows = []

    for r in summary_rows:
        rr = dict(r)
        rr["row_type"] = "remaining_iters_summary"
        summary_export_rows.append(rr)

    for r in threshold_rows:
        rr = dict(r)
        rr["row_type"] = "threshold_simulation"
        summary_export_rows.append(rr)

    if detail_out:
        write_csv_rows(
            detail_out,
            detail_rows,
            fieldnames=[
                "instance",
                "mip_backend",
                "step_no",
                "lo_before",
                "hi_before",
                "mid_idx",
                "radius",
                "remaining_iters",
                "sat_status",
                "sat_wall",
                "sat_cpu",
                "sat_centers",
                "mip_status",
                "mip_wall",
                "mip_cpu",
                "mip_centers",
                "winner",
                "status_mismatch",
            ],
            append=True,
        )

    if summary_out:
        write_csv_rows(
            summary_out,
            summary_export_rows,
            fieldnames=[
                "instance",
                "mip_backend",
                "row_type",
                "remaining_iters",
                "count",
                "sat_wins",
                "mip_wins",
                "ties",
                "sat_win_rate",
                "mip_win_rate",
                "sat_wall_mean",
                "mip_wall_mean",
                "sat_wall_median",
                "mip_wall_median",
                "sat_cpu_mean",
                "mip_cpu_mean",
                "mismatches",
                "threshold",
                "policy",
                "total_hybrid_wall",
                "sat_steps",
                "mip_steps",
            ],
            append=True,
        )

    best_radius = radii[best_sat_idx] if best_sat_idx is not None else None

    print(
        f"[THRESHOLD-RESULT] best_radius={best_radius} mip_backend={mip_backend} "
        f"best_threshold={best_threshold} best_threshold_total_wall={best_threshold_total_wall}",
        flush=True
    )

    return {
        "status": "OK" if best_sat_idx is not None else "infeasible",
        "best_radius": best_radius,
        "detail_rows": detail_rows,
        "summary_rows": summary_export_rows,
        "best_threshold": best_threshold,
        "best_threshold_total_wall": best_threshold_total_wall,
        "search_time": time.perf_counter() - search_t0,
    }
