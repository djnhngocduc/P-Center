import sys
import time
import traceback

from utils.worker import _pool_initializer, _solve_radius_worker_proc


def search_min_radius_binary(
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

    def _timeout_result():
        search_elapsed = time.perf_counter() - search_t0
        if best_sat_idx is not None:
            return (
                "timeout_with_incumbent",
                radii[best_sat_idx],
                best_nvars,
                best_nclauses,
                best_centers,
                search_elapsed,
            )
        return "timeout", None, None, None, None, search_elapsed

    radii = inst.radii
    nR = len(radii)
    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, search_elapsed

    best_sat_idx = None
    best_centers = None
    best_nvars = None
    best_nclauses = None

    decided = {}

    lo = 0
    hi = nR - 1

    print(
        f"[BINARY-INIT] encoding={encoding} solver={solver_name} p={inst.p} "
        f"lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={nR}",
        flush=True,
    )

    _pool_initializer(None, inst)

    while lo <= hi:
        step_time_limit = _remaining_time()
        if deadline is not None and step_time_limit <= 0:
            print("[BINARY-TIMEOUT] instance-level time budget exhausted", flush=True)
            return _timeout_result()

        mid = (lo + hi) // 2
        R = radii[mid]

        print(
            f"[BINARY-STEP] lo={lo} hi={hi} mid={mid} R={R}",
            flush=True,
        )

        try:
            idx, Rret, status, nvars, nclauses, centers = _solve_radius_worker_proc(
                mid,
                encoding,
                solver_name,
                R,
                step_time_limit,
            )
        except Exception:
            tb = traceback.format_exc()
            print(
                f"[BINARY-EXCEPTION] mid={mid} R={R}\n{tb}",
                file=sys.stderr,
                flush=True,
            )
            idx, Rret = mid, R
            status, nvars, nclauses, centers = (
                "error",
                None,
                None,
                None,
            )

        decided[idx] = status

        print(
            f"[BINARY-DONE] idx={idx} R={Rret} status={status} "
            f"nvars={nvars} nclauses={nclauses}",
            flush=True,
        )

        if status == "sat":
            best_sat_idx = idx
            best_centers = centers
            best_nvars = nvars
            best_nclauses = nclauses

            new_lo = idx + 1
            print(
                f"[BINARY-MOVE] SAT at idx={idx}, R={Rret} -> "
                f"move right: lo {lo} -> {new_lo}, hi stays {hi}",
                flush=True,
            )
            lo = new_lo

        elif status == "unsat":
            new_hi = idx - 1
            print(
                f"[BINARY-MOVE] UNSAT at idx={idx}, R={Rret} -> "
                f"move left: hi {hi} -> {new_hi}, lo stays {lo}",
                flush=True,
            )
            hi = new_hi

        elif status in ("timeout", "error"):
            if best_sat_idx is not None:
                print(
                    f"[BINARY-FALLBACK] stop on {status}; "
                    f"use best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                    flush=True,
                )
                final_status = f"{status}_with_incumbent"
                search_elapsed = time.perf_counter() - search_t0
                return (
                    final_status,
                    radii[best_sat_idx],
                    best_nvars,
                    best_nclauses,
                    best_centers,
                    search_elapsed,
                )

            search_elapsed = time.perf_counter() - search_t0
            return status, None, nvars, nclauses, None, search_elapsed

    if best_sat_idx is None:
        print("[BINARY-RESULT] infeasible (no SAT found)", flush=True)
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, search_elapsed

    cert_unsat_idx = best_sat_idx + 1
    is_smallest_candidate = best_sat_idx == nR - 1
    certified = is_smallest_candidate or (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        if is_smallest_candidate:
            print(
                f"[BINARY-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]} "
                f"is the smallest candidate radius",
                flush=True,
            )
        else:
            cert_unsat_radius = radii[cert_unsat_idx]
            print(
                f"[BINARY-CERT] optimality boundary: "
                f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}, "
                f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={cert_unsat_radius}",
                flush=True,
            )
        final_status = "OK"
    else:
        print(
            f"[BINARY-FALLBACK] best SAT found but UNSAT boundary not certified; "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}",
            flush=True,
        )
        final_status = "uncertified"

    print(
        f"[BINARY-RESULT] status={final_status} best_idx={best_sat_idx} "
        f"best_R={radii[best_sat_idx]}",
        flush=True,
    )

    search_elapsed = time.perf_counter() - search_t0
    return (
        final_status,
        radii[best_sat_idx],
        best_nvars,
        best_nclauses,
        best_centers,
        search_elapsed,
    )
