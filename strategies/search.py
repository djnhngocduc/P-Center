import time
import traceback
import sys

from concurrent.futures import ProcessPoolExecutor, as_completed
from utils.worker import _solve_radius_worker_proc, _pool_initializer

def search_min_radius_parallel(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    radii_workers,
    seed_idx=None,
    mgr=None,
    cancel_dict=None
):
    search_t0 = time.perf_counter()

    radii = inst.radii
    nR = len(radii)
    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    decided = {}
    sat_solutions = {}
    sat_cpu_time = {}
    sat_nvars = {}
    sat_nclauses = {}

    i = seed_idx if (seed_idx is not None and 0 <= seed_idx < nR) else 0
    best_sat_idx = None

    print(
        f"[SEARCH-INIT] encoding={encoding} solver={solver_name} p={inst.p} "
        f"seed_idx={i} seed_R={radii[i]} nR={nR} workers={radii_workers}",
        flush=True
    )

    def ensure_event(idx):
        if cancel_dict is not None:
            if idx not in cancel_dict:
                cancel_dict[idx] = mgr.Event()
            return cancel_dict[idx]
        return None

    def launch_batch(ex, idxs):
        futs = {}
        print(
            f"[BATCH-SUBMIT] idxs={idxs} radii={[radii[k] for k in idxs]}",
            flush=True
        )
        for k in idxs:
            ensure_event(k)
            R = radii[k]
            print(f"[TASK-SUBMIT] idx={k} R={R}", flush=True)

            fut = ex.submit(
                _solve_radius_worker_proc,
                k,
                encoding,
                solver_name,
                R,
                time_limit,
            )
            futs[fut] = k
        return futs

    with ProcessPoolExecutor(
        max_workers=radii_workers,
        initializer=_pool_initializer,
        initargs=(cancel_dict, inst)
    ) as ex:
        while i < nR:
            j = min(nR - 1, i + radii_workers - 1)
            need = [k for k in range(i, j + 1) if k not in decided]
            if not need:
                i = j + 1
                continue

            futs = launch_batch(ex, need)
            best_sat_in_batch = None

            for fut in as_completed(futs):
                k = futs[fut]
                try:
                    idx, R, status, cpu_sec, nvars, nclauses, centers = fut.result()
                except Exception:
                    tb = traceback.format_exc()
                    print(
                        f"[FUTURE-EXCEPTION] k={k} R={radii[k]}\n{tb}",
                        file=sys.stderr,
                        flush=True,
                    )
                    idx, R = k, radii[k]
                    status, cpu_sec, nvars, nclauses, centers = (
                        "error",
                        0.0,
                        None,
                        None,
                        None,
                        None,
                    )

                print(
                    f"[TASK-DONE] idx={idx} R={R} status={status} "
                    f"cpu={cpu_sec:.6f}s",
                    flush=True
                )

                decided[idx] = status

                if status == "sat":
                    sat_solutions[idx] = centers
                    sat_cpu_time[idx] = cpu_sec
                    sat_nvars[idx] = nvars
                    sat_nclauses[idx] = nclauses

                    if (best_sat_in_batch is None) or (idx > best_sat_in_batch):
                        best_sat_in_batch = idx

                    if cancel_dict is not None:
                        for kk in need:
                            if kk < best_sat_in_batch:
                                try:
                                    cancel_dict[kk].set()
                                    print(
                                        f"[CANCEL] cancel idx={kk} "
                                        f"(dominated by SAT at idx={best_sat_in_batch})",
                                        flush=True
                                    )
                                except Exception:
                                    pass

            if best_sat_in_batch is not None:
                k = best_sat_in_batch
                print(
                    f"[REFINE] start_from_idx={k} R={radii[k]}",
                    flush=True
                )

                while k + 1 < nR:
                    nxt = k + 1

                    if nxt in decided:
                        print(
                            f"[REFINE] nxt={nxt} already decided={decided[nxt]}",
                            flush=True
                        )
                        if decided[nxt] == "unsat":
                            best_sat_idx = k
                            break
                        elif decided[nxt] == "sat":
                            k = nxt
                            continue
                        else:
                            break

                    ensure_event(nxt)
                    futs2 = launch_batch(ex, [nxt])
                    fut2 = next(iter(futs2.keys()))
                    try:
                        idx2, R2, status2, cpu_sec2, nv2, nc2, centers2 = fut2.result()
                    except Exception:
                        idx2, R2 = nxt, radii[nxt]
                        status2, cpu_sec2, nv2, nc2, centers2 = (
                            "error",
                            0.0,
                            None,
                            None,
                            None,
                            None,
                        )

                    print(
                        f"[REFINE-DONE] idx={idx2} R={R2} status={status2} "
                        f"cpu={cpu_sec2:.6f}s",
                        flush=True
                    )

                    decided[nxt] = status2

                    if status2 == "unsat":
                        best_sat_idx = k
                        break

                    elif status2 == "sat":
                        sat_solutions[nxt] = centers2
                        sat_cpu_time[nxt] = cpu_sec2
                        sat_nvars[nxt] = nv2
                        sat_nclauses[nxt] = nc2

                        if cancel_dict is not None:
                            for kk in range(0, nxt):
                                if kk in cancel_dict:
                                    try:
                                        cancel_dict[kk].set()
                                        print(
                                            f"[CANCEL] cancel idx={kk} "
                                            f"(dominated by deeper SAT at idx={nxt})",
                                            flush=True
                                        )
                                    except Exception:
                                        pass

                        k = nxt
                        continue

                    else:
                        break

                if best_sat_idx is not None:
                    break

            i = j + 1

    if best_sat_idx is None:
        if sat_solutions:
            best_sat_idx = max(sat_solutions.keys())
            print(
                f"[PARALLEL-FALLBACK] no UNSAT boundary proved; "
                f"choose best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                flush=True
            )
        else:
            print(
                "[PARALLEL-RESULT] infeasible (no SAT found at any radius)",
                flush=True
            )
            search_elapsed = time.perf_counter() - search_t0
            return "infeasible", None, None, None, None, None, search_elapsed

    best_centers = sat_solutions.get(best_sat_idx, None)
    best_sat_cpu = sat_cpu_time.get(best_sat_idx, None)
    best_nvars = sat_nvars.get(best_sat_idx, None)
    best_nclauses = sat_nclauses.get(best_sat_idx, None)

    cert_unsat_idx = best_sat_idx + 1
    certified = (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        cert_unsat_radius = radii[cert_unsat_idx]
        print(
            f"[PARALLEL-CERT] optimality boundary: "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}, "
            f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={cert_unsat_radius}",
            flush=True
        )
        final_status = "OK"
    else:
        print(
            f"[PARALLEL-FALLBACK] best SAT found but UNSAT boundary not certified; "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}",
            flush=True
        )
        final_status = "uncertified"

    print(
        f"[PARALLEL-RESULT] status={final_status} best_idx={best_sat_idx} "
        f"best_R={radii[best_sat_idx]} cpu={best_sat_cpu}",
        flush=True
    )

    search_elapsed = time.perf_counter() - search_t0
    return (
        final_status,
        radii[best_sat_idx],
        best_nvars,
        best_nclauses,
        best_centers,
        best_sat_cpu,
        search_elapsed
    )


def search_min_radius_binary(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    radii_workers,
    seed_idx=None,
    mgr=None,
    cancel_dict=None
):
    search_t0 = time.perf_counter()

    radii = inst.radii
    nR = len(radii)
    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    best_sat_idx = None
    best_centers = None
    best_sat_cpu = None
    best_nvars = None
    best_nclauses = None

    decided = {}

    lo = 0
    hi = nR - 1

    print(
        f"[BINARY-INIT] encoding={encoding} solver={solver_name} p={inst.p} "
        f"lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={nR}",
        flush=True
    )

    def ensure_event(idx):
        if cancel_dict is not None:
            if idx not in cancel_dict:
                cancel_dict[idx] = mgr.Event()
            return cancel_dict[idx]
        return None

    with ProcessPoolExecutor(
        max_workers=1,
        initializer=_pool_initializer,
        initargs=(cancel_dict, inst)
    ) as ex:
        while lo <= hi:
            mid = (lo + hi) // 2
            R = radii[mid]

            ensure_event(mid)

            print(
                f"[BINARY-STEP] lo={lo} hi={hi} mid={mid} R={R}",
                flush=True
            )

            fut = ex.submit(
                _solve_radius_worker_proc,
                mid,
                encoding,
                solver_name,
                R,
                time_limit,
            )

            try:
                idx, Rret, status, cpu_sec, nvars, nclauses, centers = fut.result()
            except Exception:
                tb = traceback.format_exc()
                print(
                    f"[BINARY-EXCEPTION] mid={mid} R={R}\n{tb}",
                    file=sys.stderr,
                    flush=True,
                )
                idx, Rret = mid, R
                status, cpu_sec, nvars, nclauses, centers = (
                    "error",
                    0.0,
                    None,
                    None,
                    None,
                    None,
                )

            decided[idx] = status

            print(
                f"[BINARY-DONE] idx={idx} R={Rret} status={status} "
                f"cpu={cpu_sec:.6f}s nvars={nvars} nclauses={nclauses}",
                flush=True
            )

            if status == "sat":
                best_sat_idx = idx
                best_centers = centers
                best_sat_cpu = cpu_sec
                best_nvars = nvars
                best_nclauses = nclauses

                new_lo = idx + 1
                print(
                    f"[BINARY-MOVE] SAT at idx={idx}, R={Rret} -> "
                    f"move right: lo {lo} -> {new_lo}, hi stays {hi}",
                    flush=True
                )
                lo = new_lo

            elif status == "unsat":
                new_hi = idx - 1
                print(
                    f"[BINARY-MOVE] UNSAT at idx={idx}, R={Rret} -> "
                    f"move left: hi {hi} -> {new_hi}, lo stays {lo}",
                    flush=True
                )
                hi = new_hi

            elif status in ("timeout", "error"):
                if best_sat_idx is not None:
                    print(
                        f"[BINARY-FALLBACK] stop on {status}; "
                        f"use best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                        flush=True
                    )
                    final_status = f"{status}_with_incumbent"
                    search_elapsed = time.perf_counter() - search_t0
                    return (
                        final_status,
                        radii[best_sat_idx],
                        best_nvars,
                        best_nclauses,
                        best_centers,
                        best_sat_cpu,
                        search_elapsed,
                    )

                search_elapsed = time.perf_counter() - search_t0
                return status, None, nvars, nclauses, None, cpu_sec, search_elapsed

    if best_sat_idx is None:
        print("[BINARY-RESULT] infeasible (no SAT found)", flush=True)
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    cert_unsat_idx = best_sat_idx + 1
    certified = (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        cert_unsat_radius = radii[cert_unsat_idx]
        print(
            f"[BINARY-CERT] optimality boundary: "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}, "
            f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={cert_unsat_radius}",
            flush=True
        )
        final_status = "OK"
    else:
        print(
            f"[BINARY-FALLBACK] best SAT found but UNSAT boundary not certified; "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}",
            flush=True
        )
        final_status = "uncertified"

    print(
        f"[BINARY-RESULT] status={final_status} best_idx={best_sat_idx} "
        f"best_R={radii[best_sat_idx]} cpu={best_sat_cpu}",
        flush=True
    )

    search_elapsed = time.perf_counter() - search_t0
    return (
        final_status,
        radii[best_sat_idx],
        best_nvars,
        best_nclauses,
        best_centers,
        best_sat_cpu,
        search_elapsed,
    )


def search_min_radius_kary(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    radii_workers,
    seed_idx=None,
    mgr=None,
    cancel_dict=None
):
    search_t0 = time.perf_counter()

    radii = inst.radii
    nR = len(radii)
    if nR == 0:
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    k = max(1, int(radii_workers))

    best_sat_idx = None
    best_unsat_idx = None

    sat_solutions = {}
    sat_cpu_time = {}
    sat_nvars = {}
    sat_nclauses = {}

    decided = {}

    lo = 0
    hi = nR - 1

    print(
        f"[KARY-INIT] encoding={encoding} solver={solver_name} p={inst.p} "
        f"lo={lo} hi={hi} R_lo={radii[lo]} R_hi={radii[hi]} nR={nR} k={k}",
        flush=True
    )

    def ensure_event(idx):
        if cancel_dict is not None:
            if idx not in cancel_dict:
                cancel_dict[idx] = mgr.Event()
            return cancel_dict[idx]
        return None

    def choose_kary_points(lo_idx, hi_idx, parts):
        if lo_idx > hi_idx:
            return []

        span = hi_idx - lo_idx + 1
        m = min(parts, span)

        if m == 1:
            return [lo_idx]

        pts = []
        for t in range(1, m + 1):
            pos = lo_idx + (t * (span + 1)) // (m + 1) - 1
            pos = max(lo_idx, min(hi_idx, pos))
            pts.append(pos)

        pts = sorted(set(pts))

        if len(pts) < m:
            for x in range(lo_idx, hi_idx + 1):
                if x not in pts:
                    pts.append(x)
                if len(pts) == m:
                    break
            pts.sort()

        return pts

    def launch_batch(ex, idxs):
        futs = {}
        print(
            f"[KARY-SUBMIT] idxs={idxs} radii={[radii[k] for k in idxs]}",
            flush=True
        )
        for idx in idxs:
            ensure_event(idx)
            fut = ex.submit(
                _solve_radius_worker_proc,
                idx,
                encoding,
                solver_name,
                radii[idx],
                time_limit,
            )
            futs[fut] = idx
        return futs

    with ProcessPoolExecutor(
        max_workers=k,
        initializer=_pool_initializer,
        initargs=(cancel_dict, inst)
    ) as ex:
        while True:
            if (
                best_sat_idx is not None
                and best_unsat_idx is not None
                and best_unsat_idx == best_sat_idx + 1
            ):
                print(
                    f"[KARY-STOP] certified optimum: "
                    f"unsat_idx={best_unsat_idx}, sat_idx={best_sat_idx}",
                    flush=True
                )
                break

            lo = 0 if best_sat_idx is None else best_sat_idx + 1
            hi = (nR - 1) if best_unsat_idx is None else best_unsat_idx - 1

            pending = [idx for idx in range(lo, hi + 1) if idx not in decided]

            print(
                f"[KARY-RANGE] unknown_range=[{lo}, {hi}] "
                f"best_unsat_idx={best_unsat_idx} "
                f"best_sat_idx={best_sat_idx} "
                f"pending={len(pending)}",
                flush=True
            )

            if not pending:
                break

            probe_lo = pending[0]
            probe_hi = pending[-1]
            probes = choose_kary_points(probe_lo, probe_hi, k)

            print(
                f"[KARY-POINTS] chosen={probes} "
                f"radii={[radii[idx] for idx in probes]}",
                flush=True
            )

            futs = launch_batch(ex, probes)

            batch_sat = []
            batch_unsat = []

            for fut in as_completed(futs):
                idx0 = futs[fut]
                try:
                    idx, R, status, cpu_sec, nvars, nclauses, centers = fut.result()
                except Exception:
                    tb = traceback.format_exc()
                    print(
                        f"[KARY-EXCEPTION] idx={idx0} R={radii[idx0]}\n{tb}",
                        file=sys.stderr,
                        flush=True,
                    )
                    idx, R = idx0, radii[idx0]
                    status, cpu_sec, nvars, nclauses, centers = (
                        "error",
                        0.0,
                        None,
                        None,
                        None,
                        None,
                    )

                print(
                    f"[KARY-DONE] idx={idx} R={R} status={status} "
                    f"cpu={cpu_sec:.6f}s nvars={nvars} nclauses={nclauses}",
                    flush=True
                )

                decided[idx] = status

                if status == "sat":
                    batch_sat.append(idx)
                    sat_solutions[idx] = centers
                    sat_cpu_time[idx] = cpu_sec
                    sat_nvars[idx] = nvars
                    sat_nclauses[idx] = nclauses

                elif status == "unsat":
                    batch_unsat.append(idx)

                elif status in ("timeout", "error"):
                    pass

            if batch_sat:
                sat_max = max(batch_sat)
                old_best_sat = best_sat_idx
                if best_sat_idx is None or sat_max > best_sat_idx:
                    best_sat_idx = sat_max
                print(
                    f"[KARY-UPDATE] SAT batch={batch_sat}, "
                    f"best_sat_idx: {old_best_sat} -> {best_sat_idx}",
                    flush=True
                )

                if cancel_dict is not None:
                    for kk in list(cancel_dict.keys()):
                        if kk < best_sat_idx:
                            try:
                                cancel_dict[kk].set()
                            except Exception:
                                pass

            if batch_unsat:
                unsat_min = min(batch_unsat)
                old_best_unsat = best_unsat_idx
                if best_unsat_idx is None or unsat_min < best_unsat_idx:
                    best_unsat_idx = unsat_min
                print(
                    f"[KARY-UPDATE] UNSAT batch={batch_unsat}, "
                    f"best_unsat_idx: {old_best_unsat} -> {best_unsat_idx}",
                    flush=True
                )

                if cancel_dict is not None:
                    for kk in list(cancel_dict.keys()):
                        if kk > best_unsat_idx:
                            try:
                                cancel_dict[kk].set()
                            except Exception:
                                pass

    if best_sat_idx is None:
        print("[KARY-RESULT] infeasible (no SAT found)", flush=True)
        search_elapsed = time.perf_counter() - search_t0
        return "infeasible", None, None, None, None, None, search_elapsed

    best_centers = sat_solutions.get(best_sat_idx, None)
    best_sat_cpu = sat_cpu_time.get(best_sat_idx, None)
    best_nvars = sat_nvars.get(best_sat_idx, None)
    best_nclauses = sat_nclauses.get(best_sat_idx, None)

    cert_unsat_idx = best_sat_idx + 1
    certified = (
        cert_unsat_idx < nR and decided.get(cert_unsat_idx) == "unsat"
    )

    if certified:
        cert_unsat_radius = radii[cert_unsat_idx]
        print(
            f"[KARY-CERT] optimality boundary: "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}, "
            f"next_unsat_idx={cert_unsat_idx}, next_unsat_R={cert_unsat_radius}",
            flush=True
        )
        final_status = "OK"
    else:
        print(
            f"[KARY-FALLBACK] best SAT found but UNSAT boundary not certified; "
            f"best_sat_idx={best_sat_idx}, best_R={radii[best_sat_idx]}",
            flush=True
        )
        final_status = "uncertified"

    print(
        f"[KARY-RESULT] status={final_status} best_idx={best_sat_idx} "
        f"best_R={radii[best_sat_idx]} cpu={best_sat_cpu}",
        flush=True
    )

    search_elapsed = time.perf_counter() - search_t0
    return (
        final_status,
        radii[best_sat_idx],
        best_nvars,
        best_nclauses,
        best_centers,
        best_sat_cpu,
        search_elapsed,
    )