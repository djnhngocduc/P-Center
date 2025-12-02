import argparse
import json
import time
import threading
import heapq
import os
import sys
import subprocess
import tempfile
from typing import List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from encoder import PCenterSAT
import multiprocessing as mp
from collections import defaultdict
import statistics as _stats
import csv

_CANCEL_SHARED = None
INF = 10 ** 12

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

EXTERNAL_SOLVERS = {
    "maplecm": (
        os.path.join(BASE_DIR, "solvers", "Maple_CM", "bin", "starexec_run_default"),
        []
    ),

    "maplechrono": (
        os.path.join(BASE_DIR, "solvers", "MapleLCMDistChronoBT", "bin", "starexec_run_default"),
        []
    ),

    "sparrow2riss": (
        os.path.join(BASE_DIR, "solvers", "Sparrow2Riss-2018", "bin", "starexec_run_NOLIMIT"),
        []
    )
}


def _pool_initializer(cancel_proxy):
    # Mỗi worker nhận proxy dict để đọc Event theo idx
    global _CANCEL_SHARED
    _CANCEL_SHARED = cancel_proxy

    # In PID worker 1 lần để kiểm chứng ~8 PID cố định
    print(f"[POOL] worker PID={os.getpid()}", file=sys.stderr, flush=True)


# ---------------- IO helpers ----------------

def load_instance_data(inst_file: str):
    with open(inst_file, "r", encoding="utf-8", errors="ignore") as f:
        return json.load(f)


def load_instance(inst_desc):
    if "tsplib" in inst_desc:
        coords = load_tsplib_coords(inst_desc["tsplib"])
        inst = PCenterSAT.from_coordinates(coords, inst_desc["p"])

        inst.radii = sorted(set(inst.radii), reverse=True)

        sr = float(inst_desc["seed_radius"])

        seed_idx = next((i for i, r in enumerate(inst.radii) if r <= sr), None)

        if seed_idx is None:
            seed_idx = len(inst.radii) - 1

    elif "orlib" in inst_desc:
        t0 = time.time()
        n_graph, p, edges = load_orlib_edge_list(inst_desc["orlib"])
        print(f"[LOAD] {inst_desc['name']}: n={n_graph}, m={len(edges)}, p={p}", flush=True)

        t1 = time.time()
        D = compute_apsp_dijkstra(n_graph, edges)
        print(
            f"[APSP] {inst_desc['name']}: done in {time.time() - t1:.2f}s "
            f"(total load {time.time() - t0:.2f}s)",
            flush=True
        )
        inst = PCenterSAT.from_distance_matrix(D, p)

        radii_set = set()
        for i in range(n_graph):
            for j in range(i + 1, n_graph):
                dij = D[i][j]
                if dij != INF and dij > 0:
                    radii_set.add(dij)

        inst.radii = sorted(radii_set, reverse=True)

        sr = int(inst_desc["seed_radius"])
        seed_idx = next((i for i, r in enumerate(inst.radii) if r <= sr), None)

        if seed_idx is None:
            seed_idx = len(inst.radii) - 1
    else:
        raise ValueError("Instance description must contain 'tsplib' or 'orlib'")

    if inst.radii:
        seed_R = inst.radii[seed_idx]
        print(
            f"[SEED] {inst_desc['name']}: "
            f"seed_idx={seed_idx}, seed_radius={seed_R}, "
            f"radii_len={len(inst.radii)}, p={inst.p}, n={inst.n}",
            flush=True
        )
    else:
        print(f"[SEED] {inst_desc['name']}: radii empty?!", flush=True)

    return inst, seed_idx


def load_tsplib_coords(path: str):
    coords = []
    in_section = False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            u = s.upper()
            if u.startswith("NODE_COORD_SECTION"):
                in_section = True
                continue
            if u.startswith("EOF"):
                break
            if not in_section:
                continue
            parts = s.split()
            if len(parts) >= 3:
                x = float(parts[1])
                y = float(parts[2])
                coords.append((x, y))
    return coords


def load_orlib_edge_list(path: str):
    edges = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().split()

        n = int(header[0])
        m = int(header[1])
        p = int(header[2])

        for _ in range(m):
            line = f.readline()
            if not line:
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            # Giảm chỉ số từ 1 xuống 0
            u, v, w = int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2])
            edges.append((u, v, w))

    return n, p, edges


def compute_apsp_dijkstra(n: int, edges: List[Tuple[int, int, float]]) -> List[List[int]]:
    c = [[INF] * n for _ in range(n)]
    for i in range(n):
        c[i][i] = 0

    for u, v, w in edges:
        ww = int(w)
        if u == v:
            c[u][v] = 0
        else:
            c[u][v] = ww
            c[v][u] = ww

    adj = [[] for _ in range(n)]
    for u in range(n):
        row = c[u]
        for v in range(n):
            w = row[v]
            if u != v and w < INF:
                adj[u].append((v, w))

    D = [[INF] * n for _ in range(n)]
    for s in range(n):
        dist = [INF] * n
        dist[s] = 0
        pq = [(0, s)]
        while pq:
            d, u = heapq.heappop(pq)
            if d > dist[u]:
                continue
            for v, w in adj[u]:
                nd = d + w
                if nd < dist[v]:
                    dist[v] = nd
                    heapq.heappush(pq, (nd, v))
        D[s] = [int(x) if x < INF else INF for x in dist]

    return D


def search_min_radius_parallel(
    inst,
    encoding,
    solver_name,
    time_limit,
    *,
    radii_workers=8,
    seed_idx=None,
    mgr=None,
    cancel_dict=None
):
    """
    Dùng pool cố định (ProcessPoolExecutor) với max_workers = radii_workers.
    Ngắt mềm job bán kính lớn hơn bằng Event + solver.interrupt().
    """
    radii = inst.radii
    nR = len(radii)
    start_wall = time.time()
    if nR == 0:
        return "infeasible", None, 0.0, None, None, None, None

    decided = {}
    best_nvars = None
    best_nclauses = None
    sat_solutions = {}
    sat_cpu_time = {}

    # seed mặc định = 0 (bán kính lớn nhất)
    i = seed_idx if (seed_idx is not None and 0 <= seed_idx < nR) else 0
    covered = set()
    best_sat_idx = None

    print(
        f"[SEARCH-INIT] encoding={encoding} solver={solver_name} p={inst.p} "
        f"seed_idx={i} seed_R={radii[i]} nR={nR} workers={radii_workers}",
        flush=True
    )

    # helper: tạo Event cho một idx nếu chưa có
    def ensure_event(idx):
        if cancel_dict is not None:
            if idx not in cancel_dict:
                cancel_dict[idx] = mgr.Event()  # Event proxy
            return cancel_dict[idx]
        return None

    # submit batch lên pool
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
                inst,
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
        initargs=(cancel_dict,)
    ) as ex:
        while i < nR:
            j = min(nR - 1, i + radii_workers - 1)
            need = [k for k in range(i, j + 1) if k not in covered]
            if not need:
                i = j + 1
                continue

            futs = launch_batch(ex, need)
            best_sat_in_batch = None

            for fut in as_completed(futs):
                k = futs[fut]
                try:
                    idx, R, status, t_sec, nvars, nclauses, centers = fut.result()
                except Exception:
                    idx, R = k, radii[k]
                    status, t_sec, nvars, nclauses, centers = (
                        "timeout",
                        0.0,
                        None,
                        None,
                        None,
                    )
                print(
                    f"[TASK-DONE] idx={idx} R={R} status={status} "
                    f"t_cpu={(t_sec if t_sec is not None else 0.0):.6f}s",
                    flush=True
                )

                decided[idx] = status

                if nvars is not None and nclauses is not None:
                    best_nvars = nvars if best_nvars is None else min(best_nvars, nvars)
                    best_nclauses = (
                        nclauses
                        if best_nclauses is None
                        else min(best_nclauses, nclauses)
                    )

                if status == "sat":
                    sat_solutions[idx] = centers
                    sat_cpu_time[idx] = t_sec
                    if (best_sat_in_batch is None) or (idx > best_sat_in_batch):
                        best_sat_in_batch = idx
                    # NGẮT MỀM mọi job có idx < idx_sat (bán kính lớn hơn)
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

            covered.update(need)

            # Co dần nếu có SAT trong batch
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

                    # chạy đơn lẻ bán kính tiếp theo
                    ensure_event(nxt)
                    futs2 = launch_batch(ex, [nxt])
                    # không cần as_completed loop phức tạp vì chỉ 1 future
                    fut2 = next(iter(futs2.keys()))
                    try:
                        idx2, R2, status2, t2, nv2, nc2, centers2 = fut2.result()
                    except Exception:
                        status2, t2, nv2, nc2, centers2 = (
                            "timeout",
                            0.0,
                            None,
                            None,
                            None,
                        )

                    print(
                        f"[REFINE-DONE] idx={idx2} R={R2} status={status2} t_cpu={t2}",
                        flush=True
                    )

                    decided[nxt] = status2
                    if nv2 is not None and nc2 is not None:
                        best_nvars = (
                            nv2 if best_nvars is None else min(best_nvars, nv2)
                        )
                        best_nclauses = (
                            nc2
                            if best_nclauses is None
                            else min(best_nclauses, nc2)
                        )

                    if status2 == "unsat":
                        best_sat_idx = k
                        break
                    elif status2 == "sat":
                        sat_solutions[nxt] = centers2
                        sat_cpu_time[nxt] = t2
                        # NGẮT MỀM mọi job có idx < nxt
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
                        # timeout -> dừng co để tránh vòng vô hạn
                        break

                if best_sat_idx is not None:
                    break

            i = j + 1

    elapsed = time.time() - start_wall

    if best_sat_idx is None:
        if sat_solutions:
            best_sat_idx = max(sat_solutions.keys())
            print(
                f"[FALLBACK] no UNSAT boundary proved; "
                f"choose best_sat_idx={best_sat_idx} R={radii[best_sat_idx]}",
                flush=True
            )
        else:
            print(
                "[RESULT] infeasible (no SAT found at any radius)",
                flush=True
            )
            return "infeasible", None, elapsed, best_nvars, best_nclauses, None, None

    best_centers = sat_solutions.get(best_sat_idx, None)
    best_sat_cpu = sat_cpu_time.get(best_sat_idx, None)

    print(
        f"[RESULT] status=OK best_idx={best_sat_idx} best_R={radii[best_sat_idx]} "
        f"elapsed_wall={elapsed:.6f}s cpu_at_best={best_sat_cpu}",
        flush=True
    )

    return (
        "OK",
        radii[best_sat_idx],
        elapsed,
        best_nvars,
        best_nclauses,
        best_centers,
        best_sat_cpu,
    )


def _write_dimacs(cnf, path):
    with open(path, "w", encoding="utf-8") as f:
        nvars = cnf.nv
        nclauses = len(cnf.clauses)
        f.write(f"p cnf {nvars} {nclauses}\n")
        for cl in cnf.clauses:
            if not cl:
                # mệnh đề rỗng
                f.write("0\n")
            else:
                f.write(" ".join(str(l) for l in cl) + " 0\n")


def _run_external_solver(solver_name, cnf, time_limit, cancel_ev=None):
    """
    Chạy solver ngoài (Maple_CM, MapleLCMDistChronoBT, Sparrow2Riss-2018).
    Trả về: (status, solver_time, model_lits | None)

    status: "sat" | "unsat" | "timeout" | "error"
    solver_time: thời gian CPU xấp xỉ (wall-clock tính trong Python)
    model_lits: list[int] model (các literal) nếu SAT, ngược lại None
    """
    bin_path, extra_args = EXTERNAL_SOLVERS[solver_name]

    with tempfile.TemporaryDirectory() as tmpdir:
        cnf_path = os.path.join(tmpdir, "formula.cnf")
        _write_dimacs(cnf, cnf_path)

        cmd = [bin_path] + extra_args + [cnf_path]

        solver_dir = os.path.dirname(bin_path)

        t0 = time.perf_counter()
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=solver_dir,
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
            t1 = time.perf_counter()
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
            return "timeout", (time_limit if time_limit else 0.0), None

        solver_time = t1 - t0

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

        print(
            f"[SOLVER-RET] solver={solver_name} rc={proc.returncode} "
            f"raw_status={status} time={solver_time:.6f}s",
            flush=True
        )

        if proc.returncode not in (0, 10, 20):
            print(f"[SOLVER-STDERR]\n{stderr}", flush=True)
            print(f"[SOLVER-STDOUT]\n{stdout}", flush=True)

        if status is None:
            if model_lits:
                status = "sat"
            else:
                status = "error" if proc.returncode not in (10, 20, 0) else "unsat"

        if status == "sat" and not model_lits:
            print(f"[SOLVER-WARN] {solver_name} reported SAT but no model.", flush=True)
            return "error", solver_time, None

        return status, solver_time, (model_lits if status == "sat" else None)


def _solve_radius_worker_proc(idx, inst, encoding, solver_name, radius, time_limit):
    """
    Chạy trong process pool (executor). Trả về:
    (idx, radius, status, elapsed, nvars, nclauses, centers)
    """
    pid = os.getpid()
    print(
        f"[WORKER-START] pid={pid} idx={idx} R={radius} "
        f"enc={encoding} solver={solver_name} limit={time_limit}s",
        flush=True
    )

    try:
        cnf, varmap = inst._encode_cnf(radius, encoding)

        print(
            f"[ENCODE] idx={idx} R={radius} |Nc|={len(varmap.get('Nc', []))} "
            f"|Nd|={len(varmap.get('Nd', []))} "
            f"candidates={len(varmap.get('candidates', []))} "
            f"bound={varmap.get('bound')} "
            f"clauses={0 if cnf is None else len(cnf.clauses)} "
            f"vars={0 if cnf is None else cnf.nv}",
            flush=True
        )

        if cnf is None:
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> UNSAT(by reduction)",
                flush=True
            )
            return idx, radius, "unsat", 0.0, None, None, None

        if solver_name in EXTERNAL_SOLVERS:
            # cancel_ev lấy từ shared dict
            cancel_ev = None
            try:
                if _CANCEL_SHARED is not None:
                    cancel_ev = _CANCEL_SHARED.get(idx)
            except Exception:
                cancel_ev = None

            status, solver_time, model = _run_external_solver(
                solver_name=solver_name,
                cnf=cnf,
                time_limit=time_limit,
                cancel_ev=cancel_ev,
            )

            nvars = cnf.nv
            nclauses = len(cnf.clauses)

            if status in ("timeout", "error"):
                print(
                    f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                    f"-> {status.upper()} cpu={solver_time:.6f}s",
                    flush=True
                )
                return idx, radius, status, solver_time, nvars, nclauses, None

            if status == "unsat":
                print(
                    f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                    f"-> UNSAT cpu={solver_time:.6f}s",
                    flush=True
                )
                return idx, radius, "unsat", solver_time, nvars, nclauses, None

            # SAT
            model_set = set(model or [])
            y_vars = varmap.get("y", [])
            Nc = set(varmap.get("Nc", []))
            candidates = set(varmap.get("candidates", []))
            chosen = {
                j
                for j, v in enumerate(y_vars)
                if (v in model_set) and (j in candidates)
            }
            centers = sorted(chosen | Nc)
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> SAT cpu={solver_time:.6f}s centers={centers}",
                flush=True
            )
            return idx, radius, "sat", solver_time, nvars, nclauses, centers

    except Exception as e:
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius} -> EXCEPTION {e}",
            flush=True
        )
        return idx, radius, "timeout", 0.0, None, None, None


# ---------------- CLI & Runner ----------------

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--instances",
        type=str,
        required=True,
        help="Path to instances.json",
    )
    ap.add_argument(
        "--encodings",
        nargs="+",
        default=["pysat_sc", "pypb_sc", "nsc", "pb_bdd"],
        help="Encodings to test: pysat_sc pypb_sc nsc pb_bdd",
    )
    ap.add_argument(
        "--solvers",
        nargs="+",
        default=["maplecm", "maplechrono", "sparrow2riss"],
        help="SAT solvers to use",
    )
    ap.add_argument(
        "--time-limit",
        type=int,
        default=14400,
        help="Time limit per radius solve (seconds)",
    )
    ap.add_argument(
        "--radii-workers",
        type=int,
        default=8,
        help="Parallel radii workers",
    )
    ap.add_argument(
        "--out",
        type=str,
        default="results.csv",
        help="CSV output path",
    )
    ap.add_argument(
        "--strategy",
        type=str,
        choices=["parallel"],
        default="parallel",
        help="Radius search strategy (paper compares parallel only)",
    )
    return ap.parse_args()


def run_experiment(
    inst_desc,
    encodings,
    solvers,
    time_limit,
    radii_workers,
    strategy,
    *,
    mgr,
    cancel_dict
):
    results = []
    inst, seed_idx = load_instance(inst_desc)
    print(
        f"[RUN-EXPERIMENT] {inst_desc['name']} "
        f"seed_idx={seed_idx} seed_radius={inst.radii[seed_idx]} "
        f"p={inst.p} n={inst.n}",
        flush=True
    )

    for encoding in encodings:
        for solver_name in solvers:
            for run_id in range(1):  # có thể tăng lên 5 nếu muốn như paper
                print(
                    f"[RUN] instance={inst_desc['name']} run_id={run_id + 1} "
                    f"encoding={encoding} solver={solver_name}",
                    flush=True
                )
                status, best_radius, elapsed_time, nvars, nclauses, centers, best_sat_cpu = search_min_radius_parallel(
                    inst,
                    encoding,
                    solver_name,
                    time_limit,
                    radii_workers=radii_workers,
                    seed_idx=seed_idx,
                    mgr=mgr,
                    cancel_dict=cancel_dict,
                )

                print(
                    f"[RUN-RESULT] instance={inst_desc['name']} run_id={run_id + 1} "
                    f"status={status} best_radius={best_radius} "
                    f"elapsed_wall={elapsed_time:.6f}s cpu_time_at_bestR={best_sat_cpu}",
                    flush=True
                )

                results.append(
                    {
                        "instance": inst_desc["name"],
                        "n": inst.n,
                        "p": inst.p,
                        "encoding": encoding,
                        "solver": solver_name,
                        "run_id": run_id + 1,
                        "status": status,
                        "best_radius": best_radius if best_radius is not None else None,
                        "time_sec": elapsed_time,
                        "cpu_time_at_bestR": best_sat_cpu,
                        "nvars": nvars,
                        "nclauses": nclauses,
                        "centers": json.dumps(centers if centers is not None else []),
                    }
                )

    print_instance_summary_for_console(results)

    return results


def print_instance_summary_for_console(all_results_for_inst):
    # Gom theo cấu hình
    cfg_runs = defaultdict(list)  # (enc,sol) -> [rows]
    inst_name = None
    n = None
    p = None

    for r in all_results_for_inst:
        inst_name = r["instance"]
        n = r["n"]
        p = r["p"]
        cfg_runs[(r["encoding"], r["solver"])].append(r)

    # bestR của từng method
    cfg_bestR = {}
    for key, runs in cfg_runs.items():
        brs = [x["best_radius"] for x in runs if x.get("best_radius") is not None]
        if brs:
            cfg_bestR[key] = min(brs)

    if not cfg_bestR:
        print(f"instance={inst_name} n={n} p={p} : no feasible radius found")
        return

    # global best radius trên instance
    gR = min(cfg_bestR.values())

    # chúng ta muốn in theo thứ tự ổn định:
    # ưu tiên solver MapleCM -> MapleChrono -> Sparrow2Riss
    # trong từng solver thì sort encoding alphabetically
    def sort_key(enc_sol):
        enc, sol = enc_sol
        solver_rank = {"maplecm": 0, "maplechrono": 1, "sparrow2riss": 2}
        return solver_rank.get(sol, 99), sol, enc

    method_cols = sorted(cfg_runs.keys(), key=sort_key)

    # In 1 dòng/ method giống format cũ nhưng thay cpu_mean thành mean CPU@R*
    for (enc, sol) in method_cols:
        runs = cfg_runs[(enc, sol)]

        ts = [
            x.get("cpu_time_at_bestR")
            for x in runs
            if x.get("status") == "OK"
            and x.get("best_radius") == gR
            and x.get("cpu_time_at_bestR") is not None
        ]

        if ts:
            mean_cpu = _stats.mean(ts)
            print(
                f"instance={inst_name} n={n} p={p} encoding={enc} solver={sol}: "
                f"radius={gR} cpu_mean={mean_cpu:.3f}s "
            )
        else:
            print(
                f"instance={inst_name} n={n} p={p} encoding={enc} solver={sol}: "
                f"radius={gR} cpu_mean=- "
            )


def write_paper_table(all_results, out_csv_path: str):
    def method_label(solver: str, enc: str) -> str:
        return f"{solver} {enc}"

    # Gom per-run theo cấu hình
    cfg_runs = defaultdict(list)  # (inst,n,p,enc,sol) -> [rows]
    all_instances = set()
    methods_seen = set()

    for r in all_results:
        inst = r["instance"]
        n = r["n"]
        p = r["p"]
        enc = r["encoding"]
        sol = r["solver"]

        all_instances.add((inst, n, p))
        methods_seen.add((enc, sol))

        br = r.get("best_radius")
        if br is None:
            continue
        key = (inst, n, p, enc, sol)
        cfg_runs[key].append(r)

    # Với từng cấu hình: bestR_cfg = min(best_radius)
    cfg_bestR = {}  # (inst,n,p,enc,sol) -> bestR_cfg
    for key, runs in cfg_runs.items():
        brs = [x["best_radius"] for x in runs if x.get("best_radius") is not None]
        if brs:
            cfg_bestR[key] = min(brs)

    # global_bestR[(inst,n,p)] = bán kính nhỏ nhất giữa mọi method
    global_bestR = {}
    for (inst, n, p, enc, sol), R in cfg_bestR.items():
        if (inst, n, p) not in global_bestR or R < global_bestR[(inst, n, p)]:
            global_bestR[(inst, n, p)] = R

    # Thời gian theo method tại đúng global_bestR
    # key: (inst,n,p,enc,sol) -> list times khi run đạt global_bestR
    times_at_global = defaultdict(list)
    for (inst, n, p, enc, sol), runs in cfg_runs.items():
        gR = global_bestR.get((inst, n, p))
        if gR is None:
            continue
        for x in runs:
            # ta chỉ collect time nếu run đạt gR của instance
            if (
                x.get("status") == "OK"
                and x.get("best_radius") == gR
                and x.get("cpu_time_at_bestR") is not None
            ):
                times_at_global[(inst, n, p, enc, sol)].append(
                    x["cpu_time_at_bestR"]
                )

    def sort_key(enc_sol):
        enc, sol = enc_sol
        solver_rank = {"maplecm": 0, "maplechrono": 1, "sparrow2riss": 2}
        return solver_rank.get(sol, 99), sol, enc

    method_cols = sorted(methods_seen, key=sort_key)

    # Header CSV
    header = ["Instance", "n", "p", "radius"] + [
        method_label(solver=sol, enc=enc) for (enc, sol) in method_cols
    ]

    rows = []
    # Lưu lại để tính footer Num./Avg.
    solved_by = {
        method_label(solver=sol, enc=enc): [] for (enc, sol) in method_cols
    }

    # Với mỗi instance (inst,n,p) -> đúng 1 dòng
    for (inst, n, p) in sorted(
        all_instances, key=lambda x: (str(x[0]), x[1], x[2])
    ):
        gR = global_bestR.get((inst, n, p))
        row = {
            "Instance": inst,
            "n": n,
            "p": p,
            "radius": gR if gR is not None else "-",
        }

        for (enc, sol) in method_cols:
            col = method_label(solver=sol, enc=enc)
            ts = times_at_global.get((inst, n, p, enc, sol), [])
            if ts:
                mean_t = _stats.mean(ts)
                row[col] = f"{mean_t:.3f}"
                solved_by[col].append(mean_t)
            else:
                row[col] = "-"

        rows.append(row)

    # Footer Num. và Avg.
    footer_num = {"Instance": "Num.", "n": "", "p": "", "radius": ""}
    footer_avg = {"Instance": "Avg.", "n": "", "p": "", "radius": ""}

    for (enc, sol) in method_cols:
        col = method_label(solver=sol, enc=enc)
        lst = solved_by[col]
        if lst:
            footer_num[col] = str(len(lst))
            footer_avg[col] = f"{_stats.mean(lst):.3f}"
        else:
            footer_num[col] = "0"
            footer_avg[col] = "-"

    rows.append(footer_num)
    rows.append(footer_avg)

    # Ghi CSV
    with open(out_csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    MGR = mp.Manager()
    CANCEL = MGR.dict()

    args = parse_args()
    instance_data = load_instance_data(args.instances)

    # Chuẩn hoá đường dẫn tương đối theo vị trí instances.json
    base_dir = os.path.dirname(os.path.abspath(args.instances))
    for d in instance_data:
        if "orlib" in d and not os.path.isabs(d["orlib"]):
            d["orlib"] = os.path.join(base_dir, d["orlib"])

    all_results = []
    for inst_desc in instance_data:
        CANCEL.clear()

        print(f"[BEGIN] {inst_desc['name']}", flush=True)
        try:
            res = run_experiment(
                inst_desc,
                args.encodings,
                args.solvers,
                args.time_limit,
                args.radii_workers,
                args.strategy,
                mgr=MGR,
                cancel_dict=CANCEL,
            )
            all_results.extend(res)
            print(f"[END] {inst_desc['name']}", flush=True)
        except Exception as e:
            import traceback

            traceback.print_exc()
            print(
                f"[WARN] Bỏ qua {inst_desc.get('name')} do lỗi: {e}",
                flush=True
            )

    if all_results:
        write_paper_table(all_results, args.out)

