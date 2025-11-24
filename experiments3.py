import argparse
import json
import time
import heapq
import os, sys
import csv
import multiprocessing as mp
import statistics as _stats
import threading
from pysat.solvers import Solver
from encoder import PCenterSAT
from collections import defaultdict
from typing import List, Tuple

_CANCEL_SHARED = None 
INF = 10**12

def _pool_initializer(cancel_proxy):
    global _CANCEL_SHARED
    _CANCEL_SHARED = cancel_proxy

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
        seed_idx = next((i for i, r in enumerate(inst.radii) if r <= sr + 1e-12), None)

    elif "orlib" in inst_desc:
        t0 = time.time()
        n_graph, p, edges = load_orlib_edge_list(inst_desc["orlib"])
        print(f"[LOAD] {inst_desc['name']}: n={n_graph}, m={len(edges)}, p={p}", flush=True)

        t1 = time.time()
        D = compute_apsp_dijkstra(n_graph, edges)
        print(f"[APSP] {inst_desc['name']}: done in {time.time()-t1:.2f}s (total load {time.time()-t0:.2f}s)", flush=True)
        inst = PCenterSAT.from_distance_matrix(D, p)

        radii_set = set()
        for i in range(n_graph):
            for j in range(i+1, n_graph):
                dij = D[i][j]
                if dij != INF and dij > 0:
                    radii_set.add(dij)
        
        inst.radii = sorted(radii_set, reverse=True)

        sr = int(inst_desc["seed_radius"])
        seed_idx = next((i for i, r in enumerate(inst.radii) if r <= sr), None)

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
            u, v, w = int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2]) 

            edges.append((u, v, w))

    return n, p, edges


def compute_apsp_dijkstra(n: int, edges: List[Tuple[int, int, float]]):
    c = [[INF] * n for _ in range(n)]
    for i in range(n):
        c[i][i] = 0

    for u, v, w in edges:
        ww = int(w)
        if u == v:
            c[u][v] = 0
        else:
            # dòng sau cùng sẽ ghi đè (chuẩn OR-Lib)
            c[u][v] = ww
            c[v][u] = ww

    # 2) Chuyển sang adjacency (bỏ INF)
    adj = [[] for _ in range(n)]
    for u in range(n):
        row = c[u]
        for v in range(n):
            w = row[v]
            if u != v and w < INF:
                adj[u].append((v, w))

    # 3) APSP bằng Dijkstra-from-each-source
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

def _apply_async_wrapper(args):
    # tiện để Pool.apply_async nhận 1 tuple args
    return _solve_radius_worker_proc(*args)

def _iter_async_results(async_map):
    """Yield kết quả theo dạng 'as_completed' cho multiprocessing.Pool.apply_async."""
    pending = dict(async_map)  # {AsyncResult: k}
    while pending:
        done = []
        for h, k in list(pending.items()):
            if h.ready():
                try:
                    res = h.get()
                except Exception:
                    res = (k, None, "timeout", 0.0, None, None, None)
                yield res
                done.append(h)
        for h in done:
            pending.pop(h, None)
        if pending:
            time.sleep(0.01)

def search_min_radius_parallel(inst, encoding, solver_name, time_limit, *, radii_workers=8, seed_idx=None, mgr=None, cancel_dict=None):
    radii = inst.radii
    nR = len(radii)
    start_wall = time.time()
    if nR == 0:
        return "infeasible", None, 0.0, None, None, None, None

    dist_raw = inst.dist
    p_raw = inst.p

    decided = {}
    best_nvars = None
    best_nclauses = None
    sat_solutions = {}
    sat_cpu_time = {}

    # seed mặc định = 0 (bán kính lớn nhất)
    i = seed_idx if (seed_idx is not None and 0 <= seed_idx < nR) else 0
    covered = set()
    best_sat_idx = None

    # helper: tạo Event cho một idx nếu chưa có
    def ensure_event(idx):
        if cancel_dict is not None:
            if idx not in cancel_dict:
                cancel_dict[idx] = mgr.Event()
            return cancel_dict[idx]
        return None

    # submit batch lên pool
    def launch_batch(pool, idxs):
        asyncs = {}
        for k in idxs:
            ensure_event(k)
            R = radii[k]
            a = (k, dist_raw, p_raw, encoding, solver_name, R, time_limit)
            h = pool.apply_async(_apply_async_wrapper, (a,))
            asyncs[h] = k
        return asyncs

    with mp.get_context("spawn").Pool(
        processes=radii_workers,
        initializer=_pool_initializer,
        initargs=(cancel_dict,)
    ) as pool:
        while i < nR:
            j = min(nR - 1, i + radii_workers - 1)
            need = [k for k in range(i, j + 1) if k not in covered]
            if not need:
                i = j + 1
                continue

            asyncs = launch_batch(pool, need)
            best_sat_in_batch = None

            for idx, R, status, t_sec, nvars, nclauses, centers in _iter_async_results(asyncs):
                if R is None:
                    R = radii[idx]
                decided[idx] = status

                if nvars is not None and nclauses is not None:
                    best_nvars = nvars if best_nvars is None else min(best_nvars, nvars)
                    best_nclauses = nclauses if best_nclauses is None else min(best_nclauses, nclauses)

                if status == "sat":
                    sat_solutions[idx] = centers
                    sat_cpu_time[idx] = t_sec
                    if (best_sat_in_batch is None) or (idx > best_sat_in_batch):
                        best_sat_in_batch = idx
                    if cancel_dict is not None:
                        for kk in need:
                            if kk < best_sat_in_batch:
                                try:
                                    cancel_dict[kk].set()
                                except Exception:
                                    pass

            covered.update(need)

            if best_sat_in_batch is not None:
                k = best_sat_in_batch
                while k + 1 < nR:
                    nxt = k + 1
                    if nxt in decided:
                        if decided[nxt] == "unsat":
                            best_sat_idx = k
                            break
                        elif decided[nxt] == "sat":
                            k = nxt
                            continue

                    ensure_event(nxt)
                    asyncs2 = launch_batch(pool, [nxt])
                    h2 = next(iter(asyncs2.keys()))
                    try:
                        idx2, R2, status2, t2, nv2, nc2, centers2 = h2.get()
                    except Exception:
                        status2, t2, nv2, nc2, centers2 = "timeout", 0.0, None, None, None

                    decided[nxt] = status2
                    if nv2 is not None and nc2 is not None:
                        best_nvars = nv2 if best_nvars is None else min(best_nvars, nv2)
                        best_nclauses = nc2 if best_nclauses is None else min(best_nclauses, nc2)

                    if status2 == "unsat":
                        best_sat_idx = k
                        break
                    elif status2 == "sat":
                        sat_solutions[nxt] = centers2
                        sat_cpu_time[nxt] = t2
                        if cancel_dict is not None:
                            for kk in range(0, nxt):
                                if kk in cancel_dict:
                                    try: cancel_dict[kk].set()
                                    except Exception: pass
                        k = nxt
                        continue
                    else:
                        break

                if best_sat_idx is not None:
                    break

            i = j + 1

    elapsed = time.time() - start_wall

    if best_sat_idx is None:
        if sat_solutions:
            best_sat_idx = max(sat_solutions.keys())
        else:
            return "infeasible", None, elapsed, best_nvars, best_nclauses, None, None

    best_centers = sat_solutions.get(best_sat_idx, None)
    best_sat_cpu = sat_cpu_time.get(best_sat_idx, None)
    return "OK", radii[best_sat_idx], elapsed, best_nvars, best_nclauses, best_centers, best_sat_cpu



def _solve_radius_worker_proc(idx, dist, p, encoding, solver_name, radius, time_limit):
    """
    Chạy trong process riêng. Trả về:
    (idx, radius, status, solver_time, nvars, nclauses, centers)
    """
    try:
        inst = PCenterSAT.from_distance_matrix(dist, p)
        cnf, varmap = inst._encode_cnf(radius, encoding)
        if cnf is None:
            return (idx, radius, "unsat", 0.0, None, None, None)

        with Solver(name=solver_name, bootstrap_with=cnf.clauses) as solver:
            cancel_ev = None
            try:
                if _CANCEL_SHARED is not None:
                    cancel_ev = _CANCEL_SHARED.get(idx)
            except Exception:
                cancel_ev = None

            watcher = None
            if cancel_ev is not None and hasattr(solver, "interrupt"):
                def _watch():
                    cancel_ev.wait()
                    try:
                        solver.interrupt()
                    except Exception:
                        pass
                watcher = threading.Thread(target=_watch, daemon=True)
                watcher.start()

            timer = None
            try: 
                if time_limit and time_limit > 0 and hasattr(solver, "interrupt"):
                    timer = threading.Timer(time_limit, solver.interrupt)
                    timer.start()
                    t0 = time.process_time_ns()
                    sat = solver.solve_limited(expect_interrupt=True)
                    t1 = time.process_time_ns()
                else:
                    t0 = time.process_time_ns()
                    sat = solver.solve()
                    t1 = time.process_time_ns()

                solver_time = (t1 - t0) / 1e9

                nvars = cnf.nv
                nclauses = len(cnf.clauses)

                if sat is None:
                    return (idx, radius, "timeout", None, nvars, nclauses, None)
                if not sat:
                    return (idx, radius, "unsat", solver_time, nvars, nclauses, None)

                model = solver.get_model() or []
                if not model:
                    return (idx, radius, "timeout", solver_time, nvars, nclauses, None)

                model_set = set(model)
                y_vars = varmap.get("y", [])
                Nc = set(varmap.get("Nc", []))
                candidates = set(varmap.get("candidates", []))
                chosen = {j for j, v in enumerate(y_vars) if (v in model_set) and (j in candidates)}
                centers = sorted(chosen | Nc)
                return (idx, radius, "sat", solver_time, nvars, nclauses, centers)
            finally:
                if timer:
                    timer.cancel()
    except Exception:
        return (idx, radius, "timeout", 0.0, None, None, None)




# ---------------- CLI & Runner ----------------

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instances", type=str, required=True, help="Path to instances.json")
    ap.add_argument("--encodings", nargs="+", default=["seq"], help="Encodings to test: seq par")
    ap.add_argument("--solvers", nargs="+", default=["maplecm", "maplechrono", "riss"], help="SAT solvers to use")
    ap.add_argument("--time-limit", type=int, default=14400, help="Time limit per radius solve (seconds)")
    ap.add_argument("--radii-workers", type=int, default=8, help="Parallel radii workers (paper uses 8)")
    ap.add_argument("--out", type=str, default="results.csv", help="CSV output path")
    ap.add_argument("--strategy", type=str, choices=["parallel"], default="parallel",
    help="Radius search strategy (paper compares parallel only)")
    return ap.parse_args()

def run_experiment(inst_desc, encodings, solvers, time_limit, radii_workers, strategy, *, mgr, cancel_dict):
    results = []
    inst, seed_idx = load_instance(inst_desc)

    for encoding in encodings:
        for solver_name in solvers:
            for run_id in range(5):
                status, best_radius, elapsed_time, nvars, nclauses, centers, best_sat_cpu = search_min_radius_parallel(
                    inst, encoding, solver_name, time_limit,
                    radii_workers=radii_workers,
                    seed_idx=seed_idx,
                    mgr=mgr,
                    cancel_dict=cancel_dict
                )

                results.append({
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
                    "centers": json.dumps(centers if centers is not None else [])
                })

    print_instance_summary_for_console(results)

    return results

def print_instance_summary_for_console(all_results_for_inst):
    # Gom theo cấu hình
    cfg_runs = defaultdict(list)  # (enc,sol) -> [rows]
    inst_name = None; n = None; p = None
    for r in all_results_for_inst:
        inst_name = r["instance"]; n = r["n"]; p = r["p"]
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

    # Tính mean CPU time tại đúng gR cho từng method
    # method_cols theo thứ tự đẹp
    canonical = [("seq","maplecm"), ("seq","maplechrono"), ("seq","riss"),
                 ("par","maplecm"), ("par","maplechrono"), ("par","riss")]
    methods_seen = list(cfg_runs.keys())
    method_cols = [m for m in canonical if m in methods_seen] + [m for m in methods_seen if m not in canonical]

    def mlabel(enc, sol):
        sm = {"maplecm":"Maple","maplechrono":"MapleChrono","riss":"Riss"}.get(sol, sol)
        em = {"seq":"SEQ","par":"PAR"}.get(enc, enc).upper()
        return f"{sm} {em}"

    # In 1 dòng/ method giống format cũ nhưng thay cpu_mean thành mean CPU@R*
    for (enc, sol) in method_cols:
        runs = cfg_runs[(enc, sol)]
        ts = [x.get("cpu_time_at_bestR") for x in runs if x.get("status") == "OK" and x.get("best_radius") == gR and x.get("cpu_time_at_bestR") is not None]
        if ts:
            mean_cpu = _stats.mean(ts)
            print(f"instance={inst_name} n={n} p={p} encoding={enc} solver={sol}: "
                  f"radius={gR} cpu_mean={mean_cpu:.3f}s ")
        else:
            print(f"instance={inst_name} n={n} p={p} encoding={enc} solver={sol}: "
                  f"radius={gR} cpu_mean=- ")

def write_paper_table(all_results, out_csv_path: str):

    def method_label(solver: str, enc: str) -> str:
        solver_map = {
            "maplecm": "Maple",
            "maplechrono": "MapleChrono",
            "riss": "Riss",
        }
        enc_map = {"seq": "SEQ", "par": "PAR"}
        return f"{solver_map.get(solver, solver)} {enc_map.get(enc, enc).upper()}"

    # Gom per-run theo cấu hình
    cfg_runs = defaultdict(list)  # (inst,n,p,enc,sol) -> [rows]
    for r in all_results:
        br = r.get("best_radius")
        if br is None:
            continue
        key = (r["instance"], r["n"], r["p"], r["encoding"], r["solver"])
        cfg_runs[key].append(r)

    # Với từng cấu hình: bestR_cfg = min(best_radius)
    cfg_bestR = {}  # (inst,n,p,enc,sol) -> bestR_cfg
    for key, runs in cfg_runs.items():
        cfg_bestR[key] = min(x["best_radius"] for x in runs if x.get("best_radius") is not None)

    inst_methods = defaultdict(set)  # (inst,n,p) -> {(enc,sol)}
    global_bestR = {}
    for (inst, n, p, enc, sol), R in cfg_bestR.items():
        inst_methods[(inst, n, p)].add((enc, sol))
        if (inst, n, p) not in global_bestR or R < global_bestR[(inst, n, p)]:
            global_bestR[(inst, n, p)] = R

    # Thời gian theo method tại đúng global_bestR
    # key: (inst,n,p,enc,sol) -> list times khi run đạt global_bestR
    times_at_global = defaultdict(list)
    for (inst, n, p, enc, sol), runs in cfg_runs.items():
        gR = global_bestR[(inst, n, p)]
        for x in runs:
            if x.get("status") == "OK" and x.get("best_radius") == gR:
                t = x.get("cpu_time_at_bestR")
                if t is not None:
                    times_at_global[(inst, n, p, enc, sol)].append(t)

    # Xác định danh sách cột method theo thứ tự “đẹp”
    # Nếu bạn chỉ chạy maplecm+seq thì bảng sẽ chỉ có 1 cột "Maple SEQ".
    methods_seen = set((enc, sol) for (_, _, _, enc, sol) in cfg_runs.keys())
    canonical = [("seq", "maplecm"), ("seq", "maplechrono"), ("seq", "riss"),
                 ("par", "maplecm"), ("par", "maplechrono"), ("par", "riss")]
    method_cols = []
    for m in canonical:
        if m in methods_seen:
            method_cols.append(m)
    # Nếu còn method lạ khác thứ tự, thêm vào cuối
    for m in sorted(methods_seen):
        if m not in method_cols:
            method_cols.append(m)

    header = ["Instance", "n", "p", "radius"] + [method_label(solver=m[1], enc=m[0]) for m in method_cols]

    # Tạo các dòng dữ liệu: mỗi instance đúng 1 dòng
    rows = []
    solved_by = {method_label(solver=m[1], enc=m[0]): [] for m in method_cols}
    for (inst, n, p) in sorted(inst_methods.keys(), key=lambda x: (str(x[0]), x[1], x[2])):
        gR = global_bestR[(inst, n, p)]
        row = {"Instance": inst, "n": n, "p": p, "radius": gR}
        for enc, sol in method_cols:
            col = method_label(solver=sol, enc=enc)
            ts = times_at_global.get((inst, n, p, enc, sol))
            if ts:
                mean_t = _stats.mean(ts)
                row[col] = f"{mean_t:.3f}"
                solved_by[col].append(mean_t)
            else:
                row[col] = "-"
        rows.append(row)

    # Hai dòng footer
    footer_num = {"Instance": "Num.", "n": "", "p": "", "radius": ""}
    footer_avg = {"Instance": "Avg.", "n": "", "p": "", "radius": ""}
    for enc, sol in method_cols:
        col = method_label(solver=sol, enc=enc)
        lst = solved_by[col]
        footer_num[col] = str(len(lst)) if lst else "0"
        footer_avg[col] = f"{_stats.mean(lst):.3f}" if lst else "-"
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
                cancel_dict=CANCEL
            )
            all_results.extend(res)
            print(f"[END] {inst_desc['name']}", flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"[WARN] Bỏ qua {inst_desc.get('name')} do lỗi: {e}", flush=True)

    if all_results:
        write_paper_table(all_results, args.out)
