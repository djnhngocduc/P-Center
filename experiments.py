import argparse
import json
import time
import os
import multiprocessing as mp
from collections import defaultdict
import statistics as _stats
import csv
import traceback

from strategies.hybrid_race import search_min_radius_hybrid_race

from utils.instances import (
    load_instance_data,
    load_instance,
)

from strategies.incremental import search_min_radius_incremental
from strategies.maxsat import search_min_radius_maxsat

from strategies.search import (
    search_min_radius_parallel,
    search_min_radius_binary, 
    search_min_radius_kary,
)

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
        default=["pysat_totalizer", "pysat_mtotalizer", "pysat_kmtotalizer", "pypb_bdd", "nsc", "pb_bdd"],
        help="Encodings to test: pysat_totalizer pysat_mtotalizer pysat_kmtotalizer pypb_bdd nsc pb_bdd",
    )
    ap.add_argument(
        "--solvers",
        nargs="+",
        default=["maplecm", "maplechrono", "sparrow2riss", "glucose4", "kissat", "cplex_mip", "cplex_cp", "gurobi_mip"],
        help="Solvers to use (internal SAT / external SAT / cplex_mip / cplex_cp / gurobi_mip)",
    )
    ap.add_argument(
        "--time-limit",
        type=int,
        default=14400,
        help="Total time limit per instance run (seconds)",
    )
    ap.add_argument(
        "--search-mode",
        type=str,
        default="parallel",
        choices=["parallel", "binary", "kary", "incremental", "hybrid_race", "maxsat"],
        help="Search strategy: parallel, binary, kary, incremental, hybrid_race, or maxsat",
    )
    ap.add_argument(
        "--mip-backend",
        type=str,
        default="cplex_mip",
        choices=["cplex_mip", "gurobi_mip"],
        help="MIP backend used internally by hybrid_race: cplex_mip or gurobi_mip",
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
        default=os.path.join("results", "results.csv"),
        help="CSV output path",
    )
    return ap.parse_args()

def run_experiment(
    inst_desc,
    encodings,
    solvers,
    time_limit,
    search_mode,
    radii_workers,
    *,
    mip_backend="cplex_mip",
    mgr,
    cancel_dict,
):
    if "cplex_mip" in solvers and search_mode not in ("binary", "hybrid_race"):
        raise ValueError("CPLEX backend currently supports only binary search or hybrid_race.")
    if "cplex_cp" in solvers and search_mode not in ("parallel", "binary", "kary"):
        raise ValueError("CPO backend currently supports only parallel, binary, or kary search.")
    if "gurobi_mip" in solvers and search_mode not in ("parallel", "binary", "kary"):
        raise ValueError(
            "Gurobi backend in --solvers only supports parallel, binary, or kary search. "
            "For hybrid modes, use --mip-backend gurobi_mip instead."
        )
    
    results = []
    load_t0 = time.perf_counter()
    inst, seed_idx = load_instance(inst_desc, use_seed_radius=False)
    load_time = time.perf_counter() - load_t0
    print(
        f"[RUN-EXPERIMENT] {inst_desc['name']} "
        f"seed_idx={seed_idx} seed_radius={inst.radii[seed_idx]} "
        f"p={inst.p} n={inst.n} load_time={load_time:.6f}s",
        flush=True
    )

    if search_mode == "hybrid_race":
        if len(solvers) != 1:
            raise ValueError(
                "hybrid_race expects exactly one SAT solver in --solvers. "
                "The MIP backend is selected separately by --mip-backend."
            )
        if solvers[0] in ("cplex_mip", "gurobi_mip", "cplex_cp"):
            raise ValueError(
                "hybrid_race needs a SAT solver in --solvers, not a MIP backend. "
                "Use --mip-backend to choose cplex_mip or gurobi_mip."
            )

    if search_mode == "maxsat":
        if len(encodings) != 1 or encodings[0] not in ("maxsat_setcover", "maxsat_cover"):
            raise ValueError(
                "maxsat mode expects exactly one encoding: --encodings maxsat_setcover"
            )

        for s in solvers:
            if s not in ("rc2", "openwbo"):
                raise ValueError(
                    "maxsat mode supports only --solvers rc2, openwbo"
                )

    for solver_name in solvers:
        solver_encodings = encodings if solver_name not in ("cplex_mip", "cplex_cp", "gurobi_mip") else ["setcover"]
        for encoding in solver_encodings:
            for run_id in range(1):
                print(
                    f"[RUN] instance={inst_desc['name']} run_id={run_id + 1} "
                    f"encoding={encoding} solver={solver_name} search={search_mode}",
                    flush=True
                )
                if search_mode == "parallel":
                    search_fn = search_min_radius_parallel
                elif search_mode == "binary":
                    search_fn = search_min_radius_binary
                elif search_mode == "kary":
                    search_fn = search_min_radius_kary
                elif search_mode == "incremental":
                    search_fn = search_min_radius_incremental
                elif search_mode == "hybrid_race":
                    search_fn = search_min_radius_hybrid_race
                elif search_mode == "maxsat":
                    search_fn = search_min_radius_maxsat
                else:
                    raise ValueError(f"Unknown search_mode: {search_mode}")

                cancel_dict.clear()
                
                extra_kwargs = dict(
                    radii_workers=radii_workers,
                    seed_idx=seed_idx,
                    mgr=mgr,
                    cancel_dict=cancel_dict,
                )

                if search_mode == "hybrid_race":
                    extra_kwargs["mip_backend"] = mip_backend

                status, best_radius, nvars, nclauses, centers, search_elapsed = search_fn(
                    inst,
                    encoding,
                    solver_name,
                    time_limit,
                    **extra_kwargs,
                )

                total_time = load_time + search_elapsed

                print(
                    f"[RUN-RESULT] instance={inst_desc['name']} run_id={run_id + 1} "
                    f"status={status} best_radius={best_radius} "
                    f"load_time={load_time:.6f}s " 
                    f"search_time={search_elapsed:.6f}s total_time={total_time:.6f}s",
                    flush=True
                )

                results.append(
                    {
                        "instance": inst_desc["name"],
                        "n": inst.n,
                        "p": inst.p,
                        "encoding": encoding,
                        "solver": solver_name,
                        "search_mode": search_mode,
                        "run_id": run_id + 1,
                        "status": status,
                        "best_radius": best_radius if best_radius is not None else None,
                        "load_time": load_time,
                        "search_time": search_elapsed,
                        "total_time": total_time,
                        "nvars": nvars,
                        "nclauses": nclauses,
                        "centers": json.dumps(centers if centers is not None else []),
                    }
                )

    print_instance_summary_for_console(results)

    return results

def sort_key(enc_sol_mode):
    enc, sol, mode = enc_sol_mode
    solver_rank = {"maplecm": 0, "maplechrono": 1, "sparrow2riss": 2, "glucose4": 3, "kissat": 5, "rc2": 6, "openwbo": 7, "cplex_mip": 8, "cplex_cp": 9, "gurobi_mip": 10}
    mode_rank = {"parallel": 0, "binary": 1, "kary": 2, "incremental": 3, "maxsat": 4, "hybrid_race": 5}
    return solver_rank.get(sol, 99), sol, mode_rank.get(mode, 99), mode, enc

def print_instance_summary_for_console(all_results_for_inst):
    cfg_runs = defaultdict(list)
    inst_name = None
    n = None
    p = None

    for r in all_results_for_inst:
        inst_name = r["instance"]
        n = r["n"]
        p = r["p"]
        mode = r.get("search_mode", "parallel")
        cfg_runs[(r["encoding"], r["solver"], mode)].append(r)

    cfg_bestR = {}
    for key, runs in cfg_runs.items():
        brs = [x["best_radius"] for x in runs if x.get("best_radius") is not None]
        if brs:
            cfg_bestR[key] = min(brs)

    if not cfg_bestR:
        print(f"instance={inst_name} n={n} p={p} : no feasible radius found")
        return

    method_cols = sorted(cfg_runs.keys(), key=sort_key)

    for (enc, sol, mode) in method_cols:
        runs = cfg_runs[(enc, sol, mode)]
        gR = cfg_bestR.get((enc, sol, mode))

        loads = [
            x.get("load_time")
            for x in runs
            if x.get("status") in ("OK", "uncertified")
            and x.get("best_radius") == gR
            and x.get("load_time") is not None
        ]

        searches = [
            x.get("search_time")
            for x in runs
            if x.get("status") in ("OK", "uncertified")
            and x.get("best_radius") == gR
            and x.get("search_time") is not None
        ]

        totals = [
            x.get("total_time")
            for x in runs
            if x.get("status") in ("OK", "uncertified")
            and x.get("best_radius") == gR
            and x.get("total_time") is not None
        ]

        load_str = f"{_stats.mean(loads):.3f}s" if loads else "-"
        search_str = f"{_stats.mean(searches):.3f}s" if searches else "-"
        total_str = f"{_stats.mean(totals):.3f}s" if totals else "-"

        print(
            f"instance={inst_name} n={n} p={p} encoding={enc} solver={sol} mode={mode}: "
            f"radius={gR} load_mean={load_str} "
            f"search_mean={search_str} total_mean={total_str}"
        )


def write_paper_table(all_results, out_csv_path: str):
    def method_label(solver: str, enc: str, mode: str) -> str:
        return f"{solver} {enc} {mode}"

    cfg_runs = defaultdict(list)
    all_instances = set()
    methods_seen = set()
    sizes_at_global = {}

    for r in all_results:
        inst = r["instance"]
        n = r["n"]
        p = r["p"]
        enc = r["encoding"]
        sol = r["solver"]
        mode = r.get("search_mode", "parallel")

        all_instances.add((inst, n, p))
        methods_seen.add((enc, sol, mode))

        br = r.get("best_radius")
        if br is None:
            continue

        key = (inst, n, p, enc, sol, mode)
        cfg_runs[key].append(r)

    cfg_bestR = {}
    for key, runs in cfg_runs.items():
        brs = [x["best_radius"] for x in runs if x.get("best_radius") is not None]
        if brs:
            cfg_bestR[key] = min(brs)

    global_bestR = {}
    for (inst, n, p, enc, sol, mode), R in cfg_bestR.items():
        k = (inst, n, p, enc, sol, mode)
        if k not in global_bestR or R < global_bestR[k]:
            global_bestR[k] = R

    load_times_at_global = defaultdict(list)
    search_times_at_global = defaultdict(list)
    total_times_at_global = defaultdict(list)

    for (inst, n, p, enc, sol, mode), runs in cfg_runs.items():
        gR = global_bestR.get((inst, n, p, enc, sol, mode))
        if gR is None:
            continue

        for x in runs:
            if x.get("status") not in ("OK", "uncertified"):
                continue
            if x.get("best_radius") != gR:
                continue

            key_size = (inst, n, p, enc, sol, mode)
            if key_size not in sizes_at_global:
                sizes_at_global[key_size] = (x.get("nvars"), x.get("nclauses"))

            load_v = x.get("load_time")
            if load_v is not None:
                load_times_at_global[(inst, n, p, enc, sol, mode)].append(load_v)

            search_v = x.get("search_time")
            if search_v is not None:
                search_times_at_global[(inst, n, p, enc, sol, mode)].append(search_v)

            total_v = x.get("total_time")
            if total_v is not None:
                total_times_at_global[(inst, n, p, enc, sol, mode)].append(total_v)

    method_cols = sorted(methods_seen, key=sort_key)

    header = ["Instance", "n", "p"]

    for (enc, sol, mode) in method_cols:
        m = method_label(solver=sol, enc=enc, mode=mode)
        header.append(f"{m} radius")
        header.append(f"{m} #vars")
        header.append(f"{m} #clauses")
        header.append(f"{m} load")
        header.append(f"{m} search")
        header.append(f"{m} total")

    rows = []

    solved_load = {method_label(sol, enc, mode): [] for (enc, sol, mode) in method_cols}
    solved_search = {method_label(sol, enc, mode): [] for (enc, sol, mode) in method_cols}
    solved_total = {method_label(sol, enc, mode): [] for (enc, sol, mode) in method_cols}

    for (inst, n, p) in sorted(all_instances, key=lambda x: (str(x[0]), x[1], x[2])):
        row = {
            "Instance": inst,
            "n": n,
            "p": p,
        }

        for (enc, sol, mode) in method_cols:
            m = method_label(solver=sol, enc=enc, mode=mode)

            radius_col = f"{m} radius"
            vars_col = f"{m} #vars"
            clauses_col = f"{m} #clauses"
            load_col = f"{m} load"
            search_col = f"{m} search"
            total_col = f"{m} total"

            gR = global_bestR.get((inst, n, p, enc, sol, mode))
            row[radius_col] = gR if gR is not None else "-"

            nv, nc = sizes_at_global.get((inst, n, p, enc, sol, mode), ("-", "-"))
            row[vars_col] = nv
            row[clauses_col] = nc

            ts_load = load_times_at_global.get((inst, n, p, enc, sol, mode), [])
            if ts_load:
                mean_load = _stats.mean(ts_load)
                row[load_col] = f"{mean_load:.3f}"
                solved_load[m].append(mean_load)
            else:
                row[load_col] = "-"

            ts_search = search_times_at_global.get((inst, n, p, enc, sol, mode), [])
            if ts_search:
                mean_search = _stats.mean(ts_search)
                row[search_col] = f"{mean_search:.3f}"
                solved_search[m].append(mean_search)
            else:
                row[search_col] = "-"

            ts_total = total_times_at_global.get((inst, n, p, enc, sol, mode), [])
            if ts_total:
                mean_total = _stats.mean(ts_total)
                row[total_col] = f"{mean_total:.3f}"
                solved_total[m].append(mean_total)
            else:
                row[total_col] = "-"

        rows.append(row)

    footer_num = {"Instance": "Num.", "n": "", "p": ""}
    footer_avg = {"Instance": "Avg.", "n": "", "p": ""}

    for (enc, sol, mode) in method_cols:
        m = method_label(solver=sol, enc=enc, mode=mode)

        radius_col = f"{m} radius"
        vars_col = f"{m} #vars"
        clauses_col = f"{m} #clauses"
        load_col = f"{m} load"
        search_col = f"{m} search"
        total_col = f"{m} total"

        footer_num[radius_col] = ""
        footer_num[vars_col] = ""
        footer_num[clauses_col] = ""
        footer_avg[radius_col] = ""
        footer_avg[vars_col] = ""
        footer_avg[clauses_col] = ""

        lst_load = solved_load.get(m, [])
        lst_search = solved_search.get(m, [])
        lst_total = solved_total.get(m, [])

        footer_num[load_col] = str(len(lst_load)) if lst_load else "0"
        footer_avg[load_col] = f"{_stats.mean(lst_load):.3f}" if lst_load else "-"

        footer_num[search_col] = str(len(lst_search)) if lst_search else "0"
        footer_avg[search_col] = f"{_stats.mean(lst_search):.3f}" if lst_search else "-"

        footer_num[total_col] = str(len(lst_total)) if lst_total else "0"
        footer_avg[total_col] = f"{_stats.mean(lst_total):.3f}" if lst_total else "-"

    rows.append(footer_num)
    rows.append(footer_avg)

    with open(out_csv_path, "w", newline="", encoding="utf-8") as f:
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
                args.search_mode,
                args.radii_workers,
                mip_backend=args.mip_backend,
                mgr=MGR,
                cancel_dict=CANCEL,
            )
            all_results.extend(res)
            print(f"[END] {inst_desc['name']}", flush=True)
        except Exception as e:
            traceback.print_exc()
            print(
                f"[WARN] Bỏ qua {inst_desc.get('name')} do lỗi: {e}",
                flush=True
            )

    if all_results:
        write_paper_table(all_results, args.out)

