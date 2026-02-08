import argparse
import json
import csv
from collections import defaultdict

EPS = 1e-12


# --------------------------------------------------
# Compute coverage symmetry
# --------------------------------------------------

def compute_coverage_symmetry(inst, radius):
    Nc, Nd, enabled_centers, demands, _ = inst.compute_reduction(radius)

    candidates = sorted(enabled_centers - Nc - Nd)

    signature = defaultdict(list)

    for c in candidates:
        covered = []
        row = inst.dist[c]
        for u in demands:
            if row[u] <= radius + EPS:
                covered.append(u)

        signature[tuple(covered)].append(c)

    sym_classes = {
        k: v for k, v in signature.items()
        if len(v) >= 2
    }

    return sym_classes, Nc, Nd, candidates


# --------------------------------------------------
# Load instance EXACTLY like experiments.py
# --------------------------------------------------

def load_instance_and_seed_radius(inst_desc):
    from experiments import load_instance
    inst, seed_idx = load_instance(inst_desc)
    radius = inst.radii[seed_idx]
    return inst, radius


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", required=True)
    parser.add_argument("--out-prefix", default="symmetry")
    args = parser.parse_args()

    with open(args.instances, "r", encoding="utf-8") as f:
        data = json.load(f)

    summary_rows = []
    detail_rows = []

    for inst_desc in data:
        name = inst_desc["name"]
        p = inst_desc["p"]

        print(f"[PROCESS] {name}, p={p}")

        inst, radius = load_instance_and_seed_radius(inst_desc)
        print(f"  mapped radius = {radius}")

        sym_classes, Nc, Nd, candidates = compute_coverage_symmetry(inst, radius)

        total_sym_nodes = sum(len(v) for v in sym_classes.values())
        max_class = max((len(v) for v in sym_classes.values()), default=0)
        num_classes = len(sym_classes)
        num_candidates = len(candidates)

        symmetry_ratio = (
            total_sym_nodes / num_candidates
            if num_candidates > 0 else 0
        )

        # -------- Summary row
        summary_rows.append({
            "instance": name,
            "p": p,
            "radius_used": radius,
            "num_candidates": num_candidates,
            "num_symmetry_classes": num_classes,
            "total_symmetric_nodes": total_sym_nodes,
            "max_class_size": max_class,
            "symmetry_ratio": symmetry_ratio
        })

        # -------- Detail rows
        for idx, centers in enumerate(sym_classes.values(), start=1):
            detail_rows.append({
                "instance": name,
                "p": p,
                "radius_used": radius,
                "class_id": idx,
                "class_size": len(centers),
                "centers": " ".join(map(str, centers))
            })

    # -------- Write Summary CSV
    summary_path = f"{args.out_prefix}_summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    # -------- Write Detail CSV
    detail_path = f"{args.out_prefix}_details.csv"
    with open(detail_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=detail_rows[0].keys())
        writer.writeheader()
        writer.writerows(detail_rows)

    print(f"[DONE] Written:")
    print(f"  {summary_path}")
    print(f"  {detail_path}")


if __name__ == "__main__":
    main()
