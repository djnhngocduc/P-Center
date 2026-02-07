import argparse
import json
from collections import defaultdict
from openpyxl import Workbook

EPS = 1e-12


# --------------------------------------------------
# Compute coverage symmetry
# --------------------------------------------------

def compute_coverage_symmetry(inst, radius):
    Nc, Nd, enabled_centers, demands, _ = inst.compute_reduction(radius)

    candidates = sorted((enabled_centers - Nc) - Nd)

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
    """
    Use the SAME pipeline as experiments.py
    """
    from experiments import load_instance

    # load_instance returns (inst, seed_idx)
    inst, seed_idx = load_instance(inst_desc)

    # use the mapped radius from inst.radii
    radius = inst.radii[seed_idx]

    return inst, radius


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", required=True)
    parser.add_argument("--output", default="symmetry_stats.xlsx")
    args = parser.parse_args()

    with open(args.instances, "r", encoding="utf-8") as f:
        data = json.load(f)

    wb = Workbook()
    summary_ws = wb.active
    summary_ws.title = "Summary"

    summary_ws.append([
        "Instance",
        "p",
        "Radius_used",
        "#Candidates",
        "#SymmetryClasses",
        "TotalSymmetricNodes",
        "MaxClassSize",
        "SymmetryRatio"
    ])

    for inst_desc in data:
        name = inst_desc["name"]
        p = inst_desc["p"]

        print(f"[PROCESS] {name}, p={p}")

        # ---- identical pipeline to experiments
        inst, radius = load_instance_and_seed_radius(inst_desc)

        print(f"  mapped radius = {radius}")

        sym_classes, Nc, Nd, candidates = compute_coverage_symmetry(
            inst, radius
        )

        total_sym_nodes = sum(len(v) for v in sym_classes.values())
        max_class = max((len(v) for v in sym_classes.values()), default=0)
        num_classes = len(sym_classes)
        num_candidates = len(candidates)

        symmetry_ratio = (
            total_sym_nodes / num_candidates
            if num_candidates > 0 else 0
        )

        # ---- Summary row
        summary_ws.append([
            name,
            p,
            radius,
            num_candidates,
            num_classes,
            total_sym_nodes,
            max_class,
            symmetry_ratio
        ])

        # ---- Detail sheet
        sheet_name = f"{name}_p{p}"
        ws = wb.create_sheet(title=sheet_name[:31])

        ws.append(["Instance", name])
        ws.append(["p", p])
        ws.append(["Radius_used", radius])
        ws.append(["|Nc|", len(Nc)])
        ws.append(["|Nd|", len(Nd)])
        ws.append(["#Candidates", num_candidates])
        ws.append(["#SymmetryClasses", num_classes])
        ws.append(["TotalSymmetricNodes", total_sym_nodes])
        ws.append(["MaxClassSize", max_class])
        ws.append([])

        ws.append(["ClassID", "ClassSize", "Centers"])

        for idx, centers in enumerate(sym_classes.values(), start=1):
            ws.append([
                idx,
                len(centers),
                ",".join(map(str, centers))
            ])

    wb.save(args.output)
    print(f"[DONE] Excel written to {args.output}")


if __name__ == "__main__":
    main()
