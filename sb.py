import argparse
import json
import csv
from collections import defaultdict

EPS = 1e-12


# --------------------------------------------------
# Extract orbit ids safely from pynauty.autgrp
# --------------------------------------------------

def _extract_orbit_ids_from_autgrp_result(result):
    if not isinstance(result, (tuple, list)) or len(result) == 0:
        raise TypeError(f"autgrp returned unexpected: {type(result)} {result!r}")

    last = result[-1]

    if isinstance(last, int):
        if len(result) < 2:
            raise TypeError(f"autgrp malformed: {result!r}")
        orbit_ids = result[-2]
    else:
        orbit_ids = last

    if isinstance(orbit_ids, int):
        raise TypeError(f"orbit_ids malformed: {result!r}")

    return orbit_ids


# --------------------------------------------------
# CURRENT-CNF-FAITHFUL AUTOMORPHISM
# Graph = candidate centers + active demands
# edge(c,u) iff dist[c][u] <= radius
# --------------------------------------------------

def compute_automorphism_symmetry(
    inst,
    radius,
    candidates,
    active_demands,
):
    import pynauty

    candidates = sorted(candidates)
    active_demands = sorted(active_demands)

    nC = len(candidates)
    nD = len(active_demands)

    total_vertices = nC + nD
    if total_vertices == 0:
        return []

    center_index = {c: i for i, c in enumerate(candidates)}
    demand_index = {u: i for i, u in enumerate(active_demands)}

    adj = {i: [] for i in range(total_vertices)}

    # ----------------------------------
    # center -- demand edges
    # ----------------------------------
    for c in candidates:
        ci = center_index[c]
        row = inst.dist[c]

        for u in active_demands:
            if row[u] <= radius + EPS:
                ui = nC + demand_index[u]
                adj[ci].append(ui)
                adj[ui].append(ci)

    # ----------------------------------
    # Coloring: centers / demands
    # ----------------------------------
    coloring = []
    if nC > 0:
        coloring.append(set(range(0, nC)))
    if nD > 0:
        coloring.append(set(range(nC, nC + nD)))

    g = pynauty.Graph(
        number_of_vertices=total_vertices,
        adjacency_dict=adj,
        vertex_coloring=coloring,
    )

    result = pynauty.autgrp(g)
    orbit_ids = _extract_orbit_ids_from_autgrp_result(result)

    orbit_dict = defaultdict(list)
    for v, oid in enumerate(orbit_ids):
        orbit_dict[oid].append(v)

    auto_classes = []

    for orbit in orbit_dict.values():
        centers = [candidates[v] for v in orbit if v < nC]
        if len(centers) >= 2:
            auto_classes.append(sorted(centers))

    auto_classes.sort()
    return auto_classes


# --------------------------------------------------
# Symmetry breaking clauses
# --------------------------------------------------

def automorphism_symmetry_breaking(
    inst,
    cnf,
    radius,
    candidates,
    active_demands,
    mode="chain",
):
    ylit = inst.y_lit_all

    orbits = compute_automorphism_symmetry(
        inst=inst,
        radius=radius,
        candidates=candidates,
        active_demands=active_demands,
    )

    if not orbits:
        return

    for group in orbits:
        group = sorted(group)

        if mode == "leader":
            leader = group[0]
            for c in group[1:]:
                cnf.append([-ylit[c], ylit[leader]])
        else:  # chain
            for a, b in zip(group, group[1:]):
                cnf.append([-ylit[b], ylit[a]])


def load_instance_and_seed_radius(inst_desc):
    from experiments import load_instance
    inst, seed_idx = load_instance(inst_desc, use_seed_radius=True)
    radius = inst.radii[seed_idx]
    return inst, radius


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

        # match current encoder.py
        candidates = list(range(inst.n))
        active_demands = list(range(inst.n))
        num_candidates = len(candidates)

        auto_classes = compute_automorphism_symmetry(
            inst=inst,
            radius=radius,
            candidates=candidates,
            active_demands=active_demands,
        )

        total_auto_nodes = sum(len(c) for c in auto_classes)
        max_auto_class = max((len(c) for c in auto_classes), default=0)
        num_auto_classes = len(auto_classes)

        auto_ratio = (total_auto_nodes / num_candidates) if num_candidates > 0 else 0.0

        summary_rows.append({
            "instance": name,
            "p": p,
            "radius_used": radius,
            "symmetry_type": "automorphism",
            "num_candidates": num_candidates,
            "num_demands": len(active_demands),
            "num_symmetry_classes": num_auto_classes,
            "total_symmetric_nodes": total_auto_nodes,
            "max_class_size": max_auto_class,
            "symmetry_ratio": auto_ratio,
        })

        for idx, centers in enumerate(auto_classes, start=1):
            detail_rows.append({
                "instance": name,
                "p": p,
                "radius_used": radius,
                "symmetry_type": "automorphism",
                "class_id": idx,
                "class_size": len(centers),
                "centers": " ".join(map(str, centers)),
            })

    if summary_rows:
        summary_path = f"{args.out_prefix}_summary.csv"
        with open(summary_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            writer.writeheader()
            writer.writerows(summary_rows)
    else:
        summary_path = None

    if detail_rows:
        detail_path = f"{args.out_prefix}_details.csv"
        with open(detail_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(detail_rows[0].keys()))
            writer.writeheader()
            writer.writerows(detail_rows)
    else:
        detail_path = None

    print("[DONE]")
    if summary_path:
        print("  ", summary_path)
    if detail_path:
        print("  ", detail_path)


if __name__ == "__main__":
    main()