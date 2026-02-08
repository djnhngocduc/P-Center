import argparse
import json
import csv
from collections import defaultdict

EPS = 1e-12


# --------------------------------------------------
# IDENTICAL COVERAGE SYMMETRY
# --------------------------------------------------

def compute_identical_coverage_symmetry(inst, radius):
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

    sym_classes = [
        sorted(v) for v in signature.values()
        if len(v) >= 2
    ]

    return sym_classes, Nc, Nd, candidates


# --------------------------------------------------
# AUTOMORPHISM SYMMETRY (pynauty)
# --------------------------------------------------

def _extract_orbit_ids_from_autgrp_result(result):
    """
    pynauty.autgrp(g) thường trả tuple dạng:
      (generators, grpsize1, grpsize2, orbits, numorbits)
    => phần tử cuối (numorbits) là int, còn 'orbits' (list/iterable) nằm ở -2.

    Nhưng để an toàn theo nhiều version/fork:
    - Nếu result[-1] là int -> lấy result[-2]
    - Nếu result[-1] iterable -> lấy result[-1]
    - Nếu vẫn không đúng -> raise để dễ debug.
    """
    if not isinstance(result, (tuple, list)) or len(result) == 0:
        raise TypeError(f"autgrp returned unexpected type/value: {type(result)} {result!r}")

    last = result[-1]

    # Trường hợp phổ biến: last là numorbits (int)
    if isinstance(last, int):
        if len(result) < 2:
            raise TypeError(f"autgrp returned only an int: {result!r}")
        orbit_ids = result[-2]
    else:
        orbit_ids = last

    # orbit_ids phải iterable (list/tuple)
    if isinstance(orbit_ids, int):
        raise TypeError(
            "orbit_ids is still int after extraction. "
            f"autgrp returned: {result!r}"
        )

    # Một số version trả orbits dạng list length = nVertices,
    # mỗi phần tử là orbit-id (int)
    return orbit_ids


def compute_automorphism_symmetry(inst, radius):
    import pynauty

    Nc, Nd, enabled_centers, demands, _ = inst.compute_reduction(radius)
    candidates = sorted(enabled_centers - Nc - Nd)

    nC = len(candidates)
    nD = len(demands)
    total_vertices = nC + nD

    if total_vertices == 0:
        return []

    center_index = {c: i for i, c in enumerate(candidates)}
    demand_index = {u: i for i, u in enumerate(demands)}

    adj = {i: [] for i in range(total_vertices)}

    for c in candidates:
        ci = center_index[c]
        row = inst.dist[c]

        for u in demands:
            if row[u] <= radius + EPS:
                ui = demand_index[u] + nC
                adj[ci].append(ui)
                adj[ui].append(ci)

    coloring = [
        set(range(0, nC)),
        set(range(nC, nC + nD))
    ]

    g = pynauty.Graph(
        number_of_vertices=total_vertices,
        adjacency_dict=adj,
        vertex_coloring=coloring
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
# Load instance
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

        # =============================
        # IDENTICAL
        # =============================

        identical_classes, Nc, Nd, candidates = \
            compute_identical_coverage_symmetry(inst, radius)

        num_candidates = len(candidates)
        total_sym_nodes = sum(len(c) for c in identical_classes)
        max_class = max((len(c) for c in identical_classes), default=0)
        num_classes = len(identical_classes)

        symmetry_ratio = (
            total_sym_nodes / num_candidates
            if num_candidates > 0 else 0
        )

        summary_rows.append({
            "instance": name,
            "p": p,
            "radius_used": radius,
            "symmetry_type": "identical",
            "num_candidates": num_candidates,
            "num_symmetry_classes": num_classes,
            "total_symmetric_nodes": total_sym_nodes,
            "max_class_size": max_class,
            "symmetry_ratio": symmetry_ratio
        })

        for idx, centers in enumerate(identical_classes, start=1):
            detail_rows.append({
                "instance": name,
                "p": p,
                "radius_used": radius,
                "symmetry_type": "identical",
                "class_id": idx,
                "class_size": len(centers),
                "centers": " ".join(map(str, centers))
            })

        # =============================
        # AUTOMORPHISM
        # =============================

        auto_classes = compute_automorphism_symmetry(inst, radius)

        total_auto_nodes = sum(len(c) for c in auto_classes)
        max_auto_class = max((len(c) for c in auto_classes), default=0)
        num_auto_classes = len(auto_classes)

        auto_ratio = (
            total_auto_nodes / num_candidates
            if num_candidates > 0 else 0
        )

        summary_rows.append({
            "instance": name,
            "p": p,
            "radius_used": radius,
            "symmetry_type": "automorphism",
            "num_candidates": num_candidates,
            "num_symmetry_classes": num_auto_classes,
            "total_symmetric_nodes": total_auto_nodes,
            "max_class_size": max_auto_class,
            "symmetry_ratio": auto_ratio
        })

        for idx, centers in enumerate(auto_classes, start=1):
            detail_rows.append({
                "instance": name,
                "p": p,
                "radius_used": radius,
                "symmetry_type": "automorphism",
                "class_id": idx,
                "class_size": len(centers),
                "centers": " ".join(map(str, centers))
            })

    # ------------------------------
    # Write CSV
    # ------------------------------

    summary_path = f"{args.out_prefix}_summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    detail_path = f"{args.out_prefix}_details.csv"
    with open(detail_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=detail_rows[0].keys())
        writer.writeheader()
        writer.writerows(detail_rows)

    print("[DONE]")
    print("  ", summary_path)
    print("  ", detail_path)


if __name__ == "__main__":
    main()
