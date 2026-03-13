import argparse
import csv
import json
from collections import defaultdict

EPS = 1e-12


# --------------------------------------------------
# Helpers: parse pynauty.autgrp(...)
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


def _extract_generators_from_autgrp_result(result, total_vertices):
    """
    Expect result[0] to contain permutations.
    Keep only objects that are valid permutations of 0..total_vertices-1.
    """
    if not isinstance(result, (tuple, list)) or len(result) == 0:
        raise TypeError(f"autgrp returned unexpected: {type(result)} {result!r}")

    first = result[0]
    if not isinstance(first, (tuple, list)):
        raise TypeError(f"autgrp generators malformed: {result!r}")

    gens = []
    for g in first:
        if not isinstance(g, (tuple, list)):
            continue
        if len(g) != total_vertices:
            continue

        perm = [int(x) for x in g]
        if sorted(perm) != list(range(total_vertices)):
            continue
        gens.append(perm)

    return gens


# --------------------------------------------------
# Build current coverage graph
# --------------------------------------------------

def _build_coverage_graph(inst, radius, candidates, active_demands):
    candidates = sorted(candidates)
    active_demands = sorted(active_demands)

    nC = len(candidates)
    nD = len(active_demands)
    total_vertices = nC + nD

    center_index = {c: i for i, c in enumerate(candidates)}
    demand_index = {u: i for i, u in enumerate(active_demands)}

    adj = {i: [] for i in range(total_vertices)}

    for c in candidates:
        ci = center_index[c]
        row = inst.dist[c]
        for u in active_demands:
            if row[u] <= radius + EPS:
                ui = nC + demand_index[u]
                adj[ci].append(ui)
                adj[ui].append(ci)

    coloring = []
    if nC > 0:
        coloring.append(set(range(0, nC)))
    if nD > 0:
        coloring.append(set(range(nC, nC + nD)))

    return candidates, active_demands, nC, nD, total_vertices, adj, coloring


# --------------------------------------------------
# Orbit classes + center generators
# --------------------------------------------------

def compute_center_generators(
    inst,
    radius,
    candidates,
    active_demands,
):
    import pynauty

    candidates, active_demands, nC, nD, total_vertices, adj, coloring = _build_coverage_graph(
        inst, radius, candidates, active_demands
    )

    if total_vertices == 0:
        return candidates, [], []

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

    raw_gens = _extract_generators_from_autgrp_result(result, total_vertices)

    center_gens = []
    seen = set()

    for perm in raw_gens:
        cperm = []
        ok = True
        for i in range(nC):
            img = perm[i]
            if not (0 <= img < nC):
                ok = False
                break
            cperm.append(img)

        if not ok:
            continue

        tpl = tuple(cperm)
        if tpl == tuple(range(nC)):
            continue
        if tpl in seen:
            continue

        seen.add(tpl)
        center_gens.append(cperm)

    if nC > 0 and not center_gens and auto_classes:
        raise RuntimeError(
            "No nontrivial center generators parsed from pynauty.autgrp(...), "
            "but nontrivial center orbits exist. Check autgrp output format."
        )

    return candidates, auto_classes, center_gens


# --------------------------------------------------
# Lex-leader encoding
# --------------------------------------------------

def _new_var(cnf):
    v = max(getattr(cnf, "nv", 0), 0) + 1
    cnf.nv = v
    return v


def _add_equiv_same_bit(cnf, z, x, y):
    # z <-> (x == y)
    cnf.append([-z, -x, y])
    cnf.append([-z, x, -y])
    cnf.append([-x, -y, z])
    cnf.append([x, y, z])


def _add_equiv_prefix_and_same(cnf, z, prev_eq, x, y):
    # z <-> (prev_eq AND (x == y))
    cnf.append([-z, prev_eq])
    cnf.append([-z, -x, y])
    cnf.append([-z, x, -y])
    cnf.append([-prev_eq, -x, -y, z])
    cnf.append([-prev_eq, x, y, z])


def add_lex_leader_geq_clause_system(cnf, xs, ys):
    """
    Add CNF for:
        xs >=_lex ys
    with Boolean order 1 > 0.
    """
    m = len(xs)
    if m != len(ys):
        raise ValueError("xs and ys must have same length")
    if m == 0:
        return

    eq = []

    for i in range(m):
        x = xs[i]
        y = ys[i]

        if i == 0:
            cnf.append([x, -y])
        else:
            cnf.append([-eq[i - 1], x, -y])

        if i < m - 1:
            z = _new_var(cnf)
            if i == 0:
                _add_equiv_same_bit(cnf, z, x, y)
            else:
                _add_equiv_prefix_and_same(cnf, z, eq[i - 1], x, y)
            eq.append(z)


def add_lex_leader_for_generator(
    inst,
    cnf,
    candidates_sorted,
    center_perm,
):
    """
    Add:
        X >=_lex pi(X)

    where X is the y-vector over candidates_sorted.
    """
    nC = len(candidates_sorted)
    if nC == 0:
        return

    inv = [0] * nC
    for i, j in enumerate(center_perm):
        inv[j] = i

    xs = [inst.y_lit_all[c] for c in candidates_sorted]
    ys = [inst.y_lit_all[candidates_sorted[inv_pos]] for inv_pos in inv]

    add_lex_leader_geq_clause_system(cnf, xs, ys)


# --------------------------------------------------
# Public API: full lex-leader only
# --------------------------------------------------

def automorphism_full_lex_leader_symmetry_breaking(
    inst,
    cnf,
    radius,
    candidates,
    active_demands,
    verbose=False,
):
    candidates_sorted, orbits, center_gens = compute_center_generators(
        inst=inst,
        radius=radius,
        candidates=candidates,
        active_demands=active_demands,
    )

    if verbose:
        print(
            f"[SB-FULL] generators={len(center_gens)} "
            f"orbits={len(orbits)} candidates={len(candidates_sorted)}"
        )

    for gen in center_gens:
        add_lex_leader_for_generator(
            inst=inst,
            cnf=cnf,
            candidates_sorted=candidates_sorted,
            center_perm=gen,
        )


# --------------------------------------------------
# Analysis CLI
# --------------------------------------------------

def load_instance_and_seed_radius(inst_desc):
    from experiments import load_instance
    inst, seed_idx = load_instance(inst_desc)
    radius = inst.radii[seed_idx]
    return inst, radius


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", required=True)
    parser.add_argument("--out-prefix", default="symmetry_full")
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

        candidates = list(range(inst.n))
        active_demands = list(range(inst.n))

        candidates_sorted, auto_classes, center_gens = compute_center_generators(
            inst=inst,
            radius=radius,
            candidates=candidates,
            active_demands=active_demands,
        )

        total_auto_nodes = sum(len(c) for c in auto_classes)
        max_auto_class = max((len(c) for c in auto_classes), default=0)
        num_auto_classes = len(auto_classes)
        auto_ratio = (
            (total_auto_nodes / len(candidates_sorted))
            if len(candidates_sorted) > 0 else 0.0
        )

        summary_rows.append({
            "instance": name,
            "p": p,
            "radius_used": radius,
            "symmetry_type": "automorphism_full_lex",
            "num_candidates": len(candidates_sorted),
            "num_demands": len(active_demands),
            "num_symmetry_classes": num_auto_classes,
            "total_symmetric_nodes": total_auto_nodes,
            "max_class_size": max_auto_class,
            "num_generators": len(center_gens),
            "symmetry_ratio": auto_ratio,
        })

        for idx, centers in enumerate(auto_classes, start=1):
            detail_rows.append({
                "instance": name,
                "p": p,
                "radius_used": radius,
                "symmetry_type": "automorphism_full_lex",
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