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
# FULL CNF-FAITHFUL AUTOMORPHISM
# --------------------------------------------------

def compute_automorphism_symmetry(
    inst,
    radius,
    candidates,
    active_demands,
    at_least_pairs
):
    import pynauty

    cand_set = set(candidates)

    nC = len(candidates)
    nD = len(active_demands)

    # keep only relevant pairs among candidates
    pair_list = [
        (a, b)
        for (a, b) in at_least_pairs
        if a in cand_set and b in cand_set
    ]

    nP = len(pair_list)

    total_vertices = nC + nD + nP
    if total_vertices == 0:
        return []

    center_index = {c: i for i, c in enumerate(candidates)}
    demand_index = {u: i for i, u in enumerate(active_demands)}

    adj = {i: [] for i in range(total_vertices)}

    # ----------------------------------
    # 1) center — demand edges
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
    # 2) pair-gadget — center edges
    # ----------------------------------
    baseP = nC + nD

    for k, (a, b) in enumerate(pair_list):
        pv = baseP + k
        ia = center_index[a]
        ib = center_index[b]

        adj[pv].append(ia)
        adj[ia].append(pv)

        adj[pv].append(ib)
        adj[ib].append(pv)

    # ----------------------------------
    # Coloring: centers / demands / pair-nodes
    # ----------------------------------
    coloring = []

    if nC > 0:
        coloring.append(set(range(0, nC)))
    if nD > 0:
        coloring.append(set(range(nC, nC + nD)))
    if nP > 0:
        coloring.append(set(range(nC + nD, total_vertices)))

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
    self,
    cnf,
    radius,
    candidates,
    active_demands,
    at_least_pairs,
    mode="chain",
):
    ylit = self.y_lit_all

    orbits = compute_automorphism_symmetry(
        self,
        radius,
        candidates,
        active_demands,
        at_least_pairs,
    )

    if not orbits:
        return

    for group in orbits:
        group.sort()

        if mode == "leader":
            leader = group[0]
            for c in group[1:]:
                cnf.append([-ylit[c], ylit[leader]])
        else:  # chain
            for a, b in zip(group, group[1:]):
                cnf.append([-ylit[b], ylit[a]])