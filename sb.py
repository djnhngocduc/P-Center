from collections import defaultdict

EPS = 1e-12

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


def compute_automorphism_symmetry(inst, radius, Nc, Nd, enabled_centers, demands):
    import pynauty

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

# ============================================================
# Symmetry breaking using automorphism orbits
# ============================================================

def automorphism_symmetry_breaking(
    self,
    cnf,
    radius: float,
    Nc, Nd, enabled_centers, demands,
    mode: str = "chain",   # "chain" or "leader"
):
    """
    Break FULL structural automorphism symmetry of the coverage graph.
    No dominance. No identical-only logic.
    """

    ylit = self.y_lit_all

    orbits = compute_automorphism_symmetry(self, radius, Nc, Nd, enabled_centers, demands)

    if not orbits:
        return

    for group in orbits:
        group.sort()

        if mode == "leader":
            leader = group[0]
            for c in group[1:]:
                cnf.append([-ylit[c], ylit[leader]])

        else:  # chain mode
            for a, b in zip(group, group[1:]):
                cnf.append([-ylit[b], ylit[a]])