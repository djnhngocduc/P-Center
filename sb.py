# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional

try:
    import pynauty
except Exception:
    pynauty = None


def sb_coverage_order(self, cnf, radius: float, candidates: List[int], demands: List[int], Nc: Set[int]) -> None:
    """
    Dominance-based ordering only (safe).
    """
    eps = 1e-12

    active_demands = [
        u for u in demands
        if not any(self.dist[c][u] <= radius + eps for c in Nc)
    ]

    cover: Dict[int, Set[int]] = {}
    for c in candidates:
        cover[c] = {u for u in active_demands if self.dist[c][u] <= radius + eps}

    for i, ci in enumerate(candidates):
        for cj in candidates[i + 1:]:
            if cover[ci] >= cover[cj] and cover[ci] != cover[cj]:
                cnf.append([-self._y(cj), self._y(ci)])
            elif cover[cj] >= cover[ci] and cover[cj] != cover[ci]:
                cnf.append([-self._y(ci), self._y(cj)])


def orbit_sb_y_only(
    cnf_y_clauses: List[List[int]],
    y_vars: List[int],
    *,
    orbit_mode: str = "chain",   # "chain" or "leader"
) -> List[List[int]]:
    """
    Compute orbits using ONLY clauses over y-vars.
    Returns list of orbit groups (each is list of y var ids).
    """
    if pynauty is None:
        return []
    if len(y_vars) <= 1:
        return []

    y_set = set(y_vars)

    # filter clauses to those using only y-vars (and non-empty)
    clauses = []
    for cl in cnf_y_clauses:
        if not cl:
            continue
        ok = True
        for lit in cl:
            if abs(lit) not in y_set:
                ok = False
                break
        if ok:
            clauses.append(cl)

    nv = len(y_vars)
    # map y var id -> 1..nv contiguous
    y_sorted = sorted(y_vars)
    vid = {y: i + 1 for i, y in enumerate(y_sorted)}  # 1..nv

    # Graph reduction (literal/clause incidence) but ONLY for these nv vars
    # literal nodes: 0..2*nv-1
    def pos_node(i: int) -> int:  # i in 1..nv
        return 2 * (i - 1)

    def neg_node(i: int) -> int:
        return 2 * (i - 1) + 1

    m = len(clauses)
    base_clause = 2 * nv
    V = base_clause + m

    adj: Dict[int, Set[int]] = {i: set() for i in range(V)}

    # mate edges
    for i in range(1, nv + 1):
        a = pos_node(i)
        b = neg_node(i)
        adj[a].add(b)
        adj[b].add(a)

    # clause-literal edges
    for ci, cl in enumerate(clauses):
        cnode = base_clause + ci
        for lit in cl:
            y = abs(lit)
            i = vid[y]
            lnode = pos_node(i) if lit > 0 else neg_node(i)
            adj[cnode].add(lnode)
            adj[lnode].add(cnode)

    adjacency_dict = {v: sorted(list(nbs)) for v, nbs in adj.items()}

    pos_lits = set(range(0, 2 * nv, 2))
    neg_lits = set(range(1, 2 * nv, 2))
    clause_nodes = set(range(base_clause, V))
    vertex_coloring = [pos_lits, neg_lits, clause_nodes]

    try:
        g = pynauty.Graph(
            number_of_vertices=V,
            directed=False,
            adjacency_dict=adjacency_dict,
            vertex_coloring=vertex_coloring,
        )

        orbits = None
        if hasattr(pynauty, "autgrp"):
            res = pynauty.autgrp(g)
            if isinstance(res, tuple) and len(res) >= 4:
                cand = res[3]
                if isinstance(cand, list) and len(cand) == V and all(isinstance(x, int) for x in cand):
                    orbits = cand
            if orbits is None and isinstance(res, tuple):
                for item in res:
                    if isinstance(item, list) and len(item) == V and all(isinstance(x, int) for x in item):
                        orbits = item
                        break

        if orbits is None and hasattr(pynauty, "orbits"):
            tmp = pynauty.orbits(g)
            if isinstance(tmp, list) and len(tmp) == V and all(isinstance(x, int) for x in tmp):
                orbits = tmp

        if orbits is None:
            return []

        # group by orbit id of positive literal node
        groups: Dict[int, List[int]] = {}
        for y in y_sorted:
            i = vid[y]
            oid = orbits[pos_node(i)]
            groups.setdefault(oid, []).append(y)

        out = []
        for ys in groups.values():
            if len(ys) >= 2:
                out.append(sorted(ys))
        out.sort(key=lambda lst: (len(lst), lst[0]))
        return out

    except Exception:
        return []


def add_orbit_sb_clauses(cnf, orbit_groups: List[List[int]], *, orbit_mode: str = "chain") -> None:
    """
    Add SB constraints on y-vars for each orbit group.
    """
    for ys in orbit_groups:
        if len(ys) < 2:
            continue
        ys = sorted(ys)
        if orbit_mode == "leader":
            leader = ys[0]
            for y in ys[1:]:
                cnf.append([-y, leader])
        else:
            for a, b in zip(ys, ys[1:]):
                cnf.append([-b, a])
