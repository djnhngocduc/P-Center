# sb.py
from __future__ import annotations
from typing import Dict, List, Set

try:
    import pynauty
except Exception:
    pynauty = None


def sb_coverage_order(
    self,
    cnf,
    radius: float,
    candidates: List[int],
    demands: List[int],
    Nc: Set[int],
) -> None:
    """
    STRONG & SAFE dominance-based ordering:
    enforce ordering only when coverage inclusion holds.

    NOTE: Orbit-SB is NOT done here anymore.
    """
    eps = 1e-12

    active_demands = [
        u for u in demands
        if not any(self.dist[c][u] <= radius + eps for c in Nc)
    ]

    cover: Dict[int, Set[int]] = {}
    for c in candidates:
        cover[c] = {u for u in active_demands if self.dist[c][u] <= radius + eps}

    # dominance-based ordering
    for i, ci in enumerate(candidates):
        for cj in candidates[i + 1:]:
            if cover[ci] >= cover[cj] and cover[ci] != cover[cj]:
                cnf.append([-self._y(cj), self._y(ci)])
            elif cover[cj] >= cover[ci] and cover[cj] != cover[ci]:
                cnf.append([-self._y(ci), self._y(cj)])


def orbit_sb_from_cnf(
    self,
    cnf,
    y_vars: List[int],
    *,
    orbit_mode: str = "chain",   # "chain" or "leader"
) -> None:
    """
    Sound orbit-based SB computed on the FINAL CNF (including aux vars).
    Graph construction (CNF -> colored graph):
      - For each variable v in 1..nv: two literal nodes v+ and v-
      - One clause node per clause
      - Clause node connected to its literal nodes
      - Mate edge v+--v- for each v
      - Color classes: pos-lits, neg-lits, clause nodes

    We then restrict SB ONLY to y_vars based on orbit of their positive literal node.
    """
    if pynauty is None:
        return
    if not y_vars or len(y_vars) <= 1:
        return

    # Get nv robustly
    nv = getattr(cnf, "nv", 0)
    if nv <= 0:
        mx = 0
        for cl in cnf.clauses:
            for lit in cl:
                a = abs(lit)
                if a > mx:
                    mx = a
        nv = mx
        if nv <= 0:
            return

    # literal nodes: 0..2*nv-1
    # pos(v) = 2*(v-1), neg(v)=2*(v-1)+1
    def pos_node(v: int) -> int:
        return 2 * (v - 1)

    def neg_node(v: int) -> int:
        return 2 * (v - 1) + 1

    m = len(cnf.clauses)
    base_clause = 2 * nv
    V = base_clause + m

    adj: Dict[int, Set[int]] = {i: set() for i in range(V)}

    # mate edges
    for v in range(1, nv + 1):
        a = pos_node(v)
        b = neg_node(v)
        adj[a].add(b)
        adj[b].add(a)

    # clause-literal incidence
    for ci, cl in enumerate(cnf.clauses):
        cnode = base_clause + ci
        for lit in cl:
            v = abs(lit)
            if v < 1 or v > nv:
                continue
            lnode = pos_node(v) if lit > 0 else neg_node(v)
            adj[cnode].add(lnode)
            adj[lnode].add(cnode)

    adjacency_dict = {v: sorted(list(nbs)) for v, nbs in adj.items()}

    # color classes
    pos_lits = set(range(0, 2 * nv, 2))
    neg_lits = set(range(1, 2 * nv, 2))
    clauses = set(range(base_clause, V))
    vertex_coloring = [pos_lits, neg_lits, clauses]

    try:
        g = pynauty.Graph(
            number_of_vertices=V,
            directed=False,
            adjacency_dict=adjacency_dict,
            vertex_coloring=vertex_coloring,
        )

        orbits = None

        # Prefer direct orbits() if available (most reliable across builds)
        if hasattr(pynauty, "orbits"):
            tmp = pynauty.orbits(g)
            if isinstance(tmp, list) and len(tmp) == V and all(isinstance(x, int) for x in tmp):
                orbits = tmp

        # Fallback: autgrp(g)[3] if available
        if orbits is None and hasattr(pynauty, "autgrp"):
            res = pynauty.autgrp(g)
            if isinstance(res, tuple) and len(res) >= 4:
                cand = res[3]
                if isinstance(cand, list) and len(cand) == V and all(isinstance(x, int) for x in cand):
                    orbits = cand

        if orbits is None:
            return

        # Group y-vars by orbit id of their positive literal node
        groups: Dict[int, List[int]] = {}
        for y in y_vars:
            if 1 <= y <= nv:
                oid = orbits[pos_node(y)]
                groups.setdefault(oid, []).append(y)

        # Add SB inside each orbit group
        for ys in groups.values():
            ys = sorted(set(ys))
            if len(ys) < 2:
                continue

            if orbit_mode == "leader":
                leader = ys[0]
                for y in ys[1:]:
                    cnf.append([-y, leader])
            else:
                for a, b in zip(ys, ys[1:]):
                    cnf.append([-b, a])

    except Exception:
        return
