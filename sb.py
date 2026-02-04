# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional


def _add_chain_sb(self, cnf, group: List[int]) -> None:
    """
    Enforce y[g0] >= y[g1] >= ... >= y[gk-1]
    CNF: (-y[g_{t+1}] OR y[g_t])
    """
    if len(group) < 2:
        return
    g = sorted(group)
    for a, b in zip(g, g[1:]):
        cnf.append([-self._y(b), self._y(a)])


def _add_leader_sb(self, cnf, group: List[int]) -> None:
    """
    Leader form: for all j in group\{leader}: y[j] -> y[leader]
    """
    if len(group) < 2:
        return
    g = sorted(group)
    leader = g[0]
    for j in g[1:]:
        cnf.append([-self._y(j), self._y(leader)])


def sb_coverage_order(
    self,
    cnf,
    radius: float,
    candidates: List[int],
    demands: List[int],
    Nc: Set[int],
    *,
    enable_orbit_sb: bool = True,
    orbit_mode: str = "chain",  # "chain" | "leader"
    pair_clauses: Optional[List[Tuple[int, int]]] = None,
) -> None:
    """
    Dominance-based ordering (sound) + Orbit-safe equivalence SB (sound).

    Orbit-safe equivalence SB groups facilities ONLY if they are truly interchangeable
    w.r.t. the constraints currently modeled:
      - identical coverage on active demands
      - identical neighborhood in pair-clauses graph (y(a) OR y(b))
      - identical dominance in/out pattern implied by cover inclusion

    This cannot remove valid solutions => never SAT->UNSAT.
    """

    eps = 1e-12
    pair_clauses = pair_clauses or []

    # active demands = chưa được cover bởi Nc
    active_demands = [
        u for u in demands
        if not any(self.dist[c][u] <= radius + eps for c in Nc)
    ]

    # coverage sets on active_demands
    cover: Dict[int, Set[int]] = {}
    for c in candidates:
        cover[c] = {u for u in active_demands if self.dist[c][u] <= radius + eps}

    # ---------- dominance-based ordering ----------
    # record dominance edges for orbit-safe signature
    dom_out: Dict[int, Set[int]] = {c: set() for c in candidates}  # c -> set of smaller covered facility forced before?
    dom_in: Dict[int, Set[int]] = {c: set() for c in candidates}

    for i, ci in enumerate(candidates):
        for cj in candidates[i + 1:]:
            si = cover[ci]
            sj = cover[cj]

            if si >= sj and si != sj:
                # cj -> ci   i.e. y[cj] -> y[ci]
                cnf.append([-self._y(cj), self._y(ci)])
                dom_out[cj].add(ci)
                dom_in[ci].add(cj)

            elif sj >= si and sj != si:
                # ci -> cj
                cnf.append([-self._y(ci), self._y(cj)])
                dom_out[ci].add(cj)
                dom_in[cj].add(ci)

    # ---------- orbit-safe equivalence SB ----------
    if not enable_orbit_sb:
        return
    if len(candidates) <= 1:
        return

    # pair-clauses adjacency among candidates
    pair_adj: Dict[int, Set[int]] = {c: set() for c in candidates}
    cand_set = set(candidates)
    for a, b in pair_clauses:
        if a in cand_set and b in cand_set:
            pair_adj[a].add(b)
            pair_adj[b].add(a)

    # Build signatures:
    # Important: we must ensure equality implies true interchangeability.
    # Using exact sets (frozenset) is safest (slower but n<=150 ok).
    sig2group: Dict[Tuple, List[int]] = {}

    for c in candidates:
        sig = (
            frozenset(cover[c]),          # identical coverage on active demands
            frozenset(pair_adj[c]),       # identical pair-clause neighbors
            frozenset(dom_out[c]),        # identical dominance outgoing
            frozenset(dom_in[c]),         # identical dominance incoming
        )
        sig2group.setdefault(sig, []).append(c)

    # Add SB per group (only groups size>=2)
    for g in sig2group.values():
        if len(g) < 2:
            continue
        if orbit_mode == "leader":
            _add_leader_sb(self, cnf, g)
        else:
            _add_chain_sb(self, cnf, g)
