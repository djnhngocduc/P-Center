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


def _transitive_reduce_edges(
    nodes: List[int],
    out_edges: Dict[int, Set[int]],
) -> Dict[int, Set[int]]:
    """
    Given a directed graph out_edges (u -> set(v)), remove redundant edges:
      remove u->v if there exists a path u => v without using that direct edge.

    n<=150 so O(E*(V+E)) approach is fine and stable.
    """
    reduced = {u: set(vs) for u, vs in out_edges.items()}

    # Small helper BFS reachability
    def reachable(src: int, target: int, banned_edge: Tuple[int, int]) -> bool:
        bu, bv = banned_edge
        q = [src]
        seen = {src}
        while q:
            u = q.pop()
            for w in reduced.get(u, ()):
                if u == bu and w == bv:
                    continue
                if w == target:
                    return True
                if w not in seen:
                    seen.add(w)
                    q.append(w)
        return False

    # Check each edge for redundancy
    for u in nodes:
        outs = list(reduced.get(u, ()))
        for v in outs:
            # Temporarily test if v reachable without using (u,v)
            if reachable(u, v, (u, v)):
                reduced[u].discard(v)

    return reduced


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
        (we use the *reduced* implication edges that we actually add to CNF)

    This cannot remove valid solutions => never SAT->UNSAT.
    """

    eps = 1e-12
    pair_clauses = pair_clauses or []

    if len(candidates) <= 1:
        return

    # -----------------------
    # active demands (not already covered by Nc)
    # -----------------------
    active_demands = [
        u for u in demands
        if not any(self.dist[c][u] <= radius + eps for c in Nc)
    ]

    # -----------------------
    # A1: build coverage as BITSET mask over active_demands
    #   demand_index[u] in [0..k-1], cover_mask[c] is int bitmask
    # -----------------------
    k = len(active_demands)
    demand_index = {u: idx for idx, u in enumerate(active_demands)}

    cover_mask: Dict[int, int] = {}
    cover_size: Dict[int, int] = {}

    for c in candidates:
        m = 0
        row = self.dist[c]
        # iterate active_demands only (k <= n)
        for u in active_demands:
            if row[u] <= radius + eps:
                m |= (1 << demand_index[u])
        cover_mask[c] = m
        cover_size[c] = m.bit_count()

    # -----------------------
    # A2: bucketing by |cover| (sort descending)
    # Only compare i->j when size_i >= size_j
    # -----------------------
    cand_sorted = sorted(candidates, key=lambda c: (-cover_size[c], c))

    # Build implication graph edges (dominance):
    # If cover(ci) ⊋ cover(cj): add edge (cj -> ci) meaning y[cj] -> y[ci]
    out_edges: Dict[int, Set[int]] = {c: set() for c in candidates}

    # Dominance test using bitset:
    # mi dominates mj iff (mi | mj) == mi and mi != mj
    for idx_i, ci in enumerate(cand_sorted):
        mi = cover_mask[ci]
        si = cover_size[ci]

        # Only need compare with later elements (with <= size)
        for cj in cand_sorted[idx_i + 1:]:
            sj = cover_size[cj]
            if sj > si:
                # due to sort, shouldn't happen, but keep safe
                continue

            mj = cover_mask[cj]

            # If sizes equal but masks differ, neither can strictly dominate the other.
            # (Because strict dominance implies strictly larger popcount.)
            if si == sj:
                continue

            # ci dominates cj?
            if (mi | mj) == mi:
                # cj -> ci
                out_edges[cj].add(ci)
            # else: cannot have cj dominates ci because sj < si

    # -----------------------
    # A3: transitive reduction on implication graph
    # This reduces CNF size but keeps same logical strength (since implications are transitive)
    # -----------------------
    reduced_out = _transitive_reduce_edges(candidates, out_edges)

    # Emit CNF clauses for reduced edges
    # and record dom_out/dom_in for orbit-safe signature
    dom_out: Dict[int, Set[int]] = {c: set() for c in candidates}
    dom_in: Dict[int, Set[int]] = {c: set() for c in candidates}

    for u in candidates:
        for v in sorted(reduced_out.get(u, ())):
            # u -> v  => (-y[u] OR y[v])
            cnf.append([-self._y(u), self._y(v)])
            dom_out[u].add(v)
            dom_in[v].add(u)

    # -----------------------
    # Orbit-safe equivalence SB (your existing idea)
    # -----------------------
    if not enable_orbit_sb:
        return

    # pair-clauses adjacency among candidates
    pair_adj: Dict[int, Set[int]] = {c: set() for c in candidates}
    cand_set = set(candidates)
    for a, b in pair_clauses:
        if a in cand_set and b in cand_set:
            pair_adj[a].add(b)
            pair_adj[b].add(a)

    # Signature:
    # - exact same cover mask
    # - exact same pair neighbors
    # - exact same reduced dominance out/in sets
    sig2group: Dict[Tuple, List[int]] = {}

    for c in candidates:
        sig = (
            cover_mask[c],               # EXACT identical coverage on active demands (fast & safe)
            frozenset(pair_adj[c]),      # identical pair adjacency
            frozenset(dom_out[c]),       # identical reduced dominance outgoing
            frozenset(dom_in[c]),        # identical reduced dominance incoming
        )
        sig2group.setdefault(sig, []).append(c)

    for g in sig2group.values():
        if len(g) < 2:
            continue
        if orbit_mode == "leader":
            _add_leader_sb(self, cnf, g)
        else:
            _add_chain_sb(self, cnf, g)
