# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional


def _add_chain_sb(self, cnf, group: List[int], y_lit_all: List[int]) -> None:
    if len(group) < 2:
        return
    g = sorted(group)
    for a, b in zip(g, g[1:]):
        cnf.append([-y_lit_all[b], y_lit_all[a]])


def _add_leader_sb(self, cnf, group: List[int], y_lit_all: List[int]) -> None:
    if len(group) < 2:
        return
    g = sorted(group)
    leader = g[0]
    yl = y_lit_all[leader]
    for j in g[1:]:
        cnf.append([-y_lit_all[j], yl])


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
    eps = 1e-12
    pair_clauses = pair_clauses or []

    if len(candidates) <= 1:
        return

    dist = self.dist
    rad = radius + eps

    # Fast literal lookup (indices 0..n-1)
    y = self._y
    y_lit_all = [y(i) for i in range(self.n)]

    # -----------------------
    # Active demands (not covered by Nc): break early + cache Nc rows
    # -----------------------
    if Nc:
        Nc_rows = [dist[c] for c in Nc]
        active_demands: List[int] = []
        for u in demands:
            for rowc in Nc_rows:
                if rowc[u] <= rad:
                    break
            else:
                active_demands.append(u)
    else:
        active_demands = list(demands)

    k = len(active_demands)

    # -----------------------
    # A1: coverage bitmasks + group by mask
    # -----------------------
    cover_mask: Dict[int, int] = {}
    mask2cands: Dict[int, List[int]] = {}

    if k > 0:
        act = active_demands
        for c in candidates:
            row = dist[c]
            m = 0
            for bit, u in enumerate(act):
                if row[u] <= rad:
                    m |= (1 << bit)
            cover_mask[c] = m
            mask2cands.setdefault(m, []).append(c)
    else:
        for c in candidates:
            cover_mask[c] = 0
        mask2cands = {0: list(candidates)}

    for lst in mask2cands.values():
        lst.sort()

    # -----------------------
    # A2: bucket masks by popcount
    # -----------------------
    size2masks: Dict[int, List[int]] = {}
    for m in mask2cands.keys():
        size2masks.setdefault(m.bit_count(), []).append(m)

    sizes = sorted(size2masks.keys(), reverse=True)
    for s in sizes:
        size2masks[s].sort()

    dom_out: Dict[int, Set[int]] = {c: set() for c in candidates}
    dom_in: Dict[int, Set[int]] = {c: set() for c in candidates}

    # -----------------------
    # FULL dominance ordering (sound)
    # -----------------------
    if k > 0 and len(sizes) > 1:
        ALL = (1 << k) - 1
        for i, s_hi in enumerate(sizes):
            hi_masks = size2masks[s_hi]
            for s_lo in sizes[i + 1:]:
                lo_masks = size2masks[s_lo]

                for mi in hi_masks:
                    hi_nodes = mask2cands[mi]
                    nmi = ALL ^ mi  # complement in k bits

                    for mj in lo_masks:
                        # mj ⊆ mi ?
                        if (mj & nmi) != 0:
                            continue

                        lo_nodes = mask2cands[mj]

                        for cj in lo_nodes:
                            out_cj = dom_out[cj]
                            ycj = -y_lit_all[cj]
                            for ci in hi_nodes:
                                if ci in out_cj:
                                    continue
                                cnf.append([ycj, y_lit_all[ci]])
                                out_cj.add(ci)
                                dom_in[ci].add(cj)

    # -----------------------
    # orbit-safe equivalence SB (level-1 optimized, ALWAYS tries to add if groups exist)
    #  - reuse mask2cands (no extra mask_groups dict)
    #  - freeze sets once (avoid frozenset(...) per signature creation)
    # -----------------------
    if not enable_orbit_sb or len(candidates) <= 1:
        return

    # pair adjacency among candidates
    pair_adj: Dict[int, Set[int]] = {c: set() for c in candidates}
    if pair_clauses:
        for a, b in pair_clauses:
            if a in pair_adj and b in pair_adj:
                pair_adj[a].add(b)
                pair_adj[b].add(a)

    # freeze once
    pair_adj_f = {c: frozenset(pair_adj[c]) for c in candidates}
    dom_out_f = {c: frozenset(dom_out[c]) for c in candidates}
    dom_in_f  = {c: frozenset(dom_in[c])  for c in candidates}

    # iterate per mask-group using existing mask2cands
    for _, group in mask2cands.items():
        if len(group) < 2:
            continue

        sig2group: Dict[Tuple[object, object, object], List[int]] = {}
        for c in group:
            sig = (pair_adj_f[c], dom_out_f[c], dom_in_f[c])
            sig2group.setdefault(sig, []).append(c)

        for g in sig2group.values():
            if len(g) < 2:
                continue
            if orbit_mode == "leader":
                _add_leader_sb(self, cnf, g, y_lit_all)
            else:
                _add_chain_sb(self, cnf, g, y_lit_all)
