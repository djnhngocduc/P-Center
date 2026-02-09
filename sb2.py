# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional


def _add_chain_sb(self, cnf, group: List[int], y_lit_all: List[int]) -> None:
    if len(group) < 2:
        return
    group.sort()
    for a, b in zip(group, group[1:]):
        cnf.append([-y_lit_all[b], y_lit_all[a]])


def _add_leader_sb(self, cnf, group: List[int], y_lit_all: List[int]) -> None:
    if len(group) < 2:
        return
    group.sort()
    leader = group[0]
    for j in group[1:]:
        cnf.append([-y_lit_all[j], y_lit_all[leader]])


# ============================================================
# 1) COVERAGE + DOMINANCE ONLY
# ============================================================
def sb_coverage_order(
    self,
    cnf,
    radius: float,
    candidates: List[int],
    demands: List[int],
    Nc: Set[int],
) -> Tuple[
    Dict[int, List[int]],     # mask2cands
    Dict[int, Set[int]],      # dom_out
    Dict[int, Set[int]],      # dom_in
]:
    eps = 1e-12
    dist = self.dist
    rad = radius + eps
    ylit = self.y_lit_all

    # ---------- active demands ----------
    if Nc:
        active = []
        Nc_rows = [dist[c] for c in Nc]
        for u in demands:
            for row in Nc_rows:
                if row[u] <= rad:
                    break
            else:
                active.append(u)
    else:
        active = list(demands)

    k = len(active)

    # ---------- masks ----------
    mask2cands: Dict[int, List[int]] = {}
    if k > 0:
        for c in candidates:
            row = dist[c]
            m = 0
            for bit, u in enumerate(active):
                if row[u] <= rad:
                    m |= (1 << bit)
            mask2cands.setdefault(m, []).append(c)
    else:
        mask2cands = {0: list(candidates)}

    for v in mask2cands.values():
        v.sort()

    # ---------- dominance ----------
    dom_out: Dict[int, Set[int]] = {}
    dom_in: Dict[int, Set[int]] = {}

    if k > 0:
        ALL = (1 << k) - 1
        for mi, Hi in mask2cands.items():
            nmi = ALL ^ mi
            for mj, Lj in mask2cands.items():
                if mi.bit_count() <= mj.bit_count():
                    continue
                if (mj & nmi) != 0:
                    continue
                for cj in Lj:
                    for ci in Hi:
                        cnf.append([-ylit[cj], ylit[ci]])
                        dom_out.setdefault(cj, set()).add(ci)
                        dom_in.setdefault(ci, set()).add(cj)

    return mask2cands, dom_out, dom_in


# ============================================================
# 2) ORBIT SYMMETRY BREAKING
# ============================================================
def orbit_symmetry_breaking(
    self,
    cnf,
    *,
    mask2cands: Dict[int, List[int]],
    dom_out: Dict[int, Set[int]],
    dom_in: Dict[int, Set[int]],
    pair_clauses: List[Tuple[int, int]],
    orbit_mode: str = "chain",
) -> None:
    ylit = self.y_lit_all
    n = self.n

    mark = getattr(self, "_sb_mark", [0] * n)
    idx_of = getattr(self, "_sb_idx_of", [-1] * n)
    self._sb_mark = mark
    self._sb_idx_of = idx_of

    tick = getattr(self, "_sb_tick", 0) + 1
    self._sb_tick = tick

    orbit_nodes = []
    for g in mask2cands.values():
        if len(g) >= 2:
            for c in g:
                if mark[c] != tick:
                    mark[c] = tick
                    orbit_nodes.append(c)

    if len(orbit_nodes) < 2:
        return

    for i, c in enumerate(orbit_nodes):
        idx_of[c] = i
    m = len(orbit_nodes)

    pair_mask = [0] * m
    for a, b in pair_clauses:
        ia, ib = idx_of[a], idx_of[b]
        if ia >= 0 and ib >= 0:
            pair_mask[ia] |= (1 << ib)
            pair_mask[ib] |= (1 << ia)

    dom_out_mask = [0] * m
    dom_in_mask = [0] * m
    for c in orbit_nodes:
        ic = idx_of[c]
        for v in dom_out.get(c, ()):
            j = idx_of[v]
            if j >= 0:
                dom_out_mask[ic] |= (1 << j)
        for v in dom_in.get(c, ()):
            j = idx_of[v]
            if j >= 0:
                dom_in_mask[ic] |= (1 << j)

    add_sb = _add_leader_sb if orbit_mode == "leader" else _add_chain_sb

    def key(i):
        return (pair_mask[i], dom_out_mask[i], dom_in_mask[i])

    for group in mask2cands.values():
        sig = {}
        for c in group:
            i = idx_of[c]
            if i >= 0:
                sig.setdefault(key(i), []).append(c)
        for g in sig.values():
            if len(g) >= 2:
                add_sb(self, cnf, g, ylit)

    for c in orbit_nodes:
        idx_of[c] = -1
