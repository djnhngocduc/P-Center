# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional


def _add_chain_sb(self, cnf, group: List[int], y_lit_all: List[int]) -> None:
    if len(group) < 2:
        return
    group.sort()
    append = cnf.append
    ylit = y_lit_all
    for a, b in zip(group, group[1:]):
        append([-ylit[b], ylit[a]])


def _add_leader_sb(self, cnf, group: List[int], y_lit_all: List[int]) -> None:
    if len(group) < 2:
        return
    group.sort()
    leader = group[0]
    yl = y_lit_all[leader]
    append = cnf.append
    ylit = y_lit_all
    for j in group[1:]:
        append([-ylit[j], yl])


def sb_coverage_order(
    self,
    cnf,
    radius: float,
    candidates: List[int],
    demands: List[int],
    Nc: Set[int],
    *,
    enable_orbit_sb: bool = True,
    orbit_mode: str = "chain",
    pair_clauses: Optional[List[Tuple[int, int]]] = None,
) -> None:
    eps = 1e-12
    if pair_clauses is None:
        pair_clauses = []

    if len(candidates) <= 1:
        return

    dist = self.dist
    rad = radius + eps
    n = self.n

    y_lit_all = getattr(self, "y_lit_all", None)
    if y_lit_all is None:
        y = self._y
        y_lit_all = [y(i) for i in range(n)]

    append = cnf.append
    ylit = y_lit_all

    # -----------------------
    # Active demands
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
    # A1: cover bitmasks + group
    # -----------------------
    mask2cands: Dict[int, List[int]] = {}

    if k > 0:
        act = active_demands
        for c in candidates:
            row = dist[c]
            m = 0
            for bit, u in enumerate(act):
                if row[u] <= rad:
                    m |= (1 << bit)
            mask2cands.setdefault(m, []).append(c)
    else:
        mask2cands = {0: list(candidates)}

    for lst in mask2cands.values():
        lst.sort()

    # -----------------------
    # A2: bucket by popcount
    # -----------------------
    size2masks: Dict[int, List[int]] = {}
    for m in mask2cands.keys():
        size2masks.setdefault(m.bit_count(), []).append(m)

    sizes = sorted(size2masks.keys(), reverse=True)
    for s in sizes:
        size2masks[s].sort()

    # -----------------------
    # Dominance: LAZY
    # -----------------------
    dom_out: Dict[int, Set[int]] = {}
    dom_in: Dict[int, Set[int]] = {}

    # -----------------------
    # FULL dominance ordering
    # -----------------------
    if k > 0 and len(sizes) > 1:
        ALL = (1 << k) - 1
        dom_out_setdefault = dom_out.setdefault
        dom_in_setdefault = dom_in.setdefault

        for i, s_hi in enumerate(sizes):
            hi_masks = size2masks[s_hi]
            for s_lo in sizes[i + 1:]:
                lo_masks = size2masks[s_lo]

                for mi in hi_masks:
                    hi_nodes = mask2cands[mi]
                    nmi = ALL ^ mi

                    for mj in lo_masks:
                        if (mj & nmi) != 0:
                            continue
                        lo_nodes = mask2cands[mj]

                        for cj in lo_nodes:
                            out_cj = dom_out_setdefault(cj, set())
                            ycj = -ylit[cj]
                            for ci in hi_nodes:
                                if ci in out_cj:
                                    continue
                                append([ycj, ylit[ci]])
                                out_cj.add(ci)
                                dom_in_setdefault(ci, set()).add(cj)

    # -----------------------
    # orbit-SB optimized
    #   - no O(n) fresh allocations: reuse buffers on self
    #   - signature key is INT (packed), no tuple alloc
    # -----------------------
    if not enable_orbit_sb:
        return

    # reuse buffers
    mark = getattr(self, "_sb_mark", None)
    idx_of = getattr(self, "_sb_idx_of", None)
    if mark is None or len(mark) != n:
        mark = [0] * n
        self._sb_mark = mark
    if idx_of is None or len(idx_of) != n:
        idx_of = [-1] * n
        self._sb_idx_of = idx_of

    tick = getattr(self, "_sb_tick", 0) + 1
    if tick == 0x7FFFFFFF:
        # rare: reset marks
        for i in range(n):
            mark[i] = 0
        tick = 1
    self._sb_tick = tick

    orbit_nodes: List[int] = []
    for group in mask2cands.values():
        if len(group) >= 2:
            for c in group:
                if mark[c] != tick:
                    mark[c] = tick
                    orbit_nodes.append(c)

    if not orbit_nodes:
        return

    # assign orbit indices
    for i, c in enumerate(orbit_nodes):
        idx_of[c] = i
    m_orb = len(orbit_nodes)

    pair_mask_orb = [0] * m_orb
    if pair_clauses:
        for a, b in pair_clauses:
            ia = idx_of[a]
            if ia < 0:
                continue
            ib = idx_of[b]
            if ib < 0:
                continue
            pair_mask_orb[ia] |= (1 << ib)
            pair_mask_orb[ib] |= (1 << ia)

    dom_out_mask_orb = [0] * m_orb
    dom_in_mask_orb = [0] * m_orb

    dom_out_get = dom_out.get
    dom_in_get = dom_in.get

    for c in orbit_nodes:
        ic = idx_of[c]

        mo = 0
        for v in dom_out_get(c, ()):
            j = idx_of[v]
            if j >= 0:
                mo |= (1 << j)
        dom_out_mask_orb[ic] = mo

        mi = 0
        for v in dom_in_get(c, ()):
            j = idx_of[v]
            if j >= 0:
                mi |= (1 << j)
        dom_in_mask_orb[ic] = mi

    add_sb = _add_leader_sb if orbit_mode == "leader" else _add_chain_sb

    # pack 3 masks into 1 int key to avoid tuple alloc/hash
    # (Python int is arbitrary precision, safe)
    def pack_key(pm: int, om: int, im: int) -> int:
        return (pm << (2 * m_orb)) | (om << m_orb) | im

    for group in mask2cands.values():
        if len(group) < 2:
            continue

        sig2group: Dict[int, List[int]] = {}
        count = 0

        for c in group:
            ic = idx_of[c]
            if ic < 0:
                continue
            count += 1
            key = pack_key(pair_mask_orb[ic], dom_out_mask_orb[ic], dom_in_mask_orb[ic])
            sig2group.setdefault(key, []).append(c)

        if count < 2:
            continue

        for g in sig2group.values():
            if len(g) >= 2:
                add_sb(self, cnf, g, ylit)

    # cleanup idx_of for orbit nodes only (avoid O(n))
    for c in orbit_nodes:
        idx_of[c] = -1
