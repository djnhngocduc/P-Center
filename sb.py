# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional

try:
    import pynauty
except Exception:
    pynauty = None


def _facility_orbits_from_structure(
    cover: Dict[int, Set[int]],
    candidates: List[int],
    active_demands: List[int],
    pair_clauses: List[Tuple[int, int]],
    forced_true: Set[int],
) -> List[List[int]]:
    """
    Build a colored graph capturing constraints that distinguish facilities:

      Part 1: facilities (candidates)
      Part 2: active demands
      Part 3: pair-constraint nodes, one for each clause (y(a) OR y(b)) where both a,b are in candidates
      Part 4: forced-true markers, one node per forced facility in candidates

    Edges:
      - facility -- demand iff demand in cover[facility]
      - pairNode -- facility(a) and pairNode -- facility(b)
      - forceNode -- facility(f)  (marks y(f) is forced TRUE after simplification)

    Return facility orbits only (list of lists of facility ids).

    If pynauty is unavailable or fails, return [] (safe: orbit-SB disabled).
    """
    if pynauty is None:
        return []

    m = len(candidates)
    if m <= 1:
        return []

    k = len(active_demands)
    # Filter forced_true to only candidates
    cand_set = set(candidates)
    forced = sorted([f for f in forced_true if f in cand_set])
    s = len(forced)

    # Filter pair_clauses to only those with both endpoints in candidates
    pair_filtered: List[Tuple[int, int]] = []
    for a, b in pair_clauses:
        if a in cand_set and b in cand_set:
            pair_filtered.append((a, b))
    q = len(pair_filtered)

    # If nothing distinguishes facilities, all remaining facilities are symmetric
    if k == 0 and q == 0 and s == 0:
        return [sorted(candidates)]

    # Map facility id -> facility-vertex index 0..m-1
    f_index = {f: i for i, f in enumerate(candidates)}
    # Map demand id -> demand-vertex index 0..k-1
    d_index = {u: i for i, u in enumerate(active_demands)}

    # Vertex ranges:
    # facilities: [0, m-1]
    # demands   : [m, m+k-1]
    # pairs     : [m+k, m+k+q-1]
    # forces    : [m+k+q, m+k+q+s-1]
    V = m + k + q + s
    adj: Dict[int, Set[int]] = {v: set() for v in range(V)}

    # facility-demand edges
    for fi, f in enumerate(candidates):
        for u in cover.get(f, ()):
            dj = d_index.get(u)
            if dj is None:
                continue
            v_d = m + dj
            adj[fi].add(v_d)
            adj[v_d].add(fi)

    # pair constraint nodes (only between candidate facilities)
    base_pair = m + k
    for t, (a, b) in enumerate(pair_filtered):
        va = f_index.get(a)
        vb = f_index.get(b)
        if va is None or vb is None:
            continue
        vp = base_pair + t
        adj[vp].add(va); adj[va].add(vp)
        adj[vp].add(vb); adj[vb].add(vp)

    # forced-true marker nodes
    base_force = m + k + q
    for t, f in enumerate(forced):
        vf = f_index.get(f)
        if vf is None:
            continue
        vmark = base_force + t
        adj[vmark].add(vf); adj[vf].add(vmark)

    adjacency_dict = {v: sorted(list(nbs)) for v, nbs in adj.items()}

    # Color classes to forbid swapping between parts
    facilities = set(range(m))
    demands = set(range(m, m + k))
    pairs = set(range(m + k, m + k + q))
    forces = set(range(m + k + q, V))
    vertex_coloring = [facilities, demands, pairs, forces]

    # remove empty color classes (some pynauty builds dislike empty partitions)
    vertex_coloring = [cls for cls in vertex_coloring if len(cls) > 0]

    try:
        g = pynauty.Graph(
            number_of_vertices=V,
            directed=False,
            adjacency_dict=adjacency_dict,
            vertex_coloring=vertex_coloring,
        )

        orbits = None

        # Different pynauty builds expose different APIs; try common ones.
        if hasattr(pynauty, "autgrp"):
            res = pynauty.autgrp(g)
            # common: res[3] is orbits list of length V
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

        # Group facility vertices by orbit id
        groups: Dict[int, List[int]] = {}
        for fi in range(m):
            oid = orbits[fi]
            groups.setdefault(oid, []).append(fi)

        facility_orbits: List[List[int]] = []
        for idxs in groups.values():
            if len(idxs) >= 2:
                facility_orbits.append([candidates[fi] for fi in sorted(idxs)])

        facility_orbits.sort(key=lambda lst: (len(lst), lst[0]))
        return facility_orbits

    except Exception:
        return []


def _add_orbit_sb_chain(self, cnf, facility_orbits: List[List[int]]) -> None:
    """
    For each orbit [j1<j2<...<jk], enforce y[j1] >= y[j2] >= ... >= y[jk]
    CNF: (-y[j_{t+1}] OR y[j_t])
    """
    for orb in facility_orbits:
        if len(orb) < 2:
            continue
        orb = sorted(orb)
        for a, b in zip(orb, orb[1:]):
            cnf.append([-self._y(b), self._y(a)])


def _add_orbit_sb_leader(self, cnf, facility_orbits: List[List[int]]) -> None:
    """
    Leader form: for orbit orb, leader=min(orb), enforce y[j] -> y[leader] for all j != leader
    CNF: (-y[j] OR y[leader])
    """
    for orb in facility_orbits:
        if len(orb) < 2:
            continue
        leader = min(orb)
        for j in orb:
            if j != leader:
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
    orbit_mode: str = "chain",
    pair_clauses: Optional[List[Tuple[int, int]]] = None,
    forced_true: Optional[Set[int]] = None,
) -> None:
    """
    Dominance-based ordering + (optional) orbit-based SB.

    Sound orbit-SB requires symmetry computed on a structure that includes ALL constraints
    distinguishing facilities. Here we include:
      - cover constraints (facility-demand incidence)
      - pair clauses y(a) OR y(b) (both endpoints in candidates)
      - forced TRUE facilities induced by pair clauses + Nd simplification (force markers)
    """
    eps = 1e-12
    pair_clauses = pair_clauses or []
    forced_true = forced_true or set()

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

    # orbit-based SB
    if enable_orbit_sb and pynauty is not None:
        facility_orbits = _facility_orbits_from_structure(
            cover=cover,
            candidates=candidates,
            active_demands=active_demands,
            pair_clauses=pair_clauses,
            forced_true=forced_true,
        )
        if facility_orbits:
            if orbit_mode == "leader":
                _add_orbit_sb_leader(self, cnf, facility_orbits)
            else:
                _add_orbit_sb_chain(self, cnf, facility_orbits)
