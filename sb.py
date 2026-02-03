# sb.py
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Optional

try:
    import pynauty
except Exception:
    pynauty = None


def _facility_orbits_from_incidence_graph(
    cover: Dict[int, Set[int]],
    candidates: List[int],
    active_demands: List[int],
) -> List[List[int]]:
    """
    Build a bipartite incidence graph:
      left vertices  : facilities (candidates)
      right vertices : active_demands
      edge (f, d) iff d in cover[f]

    Return orbits for facility vertices only, in terms of original facility ids (candidates values).

    If pynauty is unavailable or something fails, return [].
    """
    if pynauty is None:
        return []

    m = len(candidates)
    k = len(active_demands)
    if m <= 1:
        return []

    # Map demand id -> right-vertex index 0..k-1
    d_index = {u: idx for idx, u in enumerate(active_demands)}

    # Build adjacency for graph with V = m + k vertices (0..m-1 facilities, m..m+k-1 demands)
    V = m + k
    adj: Dict[int, Set[int]] = {v: set() for v in range(V)}

    # Add edges
    for fi, f in enumerate(candidates):
        # fi in [0..m-1]
        for u in cover.get(f, ()):
            dj = d_index.get(u)
            if dj is None:
                continue
            v_d = m + dj
            adj[fi].add(v_d)
            adj[v_d].add(fi)

    # Convert adjacency sets to sorted lists as pynauty expects
    adjacency_dict = {v: sorted(list(nbs)) for v, nbs in adj.items()}

    # Vertex coloring (partitions) to forbid swapping facility<->demand
    # Two color classes: facilities and demands
    facility_vertices = set(range(m))
    demand_vertices = set(range(m, m + k))
    vertex_coloring = [facility_vertices, demand_vertices]

    try:
        g = pynauty.Graph(
            number_of_vertices=V,
            directed=False,
            adjacency_dict=adjacency_dict,
            vertex_coloring=vertex_coloring,
        )

        # Different pynauty builds expose different APIs.
        # We'll try a few common ones safely.
        orbits = None

        # 1) pynauty.autgrp(g) -> (gens, grpsize1, grpsize2, orbits, numorbits) in some versions
        if hasattr(pynauty, "autgrp"):
            res = pynauty.autgrp(g)
            # try to locate 'orbits' inside res
            # res might be tuple; orbits usually is a list of ints length V
            if isinstance(res, tuple):
                for item in res:
                    if isinstance(item, list) and len(item) == V and all(isinstance(x, int) for x in item):
                        orbits = item
                        break

        # 2) pynauty.orbits(g) -> list[int] in some versions
        if orbits is None and hasattr(pynauty, "orbits"):
            tmp = pynauty.orbits(g)
            if isinstance(tmp, list) and len(tmp) == V:
                orbits = tmp

        # If still unknown, bail out safely
        if orbits is None:
            return []

        # Group facility vertices by orbit id
        groups: Dict[int, List[int]] = {}
        for fi in range(m):
            oid = orbits[fi]
            groups.setdefault(oid, []).append(fi)

        # Keep only non-trivial orbits (size >= 2), map back to facility ids
        facility_orbits: List[List[int]] = []
        for _, idxs in groups.items():
            if len(idxs) >= 2:
                orb_facilities = [candidates[fi] for fi in sorted(idxs)]
                facility_orbits.append(orb_facilities)

        # Sort orbits for determinism
        facility_orbits.sort(key=lambda lst: (len(lst), lst[0]))
        return facility_orbits

    except Exception:
        # Must be safe: if orbit extraction fails, just skip.
        return []


def _add_orbit_sb_chain(self, cnf, facility_orbits: List[List[int]]) -> None:
    """
    Orbit-based symmetry breaking (cheap):
      for each orbit [j1<j2<...<jk], enforce y[j1] >= y[j2] >= ... >= y[jk]
    CNF: (-y[j_{t+1}] OR y[j_t])
    """
    for orb in facility_orbits:
        if len(orb) < 2:
            continue
        orb = sorted(orb)
        for a, b in zip(orb, orb[1:]):
            # y[b] -> y[a]
            cnf.append([-self._y(b), self._y(a)])


def _add_orbit_sb_leader(self, cnf, facility_orbits: List[List[int]]) -> None:
    """
    Orbit-based SB (leader form):
      for orbit orb, leader = min(orb), enforce y[j] -> y[leader] for all j != leader.
    CNF: (-y[j] OR y[leader])
    """
    for orb in facility_orbits:
        if len(orb) < 2:
            continue
        leader = min(orb)
        for j in orb:
            if j == leader:
                continue
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
    orbit_mode: str = "chain",   # "chain" or "leader"
):
    """
    STRONG & SAFE dominance-based ordering:
    enforce ordering only when coverage inclusion holds.

    Optional: orbit-based SB via automorphisms of the incidence graph (facility-demands),
    applied on remaining candidates and active demands.
    """

    eps = 1e-12

    # active demands = chưa được cover bởi Nc
    active_demands = [
        u for u in demands
        if not any(self.dist[c][u] <= radius + eps for c in Nc)
    ]

    # build coverage sets
    cover: Dict[int, Set[int]] = {}
    for c in candidates:
        cover[c] = {u for u in active_demands if self.dist[c][u] <= radius + eps}

    # -------- dominance-based ordering (your current SB) --------
    for i, ci in enumerate(candidates):
        for cj in candidates[i + 1:]:
            if cover[ci] >= cover[cj] and cover[ci] != cover[cj]:
                # cj -> ci
                cnf.append([-self._y(cj), self._y(ci)])
            elif cover[cj] >= cover[ci] and cover[cj] != cover[ci]:
                # ci -> cj
                cnf.append([-self._y(ci), self._y(cj)])

    # -------- orbit-based SB (optional) --------
    if enable_orbit_sb:
        facility_orbits = _facility_orbits_from_incidence_graph(
            cover=cover,
            candidates=candidates,
            active_demands=active_demands,
        )
        if facility_orbits:
            if orbit_mode == "leader":
                _add_orbit_sb_leader(self, cnf, facility_orbits)
            else:
                # default: chain
                _add_orbit_sb_chain(self, cnf, facility_orbits)
