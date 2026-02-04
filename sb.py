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
    dom_imps: List[Tuple[int, int]],  # (src, dst) meaning y[src] -> y[dst]
) -> List[List[int]]:
    """
    Colored undirected graph encoding:
      - cover incidence (facility-demand)
      - pair clauses y(a) OR y(b) via pair nodes
      - forced_true markers
      - dominance implications y(src) -> y(dst) via directed gadget:
          For each implication, create:
            imp core node P
            src-end node S, dst-end node T
          edges: P-S, P-T, S-f(src), T-f(dst)
          color classes separate S and T so direction is preserved.

    Return facility orbits only.
    """
    if pynauty is None:
        return []

    m = len(candidates)
    if m <= 1:
        return []

    cand_set = set(candidates)
    k = len(active_demands)

    # filter / normalize inputs to candidates only
    forced = sorted([f for f in forced_true if f in cand_set])

    pair_filtered: List[Tuple[int, int]] = []
    for a, b in pair_clauses:
        if a in cand_set and b in cand_set:
            if a <= b:
                pair_filtered.append((a, b))
            else:
                pair_filtered.append((b, a))

    dom_filtered: List[Tuple[int, int]] = []
    for src, dst in dom_imps:
        if src in cand_set and dst in cand_set and src != dst:
            dom_filtered.append((src, dst))

    q = len(pair_filtered)
    s = len(forced)
    r = len(dom_filtered)

    # If nothing distinguishes facilities, all remaining facilities are symmetric
    if k == 0 and q == 0 and s == 0 and r == 0:
        return [sorted(candidates)]

    # indices
    f_index = {f: i for i, f in enumerate(candidates)}
    d_index = {u: i for i, u in enumerate(active_demands)}

    # Vertex ranges:
    # facilities: [0, m-1]
    # demands   : [m, m+k-1]
    # pairs     : [m+k, m+k+q-1]
    # forces    : [m+k+q, m+k+q+s-1]
    # dom-core  : [base_dom, base_dom+r-1]
    # dom-src   : [base_dom_src, base_dom_src+r-1]
    # dom-dst   : [base_dom_dst, base_dom_dst+r-1]
    base_pairs = m + k
    base_forces = base_pairs + q
    base_dom = base_forces + s
    base_dom_src = base_dom + r
    base_dom_dst = base_dom_src + r
    V = base_dom_dst + r

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

    # pair nodes
    for t, (a, b) in enumerate(pair_filtered):
        va = f_index.get(a)
        vb = f_index.get(b)
        if va is None or vb is None:
            continue
        vp = base_pairs + t
        adj[vp].add(va); adj[va].add(vp)
        adj[vp].add(vb); adj[vb].add(vp)

    # forced true markers
    for t, f in enumerate(forced):
        vf = f_index.get(f)
        if vf is None:
            continue
        vm = base_forces + t
        adj[vm].add(vf); adj[vf].add(vm)

    # dominance implication gadgets (directed in undirected colored graph)
    # For each (src,dst):
    #   core = P, srcEnd = S, dstEnd = T
    #   edges: P-S, P-T, S-f(src), T-f(dst)
    for t, (src, dst) in enumerate(dom_filtered):
        vsrc = f_index.get(src)
        vdst = f_index.get(dst)
        if vsrc is None or vdst is None:
            continue
        P = base_dom + t
        S = base_dom_src + t
        T = base_dom_dst + t

        adj[P].add(S); adj[S].add(P)
        adj[P].add(T); adj[T].add(P)

        adj[S].add(vsrc); adj[vsrc].add(S)
        adj[T].add(vdst); adj[vdst].add(T)

    adjacency_dict = {v: sorted(list(nbs)) for v, nbs in adj.items()}

    # color classes
    facilities = set(range(m))
    demands = set(range(m, m + k))
    pairs = set(range(base_pairs, base_pairs + q))
    forces = set(range(base_forces, base_forces + s))
    dom_core = set(range(base_dom, base_dom + r))
    dom_src = set(range(base_dom_src, base_dom_src + r))
    dom_dst = set(range(base_dom_dst, base_dom_dst + r))

    vertex_coloring = [facilities, demands, pairs, forces, dom_core, dom_src, dom_dst]
    vertex_coloring = [cls for cls in vertex_coloring if len(cls) > 0]

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

        # group facility vertices by orbit id
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
    for orb in facility_orbits:
        if len(orb) < 2:
            continue
        orb = sorted(orb)
        for a, b in zip(orb, orb[1:]):
            cnf.append([-self._y(b), self._y(a)])


def _add_orbit_sb_leader(self, cnf, facility_orbits: List[List[int]]) -> None:
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
    Dominance-based ordering + orbit-based SB computed on a structure that ALSO includes
    the dominance implications, to avoid UNSAT from SB-composition.

    This is the key fix for "dominance + orbit-SB" soundness.
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

    # ----- dominance-based ordering + record implication edges -----
    dom_imps: List[Tuple[int, int]] = []

    for i, ci in enumerate(candidates):
        for cj in candidates[i + 1:]:
            if cover[ci] >= cover[cj] and cover[ci] != cover[cj]:
                # y[cj] -> y[ci]
                cnf.append([-self._y(cj), self._y(ci)])
                dom_imps.append((cj, ci))
            elif cover[cj] >= cover[ci] and cover[cj] != cover[ci]:
                # y[ci] -> y[cj]
                cnf.append([-self._y(ci), self._y(cj)])
                dom_imps.append((ci, cj))

    # ----- orbit-SB on FULL structure (cover + pairs + forced_true + dominance implications) -----
    if enable_orbit_sb and pynauty is not None:
        facility_orbits = _facility_orbits_from_structure(
            cover=cover,
            candidates=candidates,
            active_demands=active_demands,
            pair_clauses=pair_clauses,
            forced_true=forced_true,
            dom_imps=dom_imps,
        )
        if facility_orbits:
            if orbit_mode == "leader":
                _add_orbit_sb_leader(self, cnf, facility_orbits)
            else:
                _add_orbit_sb_chain(self, cnf, facility_orbits)
