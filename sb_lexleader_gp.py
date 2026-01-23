from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict

try:
    import pynauty  # pip install pynauty
except Exception:
    pynauty = None


def add_lex_leq(cnf, A: List[int], B: List[int], top_id: int) -> int:
    """
    Enforce boolean vector A <=lex B.
    A,B are lists of positive var ids (same length).
    Returns new top_id.
    """
    assert len(A) == len(B)

    def new_var():
        nonlocal top_id
        top_id += 1
        return top_id

    # e = True (prefix-equal so far)
    e = new_var()
    cnf.append([e])

    for a, b in zip(A, B):
        # if prefix equal then a <= b : (¬e ∨ ¬a ∨ b)
        cnf.append([-e, -a, b])

        e2 = new_var()

        # e2 -> e
        cnf.append([-e2, e])

        # e2 -> (a <-> b)
        cnf.append([-e2, -a, b])
        cnf.append([-e2, a, -b])

        # (e ∧ (a <-> b)) -> e2
        cnf.append([-e, -a, -b, e2])
        cnf.append([-e, a, b, e2])

        e = e2

    return top_id


def add_lexleader_for_perm(cnf, y_ordered: List[int], perm_img: List[int], top_id: int) -> int:
    """
    perm_img is image mapping on indices:
      i -> perm_img[i]
    Add y <=lex y∘perm  i.e. y[i] <=lex y[perm_img[i]].
    """
    B = [y_ordered[perm_img[i]] for i in range(len(y_ordered))]
    return add_lex_leq(cnf, y_ordered, B, top_id)


@dataclass
class AutoSBResult:
    perms: List[List[int]]      # list of generator perms on candidates (index mapping)
    orbits: List[List[int]]     # orbits of candidates indices (optional)


class SymmetryBreaker:
    """
    Compute automorphism generators of incidence bipartite graph and add LexLeader constraints.
    """
    def __init__(self):
        self.available = pynauty is not None

    def _build_bipartite_graph(
        self,
        cover: Dict[int, List[int]],
        uncovered_demands: List[int],
        candidates: List[int],
    ):
        """
        Graph vertices:
          0..k-1: demands
          k..k+m-1: candidates
        Edges demand_i -- cand_j if candidate covers that demand.
        Vertex coloring forces bipartite partition.
        """
        if pynauty is None:
            raise RuntimeError("pynauty not available")

        k = len(uncovered_demands)
        m = len(candidates)
        idx_d = {u: i for i, u in enumerate(uncovered_demands)}
        idx_c = {c: j for j, c in enumerate(candidates)}

        # adjacency dict: v -> list of neighbours
        adj = {i: [] for i in range(k + m)}

        for u in uncovered_demands:
            du = idx_d[u]
            for c in cover[u]:
                cj = idx_c[c]
                v_c = k + cj
                adj[du].append(v_c)
                adj[v_c].append(du)

        # vertex_coloring expects list of sets (cells)
        coloring = [set(range(0, k)), set(range(k, k + m))]
        g = pynauty.Graph(number_of_vertices=k + m, adjacency_dict=adj, vertex_coloring=coloring)
        return g, k, m

    def automorphism_generators(
        self,
        cover: Dict[int, List[int]],
        uncovered_demands: List[int],
        candidates: List[int],
        *,
        max_gens: int = 64,
    ) -> AutoSBResult:
        """
        Return generator permutations restricted to candidates.
        """
        if pynauty is None:
            return AutoSBResult(perms=[], orbits=[])

        g, k, m = self._build_bipartite_graph(cover, uncovered_demands, candidates)

        gens, size, orbits = pynauty.autgrp(g)  # gens: list of perms of length (k+m)
        cand_perms: List[List[int]] = []

        for perm in gens:
            # restrict to candidate part
            img = [0] * m
            moved = False
            for i in range(m):
                v = k + i
                v_img = perm[v]
                j = v_img - k
                img[i] = j
                if j != i:
                    moved = True
            if moved:
                cand_perms.append(img)
            if len(cand_perms) >= max_gens:
                break

        # Restrict orbits to candidates too (optional)
        cand_orbits: List[List[int]] = []
        for orb in orbits:
            xs = [v - k for v in orb if v >= k]
            if xs:
                cand_orbits.append(sorted(xs))

        return AutoSBResult(perms=cand_perms, orbits=cand_orbits)

    def add_sb(
        self,
        cnf,
        *,
        y_vars_ordered: List[int],         # y(candidates[i]) in fixed order
        cover: Dict[int, List[int]],
        uncovered_demands: List[int],
        candidates: List[int],
        top_id: int,
        max_gens: int = 64,
    ) -> int:
        """
        Add LexLeader for automorphism generators.
        """
        res = self.automorphism_generators(cover, uncovered_demands, candidates, max_gens=max_gens)
        for img in res.perms:
            top_id = add_lexleader_for_perm(cnf, y_vars_ordered, img, top_id)
        return top_id
