from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict, Optional
import json

try:
    import pynauty  # pip install pynauty
except Exception:
    pynauty = None


def add_lex_leq(cnf, A: List[int], B: List[int], top_id: int) -> int:
    """
    Enforce boolean vector A <=lex B.
    A,B are lists of positive var ids (same length).
    Returns updated top_id.
    """
    assert len(A) == len(B)

    def new_var() -> int:
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
        # (¬e ∨ ¬a ∨ b ∨ e2)  and  (¬e ∨ a ∨ ¬b ∨ e2)
        cnf.append([-e, -a, b, e2])
        cnf.append([-e, a, -b, e2])

        e = e2

    return top_id


def add_lexleader_for_perm(cnf, y_ordered: List[int], perm_img: List[int], top_id: int) -> int:
    """
    perm_img: image mapping on indices 0..m-1:
      i -> perm_img[i]
    Add y <=lex y∘perm  i.e. y[i] <=lex y[perm_img[i]].
    """
    m = len(y_ordered)
    if len(perm_img) != m:
        raise ValueError(f"perm_img length {len(perm_img)} != {m}")

    B = [y_ordered[perm_img[i]] for i in range(m)]
    return add_lex_leq(cnf, y_ordered, B, top_id)


def _validate_perm_img(img: List[int], m: int) -> None:
    if len(img) != m:
        raise ValueError(f"Permutation length {len(img)} != {m}")
    s = set(img)
    if s != set(range(m)):
        raise ValueError("perm_img is not a permutation of 0..m-1")


@dataclass
class AutoSBResult:
    perms: List[List[int]]      # generator perms on candidates (index mapping)
    orbits: List[List[int]]     # orbits of candidates indices (optional)


class SymmetryBreaker:
    """
    Two modes:
      - canonizing: use a provided canonizing set Π (list of perms) and add LexLeader for each π∈Π.
      - generators: compute automorphism generators and add LexLeader for each generator (partial SB).
    """
    def __init__(
        self,
        *,
        canonizing_perms: Optional[List[List[int]]] = None,
        canonizing_perms_path: Optional[str] = None,
        use_generators_fallback: bool = True,
    ):
        self.available = (pynauty is not None)
        self.use_generators_fallback = use_generators_fallback

        self._canonizing_perms: Optional[List[List[int]]] = None
        if canonizing_perms_path is not None:
            self._canonizing_perms = self._load_canonizing_perms(canonizing_perms_path)
        elif canonizing_perms is not None:
            self._canonizing_perms = canonizing_perms

    @staticmethod
    def _load_canonizing_perms(path: str) -> List[List[int]]:
        """
        Expected JSON formats:
          1) {"perms": [[...], [...], ...]}
          2) [[...], [...], ...]
        """
        with open(path, "r", encoding="utf-8") as f:
            obj = json.load(f)
        if isinstance(obj, dict):
            perms = obj.get("perms", None)
            if perms is None:
                raise ValueError("JSON dict must contain key 'perms'")
            return perms
        if isinstance(obj, list):
            return obj
        raise ValueError("Unsupported JSON format for canonizing perms")

    def set_canonizing_perms(self, perms: Optional[List[List[int]]]) -> None:
        self._canonizing_perms = perms

    # ----- generator-based (fallback) -----
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

        adj = {i: [] for i in range(k + m)}

        for u in uncovered_demands:
            du = idx_d[u]
            for c in cover[u]:
                cj = idx_c[c]
                v_c = k + cj
                adj[du].append(v_c)
                adj[v_c].append(du)

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
        if pynauty is None:
            return AutoSBResult(perms=[], orbits=[])

        g, k, m = self._build_bipartite_graph(cover, uncovered_demands, candidates)

        gens, size, orbits = pynauty.autgrp(g)  # perms of length (k+m)
        cand_perms: List[List[int]] = []

        for perm in gens:
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

        cand_orbits: List[List[int]] = []
        for orb in orbits:
            xs = [v - k for v in orb if v >= k]
            if xs:
                cand_orbits.append(sorted(xs))

        return AutoSBResult(perms=cand_perms, orbits=cand_orbits)

    # ----- main entry -----
    def add_sb(
        self,
        cnf,
        *,
        y_vars_ordered: List[int],
        cover: Dict[int, List[int]],
        uncovered_demands: List[int],
        candidates: List[int],
        top_id: int,
        max_gens: int = 64,
    ) -> int:
        """
        If canonizing perms Π is provided => add LexLeader for each π∈Π (this matches the paper's 'canonizing set' notion).
        Else (optional) fallback to generator perms.
        """
        m = len(y_vars_ordered)

        # 1) Prefer canonizing set Π if provided
        if self._canonizing_perms is not None:
            for img in self._canonizing_perms:
                _validate_perm_img(img, m)
                top_id = add_lexleader_for_perm(cnf, y_vars_ordered, img, top_id)
            return top_id

        # 2) Otherwise fallback to generators (partial SB)
        if not self.use_generators_fallback:
            return top_id

        if not self.available:
            return top_id

        res = self.automorphism_generators(cover, uncovered_demands, candidates, max_gens=max_gens)
        for img in res.perms:
            _validate_perm_img(img, m)
            top_id = add_lexleader_for_perm(cnf, y_vars_ordered, img, top_id)
        return top_id
