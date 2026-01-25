# symmetry_breaking_oec.py
# Orbit / Equivalence-Class Symmetry Breaking for PCenterSAT CNF
#
# Idea:
#   Centers are interchangeable (can be permuted) iff they appear identically in CNF.
#   In your encoding, a center c appears via:
#     - coverage clauses: (y(c1) v y(c2) v ...) for each uncovered demand u
#     - optional extra_pairs clauses: (y(a) v y(b))
#     - optional unit_true: (y(c)) fixing y(c)=True
#
# We build an equivalence signature for each candidate c:
#   signature(c) = (CoveredDemandsSet(c), ExtraPartnersSet(c), UnitFlag(c))
#
# Then for each class [c1 < c2 < ... < ck], add chain implications:
#   y(c2) -> y(c1), y(c3) -> y(c2), ..., y(ck) -> y(c_{k-1})
# CNF clauses: (-y(ci) v y(ci-1))
#
# This is a cheap and usually very effective SB for large TSPLIB instances.

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple, Any
from collections import defaultdict


@dataclass
class OECStats:
    enabled: bool
    n_candidates_in: int
    n_candidates_used: int
    n_classes: int
    n_classes_ge2: int
    max_class_size: int
    sb_clauses_added: int


class OrbitEquivSB:
    """
    Orbit / Equivalence-class SB.

    Usage:
        sb = OrbitEquivSB(debug=True)
        stats = sb.add_sb(
            cnf,
            candidates=candidates,
            cover=cover,  # dict: demand -> list of centers covering demand
            y_func=self._y,
            extra_pairs=sb_extra_pairs,  # optional list of (a,b) center indices
            unit_true=set(sb_units_true),  # optional set of center indices fixed true
        )
    """

    def __init__(self, *, debug: bool = False):
        self.debug = bool(debug)

    @staticmethod
    def _build_center_to_demands(
        candidates: Sequence[int],
        cover: Dict[int, Sequence[int]],
    ) -> Dict[int, List[int]]:
        """
        Build mapping center -> sorted list of uncovered demands it can cover.
        cover[u] contains centers (indices) that can cover u.
        """
        cand_set = set(candidates)
        c2d: Dict[int, List[int]] = {c: [] for c in candidates}

        # Iterate demands in deterministic order for stable signatures
        # (In your code, uncovered_demands is appended in increasing u already,
        # but here cover is dict; we sort keys for safety.)
        for u in sorted(cover.keys()):
            lst = cover.get(u, [])
            for c in lst:
                if c in cand_set:
                    c2d[c].append(u)

        # c2d[c] is already increasing because u loop is increasing
        return c2d

    @staticmethod
    def _build_extra_partner_lists(
        candidates: Sequence[int],
        extra_pairs: Optional[Sequence[Tuple[int, int]]],
    ) -> Dict[int, List[int]]:
        """
        Build mapping center -> sorted list of partners it is linked with by extra_pairs.
        extra_pairs represents clauses (y(a) v y(b)).
        """
        partners: Dict[int, List[int]] = {c: [] for c in candidates}
        if not extra_pairs:
            return partners

        cand_set = set(candidates)

        for a, b in extra_pairs:
            if a in cand_set and b in cand_set:
                partners[a].append(b)
                partners[b].append(a)

        for c in partners:
            partners[c].sort()
        return partners

    @staticmethod
    def _signature(
        c: int,
        c2d: Dict[int, List[int]],
        c2p: Dict[int, List[int]],
        unit_true: Set[int],
    ) -> Tuple[Tuple[int, ...], Tuple[int, ...], int]:
        """
        Exact (collision-free) signature:
          - tuple of covered demand IDs
          - tuple of extra partners
          - unit flag (1 if fixed-true else 0)
        """
        return (tuple(c2d.get(c, [])), tuple(c2p.get(c, [])), 1 if c in unit_true else 0)

    def add_sb(
        self,
        cnf: Any,
        *,
        candidates: Sequence[int],
        cover: Dict[int, Sequence[int]],
        y_func: Callable[[int], int],
        extra_pairs: Optional[Sequence[Tuple[int, int]]] = None,
        unit_true: Optional[Iterable[int]] = None,
        min_class_size: int = 2,
        exclude_unit_true: bool = True,
    ) -> OECStats:
        """
        Add OEC SB constraints into cnf.

        Parameters
        ----------
        cnf:
            pysat.formula.CNF or any object with .append(clause) and .clauses (optional)
        candidates:
            list of remaining candidate centers (indices in [0..n-1]) whose y-vars exist.
        cover:
            dict demand->list of candidate centers covering it (only uncovered demands).
        y_func:
            function mapping center index -> CNF var id (positive int).
        extra_pairs:
            optional list of pairs (a,b) where clause (y(a) v y(b)) exists.
            IMPORTANT: pass the filtered list only containing pairs among current candidates
                       (you already compute sb_extra_pairs for that).
        unit_true:
            optional iterable of centers fixed to True by unit clause (y(c)).
            IMPORTANT: pass the filtered list only among current candidates
                       (you already compute sb_units_true for that).
        min_class_size:
            only add SB for equivalence classes with size >= this.
        exclude_unit_true:
            if True, we do not apply SB constraints that involve centers fixed true.
            (Recommended: unit-true variables do not contribute to branching; SB there is redundant.)

        Returns
        -------
        OECStats describing what was added.
        """
        cand_in = list(candidates)
        n_in = len(cand_in)

        if n_in == 0:
            # Nothing meaningful to break
            return OECStats(
                enabled=False,
                n_candidates_in=n_in,
                n_candidates_used=0,
                n_classes=0,
                n_classes_ge2=0,
                max_class_size=0,
                sb_clauses_added=0,
            )

        unit_true_set: Set[int] = set(unit_true) if unit_true is not None else set()

        # Optionally remove fixed-true candidates from SB universe
        if exclude_unit_true and unit_true_set:
            cand = [c for c in cand_in if c not in unit_true_set]
        else:
            cand = cand_in

        if len(cand) <= 1:
            return OECStats(
                enabled=False,
                n_candidates_in=n_in,
                n_candidates_used=len(cand),
                n_classes=len(cand),
                n_classes_ge2=0,
                max_class_size=(len(cand) if cand else 0),
                sb_clauses_added=0,
            )

        # Build exact signatures
        c2d = self._build_center_to_demands(cand, cover)
        c2p = self._build_extra_partner_lists(cand, extra_pairs)

        classes: Dict[Tuple[Tuple[int, ...], Tuple[int, ...], int], List[int]] = defaultdict(list)
        for c in cand:
            sig = self._signature(c, c2d, c2p, unit_true_set)
            classes[sig].append(c)

        # Add chain implications per class
        sb_added = 0
        max_sz = 1
        n_ge2 = 0

        for sig, group in classes.items():
            if len(group) > max_sz:
                max_sz = len(group)

            if len(group) < max(min_class_size, 2):
                continue

            group.sort()
            n_ge2 += 1

            # y(group[i]) -> y(group[i-1])  for i=1..k-1
            for i in range(1, len(group)):
                ci = group[i]
                cj = group[i - 1]
                cnf.append([-int(y_func(ci)), int(y_func(cj))])
                sb_added += 1

        if self.debug:
            # Lightweight debug printing (avoid huge dumps)
            print(
                f"[OEC-SB] candidates_in={n_in} used={len(cand)} "
                f"classes={len(classes)} classes_ge{min_class_size}={n_ge2} "
                f"max_class={max_sz} sb_clauses={sb_added}",
                flush=True,
            )

        return OECStats(
            enabled=(sb_added > 0),
            n_candidates_in=n_in,
            n_candidates_used=len(cand),
            n_classes=len(classes),
            n_classes_ge2=n_ge2,
            max_class_size=max_sz,
            sb_clauses_added=sb_added,
        )
