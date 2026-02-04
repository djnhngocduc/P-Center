import math
from typing import List, Tuple, Set
from threading import RLock
from collections import deque
from pysat.formula import CNF
from pysat.card import CardEnc
from pysat.card import EncType as CardEncType

from pysat.pb import PBEnc
from pysat.pb import EncType as PBEncType

from pypblib import pblib
from pypblib.pblib import PBConfig, Pb2cnf

from sb import sb_coverage_order


_PYSAT_CNF_LOCK = RLock()
DEBUG_REDUCTION = False


class PCenterSAT:
    def __init__(self, dist: List[List[float]], p: int):
        n = len(dist)
        self.n: int = n
        self.p: int = int(p)
        self.dist: List[List[float]] = dist
        self.radii: List[float] = self._compute_radii(dist)

    def _y(self, j: int) -> int:
        # 1..n
        return 1 + j

    @staticmethod
    def from_distance_matrix(dist: List[List[float]], p: int) -> "PCenterSAT":
        return PCenterSAT(dist, p)

    @staticmethod
    def from_coordinates(coords: List[Tuple[float, float]], p: int) -> "PCenterSAT":
        n = len(coords)
        D = [[0.0] * n for _ in range(n)]
        for i, (xi, yi) in enumerate(coords):
            for j, (xj, yj) in enumerate(coords):
                dx = xi - xj
                dy = yi - yj
                D[i][j] = math.hypot(dx, dy)
        return PCenterSAT(D, p)

    @staticmethod
    def _compute_radii(D: List[List[float]]) -> List[float]:
        n = len(D)
        vals = set()
        INF = 10 ** 12
        for i in range(n):
            for j in range(i + 1, n):
                dij = D[i][j]
                if dij != 0 and dij < INF:
                    vals.add(dij)
        return sorted(vals, reverse=True)
    
    def _build_neighbours(self, radius: float) -> List[List[int]]:
        n = self.n
        neigh = [[] for _ in range(n)]
        eps = 1e-12
        for i in range(n):
            row = self.dist[i]
            for j in range(n):
                if i == j:
                    continue
                if row[j] <= radius + eps:
                    neigh[i].append(j)
        return neigh

    @staticmethod
    def _dominates(belong_to: List[Set[int]], u: int, targets: List[int]) -> bool:
        s = belong_to[u]
        for t in targets:
            if t == u:
                continue
            if t not in s:
                return False
        return True

    @staticmethod
    def _bfs_within_3(neighbours: List[List[int]], start: int) -> Set[int]:
        q = deque([(start, 0)])
        seen = {start}
        out = set()
        while q:
            u, d = q.popleft()
            if d == 3:
                continue
            for v in neighbours[u]:
                if v in seen:
                    continue
                seen.add(v)
                out.add(v)
                q.append((v, d + 1))
        return out

    @staticmethod
    def _reduction_merge(neighbours: List[List[int]], v: int, w: int) -> List[int]:
        nv = neighbours[v]
        nw = neighbours[w]
        i = j = 0
        out = []

        while i < len(nv) and j < len(nw):
            if nv[i] == w:
                i += 1
                continue
            if nw[j] == v:
                j += 1
                continue

            if nv[i] < nw[j]:
                out.append(nv[i]); i += 1
            elif nv[i] == nw[j]:
                out.append(nv[i]); i += 1; j += 1
            else:
                out.append(nw[j]); j += 1

        while i < len(nv):
            out.append(nv[i]); i += 1
        while j < len(nw):
            out.append(nw[j]); j += 1
        return out

    @staticmethod
    def _rule1(neighbours: List[List[int]], true_set: Set[int], del_set: Set[int]) -> None:
        n = len(neighbours)
        belong = [set(neighbours[v]) for v in range(n)]

        for v in range(n):
            if v in del_set:
                continue

            Nv = neighbours[v]

            N1 = []
            for u in Nv:
                dominated = False
                for t in neighbours[u]:
                    if t != v and (t not in belong[v]):
                        dominated = True
                        break
                if dominated:
                    N1.append(u)
            N1_set = set(N1)

            N2 = []
            for u in Nv:
                if u in N1_set:
                    continue
                for t in neighbours[u]:
                    if t in N1_set:
                        N2.append(u)
                        break
            N2_set = set(N2)

            if len(Nv) == (len(N1) + len(N2)):
                continue
            N3 = [u for u in Nv if (u not in N1_set and u not in N2_set)]
            if N3 and (v not in del_set):
                true_set.add(v)
                for u in N3:
                    del_set.add(u)
                for u in N2:
                    del_set.add(u)

        del_set.difference_update(true_set)

    @staticmethod
    def _rule2(
        neighbours: List[List[int]],
        at_least_pairs: Set[Tuple[int, int]],
        true_set: Set[int],
        del_set: Set[int],
    ) -> None:
        n = len(neighbours)
        belong = [set(neighbours[v]) for v in range(n)]

        for v in range(n):
            if v in del_set:
                continue

            within3 = PCenterSAT._bfs_within_3(neighbours, v)
            for w in sorted(x for x in within3 if x > v):
                if w in del_set:
                    continue

                N_vw = PCenterSAT._reduction_merge(neighbours, v, w)

                N1 = []
                N_vw_set = set(N_vw)
                for u in N_vw:
                    dominated = False
                    for t in neighbours[u]:
                        if t != v and t != w and (t not in N_vw_set):
                            dominated = True
                            break
                    if dominated:
                        N1.append(u)
                N1_set = set(N1)

                N2 = []
                for u in N_vw:
                    if u in N1_set:
                        continue
                    for t in neighbours[u]:
                        if t in N1_set:
                            N2.append(u)
                            break
                N2_set = set(N2)

                if len(N_vw) == (len(N1) + len(N2)):
                    continue

                N3 = [u for u in N_vw if (u not in N1_set and u not in N2_set)]

                if len(N3) <= 1:
                    continue

                is_dominated = False
                for u in N3:
                    if PCenterSAT._dominates(belong, u, N3):
                        is_dominated = True
                        break
                if not is_dominated:
                    for u in N2:
                        if PCenterSAT._dominates(belong, u, N3):
                            is_dominated = True
                            break
                if is_dominated:
                    continue

                dominate_v = PCenterSAT._dominates(belong, v, N3)
                dominate_w = PCenterSAT._dominates(belong, w, N3)

                if dominate_v or dominate_w:
                    if dominate_v and dominate_w:
                        if (v not in true_set) and (w not in true_set):
                            a, b = (v, w) if v < w else (w, v)
                            at_least_pairs.add((a, b))
                        for u in N3:
                            del_set.add(u)
                        for u in N2:
                            if (u in belong[v]) and (u in belong[w]):
                                del_set.add(u)
                    elif dominate_v and (not dominate_w):
                        true_set.add(v)
                        for u in N3:
                            del_set.add(u)
                        for u in N2:
                            if u in belong[v]:
                                del_set.add(u)
                    elif dominate_w and (not dominate_v):
                        true_set.add(w)
                        for u in N3:
                            del_set.add(u)
                        for u in N2:
                            if u in belong[w]:
                                del_set.add(u)
                else:
                    true_set.add(v)
                    true_set.add(w)
                    for u in N3:
                        if u not in true_set:
                            del_set.add(u)
                    for u in N2:
                        if u not in true_set:
                            del_set.add(u)

        del_set.difference_update(true_set)

    @staticmethod
    def _handle_white_nodes(
        neighbours: List[List[int]],
        is_white: List[bool],
        new_neigh: List[List[int]],
        true_set: Set[int],
        del_set: Set[int],
    ) -> None:
        n = len(neighbours)

        for i in range(n):
            if not is_white[i]:
                continue
            new_neigh[i] = [v for v in new_neigh[i] if not is_white[v]]

        for i in range(n):
            if not is_white[i]:
                continue
            if i in true_set or i in del_set:
                continue
            deg = len(new_neigh[i])
            if deg == 0:
                del_set.add(i)
            elif deg == 1:
                node = new_neigh[i][0]
                if i in new_neigh[node]:
                    new_neigh[node].remove(i)
                del_set.add(i)
                new_neigh[i].clear()

    @classmethod
    def _addition_rule(
        cls,
        neighbours: List[List[int]],
        at_least_pairs: Set[Tuple[int, int]],
        true_set: Set[int],
        del_set: Set[int],
    ) -> None:
        n = len(neighbours)
        is_white = [False] * n
        new_neigh = [[] for _ in range(n)]

        for i in range(n):
            if i in true_set:
                for node in neighbours[i]:
                    is_white[node] = True

            if (i not in true_set) and (i not in del_set):
                new_neigh[i] = [node for node in neighbours[i] if (node not in true_set and node not in del_set)]

        cls._handle_white_nodes(neighbours, is_white, new_neigh, true_set, del_set)

        find = False
        for i in range(n):
            if is_white[i]:
                continue
            if i in true_set or i in del_set:
                continue
            if len(new_neigh[i]) == 1:
                node = new_neigh[i][0]
                true_set.add(node)
                del_set.add(i)

                for neigh in list(new_neigh[node]):
                    if neigh != i:
                        find = True
                    is_white[neigh] = True
                    if node in new_neigh[neigh]:
                        new_neigh[neigh].remove(node)
                new_neigh[node].clear()

        if find:
            cls._handle_white_nodes(neighbours, is_white, new_neigh, true_set, del_set)

        for i in range(n):
            if is_white[i]:
                continue
            if i in true_set or i in del_set:
                continue
            if len(new_neigh[i]) == 0:
                true_set.add(i)

        at_least_pairs.difference_update({(a, b) for (a, b) in at_least_pairs if (a in true_set or b in true_set)})

        del_set.difference_update(true_set)

    def compute_reduction(self, radius: float):
        neighbours = self._build_neighbours(radius)

        Nc, Nd = set(), set()
        at_least_pairs: Set[Tuple[int, int]] = set()

        self._rule1(neighbours, Nc, Nd)
        self._rule2(neighbours, at_least_pairs, Nc, Nd)
        self._addition_rule(neighbours, at_least_pairs, Nc, Nd)

        cnf_extra = [[self._y(a), self._y(b)] for (a, b) in sorted(at_least_pairs)]

        enabled_centers = set(range(self.n))
        demands = list(range(self.n))
        return Nc, Nd, enabled_centers, demands, cnf_extra

    def _encode_cnf(self, radius: float, encoding: str):
        cnf = CNF()

        Nc, Nd, enabled_centers, demands, cnf_extra = self.compute_reduction(radius)

        for clause in cnf_extra:
            cnf.append(clause)

        for c in Nc:
            cnf.append([self._y(c)])

        for d in Nd:
            cnf.append([-self._y(d)])

        Npp = (enabled_centers - Nc) - Nd

        def covered_by_Nc(u: int) -> bool:
            for c in Nc:
                if self.dist[c][u] <= radius + 1e-12:
                    return True
            return False

        for u in demands:
            if covered_by_Nc(u):
                continue
            allowed = [self._y(c) for c in Npp if self.dist[c][u] <= radius + 1e-12]
            if not allowed:  
                return None, {}
            cnf.append(allowed)

        candidates = sorted(list(Npp))

        # cnf_extra contains clauses [y(a), y(b)] where y(j)=1+j
        pair_clauses = []
        for cl in cnf_extra:
            if len(cl) == 2 and cl[0] > 0 and cl[1] > 0:
                a = cl[0] - 1
                b = cl[1] - 1
                pair_clauses.append((a, b))

        sb_coverage_order(
            self, cnf, radius, candidates, demands, Nc,
            enable_orbit_sb=True,
            orbit_mode="chain",
            pair_clauses=pair_clauses
        )

        bound = self.p - len(Nc)      
        if bound < 0:
            if DEBUG_REDUCTION:
                print(f"[ENCODE-FAIL] radius={radius}: bound {bound} < 0 (p={self.p}, |Nc|={len(Nc)})")
            return None, {}

        if candidates:
            lits = [self._y(j) for j in candidates]

            existing_max = cnf.nv if hasattr(cnf, "nv") else 0
            max_y = max(self._y(j) for j in range(self.n)) if self.n > 0 else existing_max
            top_id = max(existing_max, max_y)

            if encoding == "pysat_kmtotalizer":
                enc_kind = CardEncType.kmtotalizer
                with _PYSAT_CNF_LOCK:
                    amo = CardEnc.atmost(lits=lits, bound=bound, top_id=top_id, encoding=enc_kind)
                    cnf.extend(amo.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))
            if encoding == "pysat_mtotalizer":
                enc_kind = CardEncType.mtotalizer
                with _PYSAT_CNF_LOCK:
                    amo = CardEnc.atmost(lits=lits, bound=bound, top_id=top_id, encoding=enc_kind)
                    cnf.extend(amo.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))
            if encoding == "pysat_totalizer":
                enc_kind = CardEncType.totalizer
                with _PYSAT_CNF_LOCK:
                    amo = CardEnc.atmost(lits=lits, bound=bound, top_id=top_id, encoding=enc_kind)
                    cnf.extend(amo.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))
            elif encoding == "nsc":
                self._encode_atmost_nsc(cnf, lits, bound)
            elif encoding == "pypb_bdd":
                enc_kind = PBEncType.bdd
                with _PYSAT_CNF_LOCK:
                    pbcnf = PBEnc.atmost(
                        lits=lits,
                        weights=[1] * len(lits),
                        bound=bound,
                        top_id=top_id,
                        encoding=enc_kind
                    )
                    cnf.extend(pbcnf.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(pbcnf, "nv", 0))
            elif encoding == "pb_bdd":
                self._encode_atmost_pb2cnf(cnf, lits, bound, top_id)

        info = {
            "Nc": Nc, "Nd": Nd,
            "candidates": candidates,
            "bound": bound,
            "y": [self._y(j) for j in range(self.n)]
        }
        return cnf, info

    def _encode_atmost_nsc(self, cnf, lits, bound):
        n = len(lits)
        k = bound

        if k < 0:
            cnf.append([])
            return
        if k == 0:
            for l in lits:
                cnf.append([-l])
            cnf.nv = max(getattr(cnf, "nv", 0), max(lits) if lits else 0)
            return
        if n == 0 or k >= n:
            cnf.nv = max(getattr(cnf, "nv", 0), max(lits) if lits else 0)
            return

        existing_max = getattr(cnf, "nv", 0)
        max_input = max(lits) if lits else 0
        next_var = max(existing_max, max_input) + 1

        def new_var():
            nonlocal next_var
            v = next_var
            next_var += 1
            return v

        # --- Tạo ma trận R_{i,j} ---
        # R[i][j] ~ "sau khi xem lits[0..i] đã có ít nhất (j+1) biến TRUE"
        # độ dài hàng i: min(k, i+1)
        # chỉ tạo cho i = 0 .. n-1
        R = []
        for i in range(n - 1):
            row_len = min(i + 1, k)
            row = [new_var() for _ in range(row_len)]
            R.append(row)

        # --- (1) X_i -> R_{i,1}  ===  (¬X_i ∨ R_{i,0})
        for i in range(1, n):
            cnf.append([-lits[i-1], R[i-1][0]])

        # --- (2) R_{i-1,j} -> R_{i,j}  ===  (¬R_{i-1,j} ∨ R_{i,j})
        for i in range(2, n):
            prev = R[i-2]
            curr = R[i-1]
            for j in range(len(prev)):
                cnf.append([-prev[j], curr[j]])

        # --- (3) (X_i ∧ R_{i-1,j-1}) -> R_{i,j}
        # === (¬X_i ∨ ¬R_{i-1,j-1} ∨ R_{i,j})
        for i in range(2, n):
            prev = R[i-2]
            curr = R[i-1]
            for j in range(1, len(curr)):
                cnf.append([
                    -lits[i-1],
                    -prev[j-1],
                    curr[j]
                ])


        # --- (8) X_i -> ¬R_{i-1,k} === (¬X_i ∨ ¬R_{i-1,k-1})
        for i in range(k + 1, n + 1):
            cnf.append([
                -lits[i-1],
                -R[i-2][k-1]
            ])

        cnf.nv = max(cnf.nv, next_var - 1)

    def _encode_atmost_pb2cnf(self, cnf, lits, bound, top_id_start):
        lits = [int(x) for x in lits]
        n = len(lits)
        k = int(bound)
        first_free = int(max(getattr(cnf, "nv", 0), max(lits, default=0), int(top_id_start)) + 1)

        if k < 0:
            cnf.append([])
            return
        if k == 0:
            for l in lits:
                cnf.append([-l])
            cnf.nv = max(getattr(cnf, "nv", 0), max(lits, default=0))
            return
        if n == 0 or k >= n:
            cnf.nv = max(getattr(cnf, "nv", 0), max(lits, default=0))
            return

        cfg = PBConfig()
        cfg.set_PB_Encoder(pblib.PB_BDD)

        pb2 = Pb2cnf(cfg)

        formula = [] 
        max_var = pb2.encode_at_most_k(lits, int(k), formula, int(first_free))

        for cls in formula:
            cnf.append([int(v) for v in cls])

        cnf.nv = max(getattr(cnf, "nv", 0), int(max_var))




