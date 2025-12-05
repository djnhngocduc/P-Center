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


_PYSAT_CNF_LOCK = RLock()
DEBUG_REDUCTION = False


class PCenterSAT:
    def __init__(self, dist: List[List[float]], p: int):
        n = len(dist)
        self.n: int = n
        self.p: int = int(p)
        self.dist: List[List[float]] = dist
        self.radii: List[float] = self._compute_radii(dist)

    # ----- variable map -----
    def _y(self, j: int) -> int:
        # 1..n
        return 1 + j

    # ----- factories -----
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

    # ----- radii -----
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
            # i != j; j tăng dần -> list tự sort như C++
            for j in range(n):
                if i == j:
                    continue
                if row[j] <= radius + eps:
                    neigh[i].append(j)
        return neigh

    # --------- helpers ----------
    @staticmethod
    def _dominates(belong_to: List[Set[int]], u: int, targets: List[int]) -> bool:
        # all t in targets are neighbors of u
        s = belong_to[u]
        for t in targets:
            if t == u:
                continue
            if t not in s:
                return False
        return True

    @staticmethod
    def _bfs_within_3(neighbours: List[List[int]], start: int) -> Set[int]:
        # nodes with graph distance <= 3 from start (excluding start)
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
        # mimic Solver.cpp reduction_merge (including its leftover behavior)
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

        # C++ version appends leftovers without skipping v/w
        while i < len(nv):
            out.append(nv[i]); i += 1
        while j < len(nw):
            out.append(nw[j]); j += 1
        return out

    # --------- Rule 1 (match Solver.cpp) ----------
    @staticmethod
    def _rule1(neighbours: List[List[int]], true_set: Set[int], del_set: Set[int]) -> None:
        n = len(neighbours)
        belong = [set(neighbours[v]) for v in range(n)]

        # compute N1/N2/N3 then apply (like C++, but we can apply on the fly safely)
        for v in range(n):
            if v in del_set:
                continue

            Nv = neighbours[v]

            # N1(v): u in N(v) that can dominate v
            N1 = []
            for u in Nv:
                dominated = False
                for t in neighbours[u]:
                    # t not in N(v) and t != v
                    if t != v and (t not in belong[v]):
                        dominated = True
                        break
                if dominated:
                    N1.append(u)
            N1_set = set(N1)

            # N2(v): u in N(v)\N1(v) that has a neighbour in N1(v)
            N2 = []
            for u in Nv:
                if u in N1_set:
                    continue
                for t in neighbours[u]:
                    if t in N1_set:
                        N2.append(u)
                        break
            N2_set = set(N2)

            # N3(v) = N(v) \ (N1 ∪ N2)
            if len(Nv) == (len(N1) + len(N2)):
                continue
            N3 = [u for u in Nv if (u not in N1_set and u not in N2_set)]
            if N3 and (v not in del_set):
                true_set.add(v)
                for u in N3:
                    del_set.add(u)
                for u in N2:
                    del_set.add(u)

        # keep invariant: true ∩ del = ∅
        del_set.difference_update(true_set)

    # --------- Rule 2 (match Solver.cpp behavior) ----------
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
            # mimic (v < w)
            for w in sorted(x for x in within3 if x > v):
                if w in del_set:
                    continue

                N_vw = PCenterSAT._reduction_merge(neighbours, v, w)

                # N1(v,w): u in N(v,w) can be dominated (per C++ ifCanBeDominated)
                N1 = []
                N_vw_set = set(N_vw)
                for u in N_vw:
                    dominated = False
                    for t in neighbours[u]:
                        # t not in N(v,w) and t != v and t != w
                        if t != v and t != w and (t not in N_vw_set):
                            dominated = True
                            break
                    if dominated:
                        N1.append(u)
                N1_set = set(N1)

                # N2(v,w): u in N(v,w)\N1 that has neighbour in N1
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

                # N3 = N_vw \ (N1 ∪ N2)
                N3 = [u for u in N_vw if (u not in N1_set and u not in N2_set)]

                if len(N3) <= 1:
                    continue

                # dominated check: if exists u in N3 or N2 that dominates N3 => skip
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
                        # add (v OR w) only if both not already true
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
                    # neither dominates: force both v,w true; delete N3 and N2 but don't delete any true nodes
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

        # delete edges between white nodes
        for i in range(n):
            if not is_white[i]:
                continue
            new_neigh[i] = [v for v in new_neigh[i] if not is_white[v]]

        # white nodes degree 0 or 1
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

        # mark white nodes = neighbours of forced centers
        for i in range(n):
            if i in true_set:
                for node in neighbours[i]:
                    is_white[node] = True

            if (i not in true_set) and (i not in del_set):
                new_neigh[i] = [node for node in neighbours[i] if (node not in true_set and node not in del_set)]

        cls._handle_white_nodes(neighbours, is_white, new_neigh, true_set, del_set)

        find = False
        # degree-1 nodes (non-white)
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
                    # erase node from neigh adjacency
                    if node in new_neigh[neigh]:
                        new_neigh[neigh].remove(node)
                new_neigh[node].clear()

        if find:
            cls._handle_white_nodes(neighbours, is_white, new_neigh, true_set, del_set)

        # degree-0 nodes (non-white) => must be center
        for i in range(n):
            if is_white[i]:
                continue
            if i in true_set or i in del_set:
                continue
            if len(new_neigh[i]) == 0:
                true_set.add(i)

        # drop satisfied (v OR w) where any endpoint already true
        at_least_pairs.difference_update({(a, b) for (a, b) in at_least_pairs if (a in true_set or b in true_set)})

        del_set.difference_update(true_set)

    def compute_reduction(self, radius: float):
        neighbours = self._build_neighbours(radius)

        Nc, Nd = set(), set()
        at_least_pairs: Set[Tuple[int, int]] = set()

        # dataReduction = Rule1 -> Rule2 -> AdditionRule  (exact order)
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

        # clauses from Rule 2
        for clause in cnf_extra:
            cnf.append(clause)

        # force centers in Nc
        for c in Nc:
            cnf.append([self._y(c)])

        # forbid centers in Nd
        for d in Nd:
            cnf.append([-self._y(d)])

        # remaining candidates (can be opened)
        Npp = (enabled_centers - Nc) - Nd

        # cover constraints
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

        # at-most p - |Nc|
        candidates = sorted(list(Npp))
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

            if encoding == "pysat_sc":
                enc_kind = CardEncType.seqcounter
                with _PYSAT_CNF_LOCK:
                    amo = CardEnc.atmost(lits=lits,
                                         bound=bound, top_id=top_id, encoding=enc_kind)
                    cnf.extend(amo.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))
            elif encoding == "nsc":
                self._encode_atmost_nsc(cnf, lits, bound)
            elif encoding == "pypb_sc":
                enc_kind = PBEncType.seqcounter
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
        """
        NSC<=k (new sequential counter for AtMostK)
        ràng buộc: sum(lits) <= bound

        lits: list[int] các biến dương (y_j)
        bound: int k

        Sẽ:
        - tạo biến phụ R_{i,j} với i = 0..n-1 (prefix đến i), j = 0..row_len-1 (ít nhất j+1 true)
        - thêm các clause (1)(2)(3)(8) như trong NSC≤k
        - cập nhật cnf.nv
        """
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
        for i in range(n):
            row_len = min(k, i + 1)  # NSC giảm số cột cho prefix ngắn
            row = [new_var() for _ in range(row_len)]
            R.append(row)

        # --- (1) X_i -> R_{i,1}  ===  (¬X_i ∨ R_{i,0})
        # áp cho i = 0..n-1 vì R chỉ định nghĩa tới n-1
        for i in range(n):
            if len(R[i]) >= 1:
                cnf.append([-lits[i], R[i][0]])

        # --- (2) R_{i-1,j} -> R_{i,j}  ===  (¬R_{i-1,j} ∨ R_{i,j})
        for i in range(1, n):
            # j chạy qua những cột chung giữa R[i-1] và R[i]
            lim = min(len(R[i - 1]), len(R[i]))
            for j in range(lim):
                cnf.append([-R[i - 1][j], R[i][j]])

        # --- (3) (X_i ∧ R_{i-1,j-1}) -> R_{i,j}
        # === (¬X_i ∨ ¬R_{i-1,j-1} ∨ R_{i,j})
        for i in range(1, n):
            for j in range(1, len(R[i])):
                if j - 1 < len(R[i - 1]):
                    cnf.append([
                        -lits[i],
                        -R[i - 1][j - 1],
                        R[i][j]
                    ])

        # --- (8) X_i -> ¬R_{i-1,k} === (¬X_i ∨ ¬R_{i-1,k-1})
        # chạy với i = k .. n-1 (chú ý index)
        # i ở đây là index trong lits, R prefix dùng i-1
        for i in range(k, n):
            if (k - 1) < len(R[i - 1]):
                cnf.append([
                    -lits[i],
                    -R[i - 1][k - 1]
                ])

        # cập nhật cnf.nv
        cnf.nv = max(getattr(cnf, "nv", 0), next_var - 1)

    def _encode_atmost_pb2cnf(self, cnf, lits, bound, top_id_start):
        """
        Encode ràng buộc sum(lits) <= bound bằng PyPBLib (Pb2cnf) với encoder BDD.
        Dùng đúng API: encode_at_most_k(literals, k, formula, first_free_var).
        """
        # Chuẩn hoá & các trường hợp biên
        lits = [int(x) for x in lits]
        n = len(lits)
        k = int(bound)
        first_free = int(max(getattr(cnf, "nv", 0), max(lits, default=0), int(top_id_start)) + 1)

        if k < 0:
            # vô nghiệm
            cnf.append([])  # mệnh đề rỗng
            return
        if k == 0:
            # at-most 0 -> tất cả đều âm
            for l in lits:
                cnf.append([-l])
            cnf.nv = max(getattr(cnf, "nv", 0), max(lits, default=0))
            return
        if n == 0 or k >= n:
            # không cần mã hoá gì thêm
            cnf.nv = max(getattr(cnf, "nv", 0), max(lits, default=0))
            return

        cfg = PBConfig()
        cfg.set_PB_Encoder(pblib.PB_BDD) 

        pb2 = Pb2cnf(cfg)

        formula = [] 
        max_var = pb2.encode_at_most_k(lits, int(k), formula, int(first_free))

        # nối vào CNF của PySAT
        for cls in formula:
            cnf.append([int(v) for v in cls])

        # cập nhật số biến tối đa
        cnf.nv = max(getattr(cnf, "nv", 0), int(max_var))




