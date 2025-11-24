import math
from typing import List, Tuple, Set
from threading import RLock
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
        INF = 10 ** 14
        for i in range(n):
            for j in range(i + 1, n):
                dij = D[i][j]
                if dij != 0 and dij < INF:
                    vals.add(dij)
        return sorted(vals, reverse=True)

    # ----- graph induced by radius -----
    def _induced_graph(self, radius: float):
        n = self.n
        adj = [set() for _ in range(n)]
        for v in range(n):
            row = self.dist[v]
            for u in range(n):
                if v != u and row[u] <= radius + 1e-12:
                    adj[v].add(u)
        return adj

    # ----- O-set partitions (Alber) -----
    def _partition_O_sets_single(self, adj, v):
        """For Rule 1: O-sets wrt N'(v) on G_R."""
        Nv = set(adj[v])              # N'(v)
        Nclosed = Nv | {v}            # N'[v]
        O1 = {u for u in Nv if any((w not in Nclosed) for w in adj[u])}
        O2 = set()
        for u in Nv:
            if u in O1:
                continue
            adj_u = adj[u]
            if any((w in O1) for w in adj_u) and all((w in Nclosed) for w in adj_u):
                O2.add(u)
        O3 = Nv - O1 - O2
        return O1, O2, O3

    def _partition_O_sets_pair(self, adj, v, w):
        """For Rule 2: O-sets wrt N'(v) ∪ N'(w) on G_R."""
        Nv = set(adj[v])
        Nw = set(adj[w])
        U = Nv | Nw                   # (no v,w inside)
        Nvw_closed = U | {v, w}       # N'[v,w]

        # 1) Tính O1 trước (độc lập)
        O1 = {u for u in U if any((x not in Nvw_closed) for x in adj[u])}

        # 2) Tính O2 sau khi O1 đã cố định
        O2 = set()
        for u in U:
            if u in O1:
                continue
            adj_u = adj[u]
            if any((x in O1) for x in adj_u) and all((x in Nvw_closed) for x in adj_u):
                O2.add(u)

        # 3) Phần còn lại là O3
        O3 = U - O1 - O2
        return O1, O2, O3, Nv, Nw


    # ----- reduction: Alber Rule 1 & 2 -----
    def compute_reduction(self, radius: float):
        n = self.n
        Nc, Nd = set(), set()      # forced centers; forbidden centers
        cnf_extra = []             # clauses from Rule 2 (“yv or yw”)

        # G_R: induced by radius
        adj = self._induced_graph(radius)

        def forbid_many(nodes: Set[int]):
            changed_local = False
            for u in nodes:
                if u not in Nc and u not in Nd:
                    Nd.add(u)
                    changed_local = True
            return changed_local

        def covers_all_closed(u: int, subset: Set[int]) -> bool:
            for x in subset:
                if x != u and x not in adj[u]:
                    return False
            return True

        changed = True
        while changed:
            changed = False

            # ---- Rule 1 ----
            for v in range(n):
                if v in Nd or v in Nc:
                    continue
                O1, O2, O3 = self._partition_O_sets_single(adj, v)
                if O3:
                    Nc.add(v)
                    forbid_many(O2 | O3)
                    changed = True
                    break  # restart

            if changed:
                continue

            # ---- Rule 2 ----
            for v in range(n):
                if v in Nd:
                    continue
                for w in range(v + 1, n):
                    if w in Nd:
                        continue
                    O1, O2, O3, Nv, Nw = self._partition_O_sets_pair(adj, v, w)
                    if not O3:
                        continue

                    any_u_covers = False
                    for u in (O2 | O3):
                        if covers_all_closed(u, O3):
                            any_u_covers = True
                            break
                    if any_u_covers:
                        continue

                    can_v = O3.issubset(Nv)
                    can_w = O3.issubset(Nw)

                    if can_v and can_w:
                        # at least one of v,w must be open
                        cnf_extra.append([self._y(v), self._y(w)])
                        # forbid O3 and (O2 ∩ N'(v) ∩ N'(w))
                        changed |= forbid_many(O3 | (O2 & Nv & Nw))
                        if changed:
                            break

                    elif can_v and not can_w:
                        # force v; forbid O3 and O2 ∩ N'(v)
                        if v not in Nc:
                            Nc.add(v)
                            changed = True
                        changed |= forbid_many(O3 | (O2 & Nv))
                        if changed:
                            break

                    elif can_w and not can_v:
                        # symmetric
                        if w not in Nc:
                            Nc.add(w)
                            changed = True
                        changed |= forbid_many(O3 | (O2 & Nw))
                        if changed:
                            break

                    else:
                        # none covers O3 alone -> forbid ONLY O2
                        changed |= forbid_many(O2)
                        if changed:
                            break
                if changed:
                    break

        enabled_centers = set(range(n))
        demands = list(range(n))
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
            # elif encoding == "par":
            #     self._encode_atmost_adder(cnf, [self._y(j) for j in candidates], bound)
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

    # # ----- parallel counter (optional) -----
    # def _ha_encode(self, cnf, a, b, s, c):
    #     cnf.append([ a,  b,  s]); cnf.append([-a,  b, -s])
    #     cnf.append([ a, -b, -s]); cnf.append([-a, -b,  s])
    #     cnf.append([-a, -b,  c]); cnf.append([-c, a]); cnf.append([-c, b])

    # def _fa_encode(self, cnf, a, b, cin, s, cout):
    #     cnf.append([ a,  b, -cin,  s]); cnf.append([-a,  b,  cin,  s])
    #     cnf.append([ a, -b,  cin,  s]); cnf.append([-a, -b, -cin,  s])
    #     cnf.append([-a, -b,  cout]);   cnf.append([-a, -cin, cout]); cnf.append([-b, -cin, cout])
    #     cnf.append([ -cout,  a,  b]);  cnf.append([ -cout,  a,  cin]); cnf.append([ -cout,  b,  cin])

    # def _encode_atmost_adder(self, cnf, lits, bound):
    #     if bound >= len(lits):
    #         return
    #     existing_max = cnf.nv if hasattr(cnf, "nv") else 0
    #     max_y = max(lits) if lits else 0
    #     next_var = max(existing_max, max_y) + 1

    #     def new_var():
    #         nonlocal next_var
    #         v = next_var
    #         next_var += 1
    #         return v

    #     from collections import defaultdict
    #     bucket = defaultdict(list)
    #     for l in lits:
    #         bucket[1].append(l)

    #     w = 1
    #     while True:
    #         B = bucket[w]
    #         if len(B) > 1:
    #             newB = []
    #             while len(B) >= 2:
    #                 a = B.pop()
    #                 b = B.pop()
    #                 s = new_var()
    #                 c = new_var()
    #                 self._ha_encode(cnf, a, b, s, c)
    #                 newB.append(s)
    #                 bucket[2 * w].append(c)
    #             if B:
    #                 newB.append(B.pop())
    #             bucket[w] = newB
    #         w *= 2
    #         if w > 2 * len(lits):
    #             break

    #     bits = sorted([(bucket[p][0], p) for p in bucket if bucket[p]], key=lambda x: x[1], reverse=True)

    #     b = bound

    #     def bit_of(val, pw):
    #         return 1 if (val // pw) % 2 == 1 else 0

    #     higher = []
    #     for (lit, pw) in bits:
    #         if bit_of(b, pw) == 0:
    #             clause = []
    #             for (hlit, hpw) in higher:
    #                 clause.append(hlit if bit_of(b, hpw) == 1 else -hlit)
    #             clause.append(-lit)
    #             cnf.append(clause)
    #         higher.append((lit, pw))

    #     cnf.nv = max(cnf.nv, next_var - 1)

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




