import math
from typing import List, Tuple
from threading import RLock
from pysat.formula import CNF
from pysat.card import CardEnc
from pysat.card import EncType as CardEncType

from pysat.pb import PBEnc
from pysat.pb import EncType as PBEncType

from pypblib import pblib
from pypblib.pblib import PBConfig, Pb2cnf

# from reduction import compute_reduction
# from sb import automorphism_symmetry_breaking

_PYSAT_CNF_LOCK = RLock()
DEBUG_REDUCTION = False

class PCenterSAT:
    def __init__(self, dist: List[List[float]], p: int):
        n = len(dist)
        self.n: int = n
        self.p: int = int(p)
        self.dist: List[List[float]] = dist
        self.radii: List[float] = self._compute_radii(dist)
        self.y_lit_all = [self._y(i) for i in range(self.n)] 

        self._incremental_cover_order = None

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
                d = math.hypot(dx, dy)
                D[i][j] = math.floor(100.0 * d + 0.5) / 100.0
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

    def _encode_cnf(self, radius: float, encoding: str):
        cnf = CNF()

        # Nc, Nd, enabled_centers, demands, at_least_pairs = compute_reduction(self.dist, self.n, radius)

        # cnf_extra = [[self._y(a), self._y(b)] for (a, b) in sorted(at_least_pairs)]

        # for clause in cnf_extra:
        #     cnf.append(clause)

        # for c in Nc:
        #     cnf.append([self._y(c)])

        # for d in Nd:
        #     cnf.append([-self._y(d)])

        # Npp = (enabled_centers - Nc) - Nd
        # candidates = sorted(list(Npp))

        # covered = [False] * self.n
        # if Nc:
        #     # only check u in demands
        #     for c in Nc:
        #         row = self.dist[c]
        #         for u in demands:
        #             if not covered[u] and row[u] <= radius + 1e-12:
        #                 covered[u] = True
        
        # active_demands = [u for u in demands if not covered[u]]
        active_demands = list(range(self.n))
        candidates = list(range(self.n))

        # ---- demand coverage clauses using candidates list (avoid set scans)
        for u in active_demands:
            allowed = []
            # scan candidates once
            for c in candidates:
                if self.dist[c][u] <= radius + 1e-12:
                    allowed.append(self.y_lit_all[c])
            if not allowed:
                return None, {}
            cnf.append(allowed)

        # bound = self.p - len(Nc)      
        bound = self.p
        # if bound < 0:
        #     if DEBUG_REDUCTION:
        #         print(f"[ENCODE-FAIL] radius={radius}: bound {bound} < 0 (p={self.p}, |Nc|={len(Nc)})")
        #     return None, {}
        
        # automorphism_symmetry_breaking(self, cnf, radius, candidates, active_demands, mode="chain")

        if candidates:
            lits = [self.y_lit_all[j] for j in candidates]

            existing_max = getattr(cnf, "nv", 0)
            max_y = self.y_lit_all[-1] if self.n > 0 else existing_max
            top_id = max(existing_max, max_y)

            if encoding == "pysat_kmtotalizer":
                enc_kind = CardEncType.kmtotalizer
                with _PYSAT_CNF_LOCK:
                    amo = CardEnc.atmost(lits=lits, bound=bound, top_id=top_id, encoding=enc_kind)
                    cnf.extend(amo.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))
            elif encoding == "pysat_mtotalizer":
                enc_kind = CardEncType.mtotalizer
                with _PYSAT_CNF_LOCK:
                    amo = CardEnc.atmost(lits=lits, bound=bound, top_id=top_id, encoding=enc_kind)
                    cnf.extend(amo.clauses)
                    cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))
            elif encoding == "pysat_totalizer":
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
            # "Nc": Nc, "Nd": Nd,
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
    
    def _new_aux_var(self, cnf):
        base = self.y_lit_all[-1] if self.n > 0 else 0
        cnf.nv = max(getattr(cnf, "nv", 0), base) + 1
        return cnf.nv
    
    def _prepare_incremental_cover_order(self):
        """
        Precompute, for each demand u, the centers c sorted by distance dist[c][u].
        Since reduction is OFF, candidates = all vertices and active_demands = all vertices.
        """
        if self._incremental_cover_order is not None:
            return

        cover_order = []
        for u in range(self.n):
            arr = []
            for c in range(self.n):
                arr.append((self.dist[c][u], c))
            arr.sort(key=lambda x: x[0])
            cover_order.append(arr)

        self._incremental_cover_order = cover_order

    def _build_incremental_base_cnf(self, encoding: str):
        """
        Build the permanent part of incremental SAT:
          - y variables
          - one global AtMost-p over all y
        No coverage clauses here.
        """
        cnf = CNF()

        candidates = list(range(self.n))
        active_demands = list(range(self.n))
        lits = [self.y_lit_all[j] for j in candidates]
        bound = self.p

        max_y = self.y_lit_all[-1] if self.n > 0 else 0
        cnf.nv = max(getattr(cnf, "nv", 0), max_y)
        top_id = cnf.nv

        if encoding == "pysat_kmtotalizer":
            enc_kind = CardEncType.kmtotalizer
            with _PYSAT_CNF_LOCK:
                amo = CardEnc.atmost(
                    lits=lits, bound=bound, top_id=top_id, encoding=enc_kind
                )
                cnf.extend(amo.clauses)
                cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))

        elif encoding == "pysat_mtotalizer":
            enc_kind = CardEncType.mtotalizer
            with _PYSAT_CNF_LOCK:
                amo = CardEnc.atmost(
                    lits=lits, bound=bound, top_id=top_id, encoding=enc_kind
                )
                cnf.extend(amo.clauses)
                cnf.nv = max(getattr(cnf, "nv", 0), getattr(amo, "nv", 0))

        elif encoding == "pysat_totalizer":
            enc_kind = CardEncType.totalizer
            with _PYSAT_CNF_LOCK:
                amo = CardEnc.atmost(
                    lits=lits, bound=bound, top_id=top_id, encoding=enc_kind
                )
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
                    encoding=enc_kind,
                )
                cnf.extend(pbcnf.clauses)
                cnf.nv = max(getattr(cnf, "nv", 0), getattr(pbcnf, "nv", 0))

        elif encoding == "pb_bdd":
            self._encode_atmost_pb2cnf(cnf, lits, bound, top_id)

        else:
            raise ValueError(f"Unsupported encoding for incremental mode: {encoding}")

        info = {
            "candidates": candidates,
            "active_demands": active_demands,
            "bound": bound,
            "y": [self._y(j) for j in range(self.n)],
            "radii": list(self.radii),
        }
        return cnf, info

    def _build_incremental_guarded_clauses(self, radius: float, selector_lit: int):
        """
        Build guarded coverage clauses for ONE radius:
            (-selector_lit OR y_c1 OR y_c2 OR ...)
        using precomputed sorted cover order.

        This matches the repo idea you sent:
          - AtMost is permanent
          - coverage is added lazily per tested radius
        """
        self._prepare_incremental_cover_order()

        clauses = []
        for u in range(self.n):
            allowed = [-selector_lit]
            for dcu, c in self._incremental_cover_order[u]:
                if dcu <= radius + 1e-12:
                    allowed.append(self.y_lit_all[c])
                else:
                    break
            clauses.append(allowed)

        return clauses
    
    def _build_incremental2_base_cnf(self, encoding: str):
        """
        Build the permanent part of incremental2 SAT:
          - one global AtMost-p over all y
          - unguarded coverage clauses for the largest radius radii[0]

        Later, smaller radii are enforced lazily by guarded coverage clauses
        under assumptions, so this mode keeps the user's idea of having the
        upper-bound coverage already loaded in the solver.
        """
        cnf, info = self._build_incremental_base_cnf(encoding=encoding)

        if not self.radii:
            return cnf, info

        radius_ub = self.radii[0]
        for u in range(self.n):
            allowed = []
            for c in range(self.n):
                if self.dist[c][u] <= radius_ub + 1e-12:
                    allowed.append(self.y_lit_all[c])
            if not allowed:
                cnf.append([])
            else:
                cnf.append(allowed)

        info = dict(info)
        info["base_radius"] = radius_ub
        return cnf, info
        




