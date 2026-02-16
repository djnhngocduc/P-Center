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

from reduction import compute_reduction

# from sb import sb_coverage_order, orbit_symmetry_breaking
from sb import automorphism_symmetry_breaking

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

    def _encode_cnf(self, radius: float, encoding: str):
        cnf = CNF()

        Nc, Nd, enabled_centers, demands, at_least_pairs = compute_reduction(self.dist, self.n, radius)

        cnf_extra = [[self._y(a), self._y(b)] for (a, b) in sorted(at_least_pairs)]

        for clause in cnf_extra:
            cnf.append(clause)

        for c in Nc:
            cnf.append([self._y(c)])

        for d in Nd:
            cnf.append([-self._y(d)])

        Npp = (enabled_centers - Nc) - Nd
        candidates = sorted(list(Npp))

        covered = [False] * self.n
        if Nc:
            # only check u in demands
            for c in Nc:
                row = self.dist[c]
                for u in demands:
                    if not covered[u] and row[u] <= radius + 1e-12:
                        covered[u] = True

        # ---- demand coverage clauses using candidates list (avoid set scans)
        for u in demands:
            if covered[u]:
                continue
            row_u = u  # just name
            allowed = []
            # scan candidates once
            for c in candidates:
                if self.dist[c][row_u] <= radius + 1e-12:
                    allowed.append(self.y_lit_all[c])
            if not allowed:
                return None, {}
            cnf.append(allowed)

        # cnf_extra contains clauses [y(a), y(b)] where y(j)=1+j
        # pair_clauses = []
        # for cl in cnf_extra:
        #     if len(cl) == 2 and cl[0] > 0 and cl[1] > 0:
        #         a = cl[0] - 1
        #         b = cl[1] - 1
        #         pair_clauses.append((a, b))


        # mask2cands, dom_out, dom_in = sb_coverage_order(
        #     self, cnf, radius, candidates, demands, Nc
        # )

        # orbit_symmetry_breaking(
        #     self,
        #     cnf,
        #     mask2cands=mask2cands,
        #     dom_out=dom_out,
        #     dom_in=dom_in,
        #     pair_clauses=pair_clauses,
        #     orbit_mode="chain",   # hoặc "leader"
        # )
        automorphism_symmetry_breaking(
            self,
            cnf,
            radius,
            Nc, Nd, enabled_centers, demands,
            mode="leader",  # hoặc "leader"
        )

        bound = self.p - len(Nc)      
        if bound < 0:
            if DEBUG_REDUCTION:
                print(f"[ENCODE-FAIL] radius={radius}: bound {bound} < 0 (p={self.p}, |Nc|={len(Nc)})")
            return None, {}

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




