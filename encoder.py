import math
from typing import List, Tuple
from pysat.formula import CNF, WCNF

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

        active_demands = list(range(self.n))
        candidates = list(range(self.n))

        self._prepare_incremental_cover_order()

        for u in active_demands:
            allowed = []
            for dcu, c in self._incremental_cover_order[u]:
                if dcu <= radius + 1e-12:
                    allowed.append(self.y_lit_all[c])
                else:
                    break
            if not allowed:
                return None, {}
            cnf.append(allowed)

        bound = self.p

        if candidates:
            lits = [self.y_lit_all[j] for j in candidates]

            if encoding == "nsc":
                self._encode_atmost_nsc(cnf, lits, bound)
            elif encoding == "sc":
                self._encode_atmost_sc(cnf, lits, bound)
            else:
                raise ValueError(f"Unsupported SAT encoding: {encoding}")
            

        info = {
            "candidates": candidates,
            "bound": bound,
            "y": [self._y(j) for j in range(self.n)]
        }
        return cnf, info
    
    def _encode_wcnf_maxsat(self, radius: float):
        wcnf = WCNF()

        active_demands = list(range(self.n))
        candidates = list(range(self.n))

        self._prepare_incremental_cover_order()

        for u in active_demands:
            allowed = []
            for dcu, c in self._incremental_cover_order[u]:
                if dcu <= radius + 1e-12:
                    allowed.append(self.y_lit_all[c])
                else:
                    break

            if not allowed:
                return None, {}

            wcnf.append(allowed)

        for j in candidates:
            wcnf.append([-self.y_lit_all[j]], weight=1)

        wcnf.nv = max(
            getattr(wcnf, "nv", 0),
            self.y_lit_all[-1] if self.n > 0 else 0,
        )

        info = {
            "candidates": candidates,
            "active_demands": active_demands,
            "bound": self.p,
            "y": [self._y(j) for j in range(self.n)],
            "radius": radius,
        }
        return wcnf, info

    def _build_setcover_data(self, radius: float):
        active_demands = list(range(self.n))
        candidates = list(range(self.n))

        self._prepare_incremental_cover_order()

        cover_rows = []
        for u in active_demands:
            allowed = []
            for dcu, c in self._incremental_cover_order[u]:
                if dcu <= radius + 1e-12:
                    allowed.append(c)
                else:
                    break
                
            if not allowed:
                return None
                
            cover_rows.append((u, allowed))

        return {
            "n": self.n,
            "p": self.p,
            "candidates": candidates,
            "active_demands": active_demands,
            "cover_rows": cover_rows,
        }

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

        R = []
        for i in range(n - 1):
            row_len = min(i + 1, k)
            row = [new_var() for _ in range(row_len)]
            R.append(row)

        for i in range(1, n):
            cnf.append([-lits[i-1], R[i-1][0]])

        for i in range(2, n):
            prev = R[i-2]
            curr = R[i-1]
            for j in range(len(prev)):
                cnf.append([-prev[j], curr[j]])

        for i in range(2, n):
            prev = R[i-2]
            curr = R[i-1]
            for j in range(1, len(curr)):
                cnf.append([
                    -lits[i-1],
                    -prev[j-1],
                    curr[j]
                ])

        for i in range(k + 1, n + 1):
            cnf.append([
                -lits[i-1],
                -R[i-2][k-1]
            ])

        cnf.nv = max(cnf.nv, next_var - 1)

    def _encode_atmost_sc(self, cnf, lits, bound):
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

        S = []
        for _ in range(n - 1):
            S.append([new_var() for _ in range(k)])

        cnf.append([-lits[0], S[0][0]])
        for j in range(1, k):
            cnf.append([-S[0][j]])

        for i in range(1, n - 1):
            prev = S[i - 1]
            curr = S[i]
            x = lits[i]

            cnf.append([-x, curr[0]])
            cnf.append([-prev[0], curr[0]])

            for j in range(1, k):
                cnf.append([-x, -prev[j - 1], curr[j]])
                cnf.append([-prev[j], curr[j]])

            cnf.append([-x, -prev[k - 1]])

        cnf.append([-lits[n - 1], -S[n - 2][k - 1]])
        cnf.nv = max(getattr(cnf, "nv", 0), next_var - 1)

    def _prepare_incremental_cover_order(self):
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
        cnf = CNF()

        candidates = list(range(self.n))
        active_demands = list(range(self.n))
        lits = [self.y_lit_all[j] for j in candidates]
        bound = self.p

        max_y = self.y_lit_all[-1] if self.n > 0 else 0
        cnf.nv = max(getattr(cnf, "nv", 0), max_y)

        if encoding == "nsc":
            self._encode_atmost_nsc(cnf, lits, bound)

        elif encoding == "sc":
            self._encode_atmost_sc(cnf, lits, bound)

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
        




