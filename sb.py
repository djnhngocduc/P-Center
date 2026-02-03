from collections import defaultdict

def sb_identical_neighbourhood(self, cnf, radius, candidates, demands, Nc):
    """
    Safe symmetry breaking:
    apply ordering only to centers with truly identical coverage
    AFTER reduction (Nc fixed).
    """

    signature = defaultdict(list)

    # chỉ xét demand CHƯA được cover bởi Nc
    active_demands = []
    for u in demands:
        covered = False
        for c in Nc:
            if self.dist[c][u] <= radius + 1e-12:
                covered = True
                break
        if not covered:
            active_demands.append(u)

    # build exact coverage signature
    for c in candidates:
        covered = []
        for u in active_demands:
            if self.dist[c][u] <= radius + 1e-12:
                covered.append(u)
        signature[tuple(covered)].append(c)

    # lex-chain ONLY inside true equivalence classes
    for group in signature.values():
        if len(group) <= 1:
            continue

        group.sort()
        for i in range(1, len(group)):
            # y[group[i]] -> y[group[i-1]]
            cnf.append([
                -self._y(group[i]),
                 self._y(group[i - 1])
            ])
