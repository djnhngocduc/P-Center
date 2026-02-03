from collections import defaultdict

def _sb_identical_neighbourhood(self, cnf, radius, Npp, demands):
    """
    Symmetry breaking for centers with identical neighbourhoods.
    Apply local lex ordering inside each equivalence class.
    """
    signature = defaultdict(list)

    # build neighbourhood signature
    for c in Npp:
        covered = []
        for u in demands:
            if self.dist[c][u] <= radius + 1e-12:
                covered.append(u)
        signature[tuple(covered)].append(c)

    # local lex-leader inside each group
    for group in signature.values():
        if len(group) <= 1:
            continue
        group.sort()
        for i in range(1, len(group)):
            # ¬y_group[i] ∨ y_group[i-1]
            cnf.append([
                -self._y(group[i]),
                 self._y(group[i - 1])
            ])

def _sb_component_lex(self, cnf, radius, candidates):
    """
    Symmetry breaking by connected components:
    apply lex ordering only inside each connected component.
    """
    visited = set()

    for v in candidates:
        if v in visited:
            continue

        # BFS to build component
        stack = [v]
        comp = []

        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            comp.append(u)

            for w in candidates:
                if w not in visited and self.dist[u][w] <= radius + 1e-12:
                    stack.append(w)

        if len(comp) <= 1:
            continue

        comp.sort()

        # adjacent lex: y_i -> y_{i-1}
        for i in range(1, len(comp)):
            cnf.append([
                -self._y(comp[i]),
                 self._y(comp[i - 1])
            ])

def _sb_prefix_bound(self, cnf, candidates, bound):
    """
    Cardinality-aware symmetry breaking:
    forbid choosing y_i before y_{i-1} when i > bound.
    """
    if bound <= 0:
        return

    for i in range(bound + 1, len(candidates)):
        cnf.append([
            -self._y(candidates[i]),
             self._y(candidates[i - 1])
        ])

