def sb_coverage_order(self, cnf, radius, candidates, demands, Nc):
    """
    STRONG & SAFE dominance-based symmetry breaking:
    enforce ordering only when coverage inclusion holds.
    """

    # active demands = chưa được cover bởi Nc
    active_demands = [
        u for u in demands
        if not any(self.dist[c][u] <= radius + 1e-12 for c in Nc)
    ]

    # build coverage sets
    cover = {}
    for c in candidates:
        cover[c] = {
            u for u in active_demands
            if self.dist[c][u] <= radius + 1e-12
        }

    # dominance-based ordering
    for i, ci in enumerate(candidates):
        for cj in candidates[i+1:]:
            if cover[ci] >= cover[cj] and cover[ci] != cover[cj]:
                # cj -> ci
                cnf.append([
                    -self._y(cj),
                     self._y(ci)
                ])
            elif cover[cj] >= cover[ci] and cover[cj] != cover[ci]:
                # ci -> cj
                cnf.append([
                    -self._y(ci),
                     self._y(cj)
                ])
