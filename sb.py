def sb_smallest_k(self, cnf, candidates, bound):
    """
    SAFE & STRONG symmetry breaking:
    only allow selecting among the 'bound' smallest candidate indices.
    """
    if bound <= 0:
        return

    candidates = sorted(candidates)

    # forbid selecting larger indices
    for c in candidates[bound:]:
        cnf.append([-self._y(c)])
