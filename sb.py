from collections import defaultdict

def _build_prefix_clauses(self, Nc, Nd, at_least_pairs, active_demands, allowed_map):
    """
    Build the exact CNF-prefix used by encoder.py BEFORE AtMost:
      - pair clauses: (y[a] OR y[b])
      - unit clauses: y[c] for c in Nc
      - unit clauses: ¬y[d] for d in Nd
      - coverage clauses for each u in active_demands: allowed_map[u]
    """
    typed = []

    # (1) pair clauses
    for (a, b) in at_least_pairs:
        typed.append(("pair", (self.y_lit_all[a], self.y_lit_all[b])))

    # (2) Nc units
    for c in Nc:
        typed.append(("unit_pos", (self.y_lit_all[c],)))

    # (3) Nd units
    for d in Nd:
        typed.append(("unit_neg", (-self.y_lit_all[d],)))

    # (4) coverage clauses: already computed in encoder
    for u in active_demands:
        lits = allowed_map.get(u)
        if not lits:
            # encoder would have returned None
            return None
        typed.append(("cover", tuple(lits)))

    return typed


def _cnf_prefix_to_pynauty_graph(typed_clauses, var_lits):
    import pynauty

    pos_id, neg_id = {}, {}
    nodes = 0
    for v in var_lits:
        pos_id[v] = nodes; nodes += 1
        neg_id[v] = nodes; nodes += 1

    clause_nodes = []
    for _ in typed_clauses:
        clause_nodes.append(nodes)
        nodes += 1

    adj = {i: [] for i in range(nodes)}

    def add_edge(a, b):
        adj[a].append(b)
        adj[b].append(a)

    # bind pos/neg
    for v in var_lits:
        add_edge(pos_id[v], neg_id[v])

    # clause-literal incidence
    for cnode, (ctype, cls) in zip(clause_nodes, typed_clauses):
        for lit in cls:
            v = abs(lit)
            lnode = pos_id[v] if lit > 0 else neg_id[v]
            add_edge(cnode, lnode)

    # coloring: pos, neg, and clause groups by (type, length)
    coloring = []
    coloring.append(set(pos_id.values()))
    coloring.append(set(neg_id.values()))

    by_kind = defaultdict(list)
    for cnode, (ctype, cls) in zip(clause_nodes, typed_clauses):
        by_kind[(ctype, len(cls))].append(cnode)

    for key in sorted(by_kind.keys(), key=lambda x: (x[0], x[1])):
        coloring.append(set(by_kind[key]))

    g = pynauty.Graph(
        number_of_vertices=nodes,
        adjacency_dict=adj,
        vertex_coloring=coloring,
    )
    return g, pos_id, nodes


def _get_orbit_ids_from_autgrp(result, n_vertices: int):
    """
    Robustly locate orbit-ids array in pynauty.autgrp result by signature:
      - is list/tuple
      - length == n_vertices
      - all elements are int
    This avoids relying on result[-1] / result[-2] which differs across builds.
    """
    if not isinstance(result, (tuple, list)) or len(result) == 0:
        raise TypeError(f"autgrp returned unexpected: {type(result)} {result!r}")

    # Search any element matching the orbit-id vector shape
    for item in result:
        if isinstance(item, (list, tuple)) and len(item) == n_vertices:
            if all(isinstance(x, int) for x in item):
                return item

    # Fallback: sometimes orbit ids can be nested as last element
    last = result[-1]
    if isinstance(last, (list, tuple)) and len(last) == n_vertices and all(isinstance(x, int) for x in last):
        return last

    raise TypeError(f"orbit_ids not found (n_vertices={n_vertices}): {result!r}")


def compute_automorphism_symmetry(
    self,
    radius,
    Nc,
    Nd,
    enabled_centers,
    active_demands,
    at_least_pairs,
    allowed_map,
):
    import pynauty  # noqa: F401 (import required for autgrp)

    typed_clauses = _build_prefix_clauses(
        self, Nc, Nd, at_least_pairs, active_demands, allowed_map
    )
    if typed_clauses is None:
        return []

    # candidates exactly like encoder.py (only used to filter groups)
    candidates = sorted(((enabled_centers - Nc) - Nd))
    cand_set = set(candidates)

    # variables appearing in the prefix (SAT var numbers 1..n)
    var_set = set(self.y_lit_all[c] for c in candidates)
    for (a, b) in at_least_pairs:
        var_set.add(self.y_lit_all[a])
        var_set.add(self.y_lit_all[b])
    for c in Nc:
        var_set.add(self.y_lit_all[c])
    for d in Nd:
        var_set.add(self.y_lit_all[d])

    # include vars from coverage clauses (allowed_map)
    for u in active_demands:
        for lit in allowed_map.get(u, ()):
            var_set.add(abs(lit))

    var_lits = sorted(var_set)
    if len(var_lits) <= 1:
        return []

    g, pos_id, n_vertices = _cnf_prefix_to_pynauty_graph(typed_clauses, var_lits)
    result = pynauty.autgrp(g)
    orbit_ids = _get_orbit_ids_from_autgrp(result, n_vertices)

    # orbit on pos literal node => orbit of variable
    orbit_dict = defaultdict(list)
    for lit in var_lits:
        oid = orbit_ids[pos_id[lit]]
        orbit_dict[oid].append(lit)

    # return classes as CENTER INDICES among candidates
    classes = []
    for orbit_vars in orbit_dict.values():
        grp = [v - 1 for v in orbit_vars if (v - 1) in cand_set]
        if len(grp) >= 2:
            classes.append(sorted(grp))

    classes.sort()
    return classes


def automorphism_symmetry_breaking(
    self,
    cnf,
    radius,
    Nc,
    Nd,
    enabled_centers,
    active_demands,
    at_least_pairs,
    allowed_map,
    mode="chain",
):
    ylit = self.y_lit_all

    orbits = compute_automorphism_symmetry(
        self,
        radius,
        Nc,
        Nd,
        enabled_centers,
        active_demands,
        at_least_pairs,
        allowed_map,
    )

    if not orbits:
        return

    for group in orbits:
        group.sort()
        if mode == "leader":
            leader = group[0]
            for c in group[1:]:
                cnf.append([-ylit[c], ylit[leader]])
        else:  # chain
            for a, b in zip(group, group[1:]):
                cnf.append([-ylit[b], ylit[a]])