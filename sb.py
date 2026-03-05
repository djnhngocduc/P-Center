from collections import defaultdict


def _build_prefix_clauses(self, Nc, Nd, at_least_pairs, active_demands, allowed_map):
    typed = []

    for (a, b) in at_least_pairs:
        typed.append(("pair", (self.y_lit_all[a], self.y_lit_all[b])))

    for c in Nc:
        typed.append(("unit_pos", (self.y_lit_all[c],)))

    for d in Nd:
        typed.append(("unit_neg", (-self.y_lit_all[d],)))

    for u in active_demands:
        lits = allowed_map.get(u)
        if not lits:
            return None
        typed.append(("cover", tuple(lits)))

    return typed


def _cnf_prefix_to_pynauty_graph(typed_clauses, var_lits, cand_var_lits):
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

    # -------- Coloring (CRITICAL FIX) --------
    # Split literals by whether the underlying variable is a candidate (counted by AtMost) or not.
    pos_cand = set()
    pos_other = set()
    neg_cand = set()
    neg_other = set()

    for v in var_lits:
        if v in cand_var_lits:
            pos_cand.add(pos_id[v])
            neg_cand.add(neg_id[v])
        else:
            pos_other.add(pos_id[v])
            neg_other.add(neg_id[v])

    coloring = []
    if pos_cand:  coloring.append(pos_cand)
    if pos_other: coloring.append(pos_other)
    if neg_cand:  coloring.append(neg_cand)
    if neg_other: coloring.append(neg_other)

    # Clause colors by (type, length)
    by_kind = defaultdict(list)
    for cnode, (ctype, cls) in zip(clause_nodes, typed_clauses):
        by_kind[(ctype, len(cls))].append(cnode)
    for key in sorted(by_kind.keys(), key=lambda x: (x[0], x[1])):
        coloring.append(set(by_kind[key]))
    # ----------------------------------------

    g = pynauty.Graph(
        number_of_vertices=nodes,
        adjacency_dict=adj,
        vertex_coloring=coloring,
    )
    return g, pos_id, nodes


def _get_orbit_ids_from_autgrp(result, n_vertices: int):
    if not isinstance(result, (tuple, list)) or not result:
        raise TypeError(f"autgrp returned unexpected: {type(result)} {result!r}")

    for item in result:
        if isinstance(item, (list, tuple)) and len(item) == n_vertices and all(isinstance(x, int) for x in item):
            return item

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
    import pynauty  # noqa: F401

    typed_clauses = _build_prefix_clauses(self, Nc, Nd, at_least_pairs, active_demands, allowed_map)
    if typed_clauses is None:
        return []

    # candidates exactly like encoder.py (these are the ones constrained by AtMost)
    candidates = sorted(((enabled_centers - Nc) - Nd))
    cand_set = set(candidates)
    cand_var_lits = set(self.y_lit_all[c] for c in candidates)  # SAT var numbers (1..n)

    # vars in prefix
    var_set = set(cand_var_lits)

    for (a, b) in at_least_pairs:
        var_set.add(self.y_lit_all[a])
        var_set.add(self.y_lit_all[b])
    for c in Nc:
        var_set.add(self.y_lit_all[c])
    for d in Nd:
        var_set.add(self.y_lit_all[d])

    for u in active_demands:
        for lit in allowed_map.get(u, ()):
            var_set.add(abs(lit))

    var_lits = sorted(var_set)
    if len(var_lits) <= 1:
        return []

    g, pos_id, n_vertices = _cnf_prefix_to_pynauty_graph(typed_clauses, var_lits, cand_var_lits)
    result = pynauty.autgrp(g)
    orbit_ids = _get_orbit_ids_from_autgrp(result, n_vertices)

    orbit_dict = defaultdict(list)
    for lit in var_lits:
        oid = orbit_ids[pos_id[lit]]
        orbit_dict[oid].append(lit)

    # Return classes as CENTER INDICES among candidates only
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