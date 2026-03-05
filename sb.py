from collections import defaultdict

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

def _cnf_to_pynauty_graph_full(cnf, cand_vars):
    import pynauty

    # Collect all vars appearing in CNF (including aux vars)
    var_set = set()
    for cls in cnf.clauses:
        for lit in cls:
            var_set.add(abs(lit))
    var_lits = sorted(var_set)

    # literal nodes
    pos_id, neg_id = {}, {}
    nodes = 0
    for v in var_lits:
        pos_id[v] = nodes; nodes += 1
        neg_id[v] = nodes; nodes += 1

    # clause nodes
    clause_nodes = []
    for _ in cnf.clauses:
        clause_nodes.append(nodes)
        nodes += 1

    adj = {i: [] for i in range(nodes)}
    def add_edge(a, b):
        adj[a].append(b); adj[b].append(a)

    # bind pos/neg
    for v in var_lits:
        add_edge(pos_id[v], neg_id[v])

    # incidence edges
    for cnode, cls in zip(clause_nodes, cnf.clauses):
        for lit in cls:
            v = abs(lit)
            lnode = pos_id[v] if lit > 0 else neg_id[v]
            add_edge(cnode, lnode)

    # ---- Coloring (CRITICAL) ----
    # Split literals by candidate vars vs others, and by sign
    pos_cand, pos_other, neg_cand, neg_other = set(), set(), set(), set()
    for v in var_lits:
        if v in cand_vars:
            pos_cand.add(pos_id[v]); neg_cand.add(neg_id[v])
        else:
            pos_other.add(pos_id[v]); neg_other.add(neg_id[v])

    coloring = []
    if pos_cand:  coloring.append(pos_cand)
    if pos_other: coloring.append(pos_other)
    if neg_cand:  coloring.append(neg_cand)
    if neg_other: coloring.append(neg_other)

    # Clause coloring by length only (full CNF already distinguishes structure via neighbors)
    by_len = defaultdict(list)
    for cnode, cls in zip(clause_nodes, cnf.clauses):
        by_len[len(cls)].append(cnode)
    for L in sorted(by_len):
        coloring.append(set(by_len[L]))
    # -----------------------------

    g = pynauty.Graph(
        number_of_vertices=nodes,
        adjacency_dict=adj,
        vertex_coloring=coloring,
    )
    return g, pos_id, nodes

def compute_candidate_orbits_from_full_cnf(cnf, candidate_center_indices, y_lit_all):
    import pynauty

    cand_vars = set(y_lit_all[c] for c in candidate_center_indices)  # SAT var ids
    if len(cand_vars) <= 1:
        return []

    g, pos_id, n_vertices = _cnf_to_pynauty_graph_full(cnf, cand_vars)
    result = pynauty.autgrp(g)
    orbit_ids = _get_orbit_ids_from_autgrp(result, n_vertices)

    # orbit by variable using pos literal node id
    orbit_dict = defaultdict(list)
    for v in cand_vars:
        oid = orbit_ids[pos_id[v]]
        orbit_dict[oid].append(v)

    # map back to center indices: v = 1 + c => c = v - 1
    orbits = []
    for vars_ in orbit_dict.values():
        grp = sorted([v - 1 for v in vars_])
        if len(grp) >= 2:
            orbits.append(grp)
    orbits.sort()
    return orbits

def automorphism_symmetry_breaking(self, cnf, candidates, mode="chain"):
    ylit = self.y_lit_all
    orbits = compute_candidate_orbits_from_full_cnf(cnf, candidates, ylit)
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