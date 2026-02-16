from typing import List, Tuple, Set
from collections import deque
import csv
import os

def _build_neighbours(dist_matrix: List[List[float]], n: int, radius: float) -> List[List[int]]:
    neigh = [[] for _ in range(n)]
    eps = 1e-12
    for i in range(n):
        row = dist_matrix[i]
        for j in range(n):
            if i == j:
                continue
            if row[j] <= radius + eps:
                neigh[i].append(j)
    return neigh

def _dominates(belong_to: List[Set[int]], u: int, targets: List[int]) -> bool:
    s = belong_to[u]
    for t in targets:
        if t == u:
            continue
        if t not in s:
            return False
    return True

def _bfs_within_3(neighbours: List[List[int]], start: int) -> Set[int]:
    q = deque([(start, 0)])
    seen = {start}
    out = set()
    while q:
        u, d = q.popleft()
        if d == 3:
            continue
        for v in neighbours[u]:
            if v in seen:
                continue
            seen.add(v)
            out.add(v)
            q.append((v, d + 1))
    return out

def _reduction_merge(neighbours: List[List[int]], v: int, w: int) -> List[int]:
    nv = neighbours[v]
    nw = neighbours[w]
    i = j = 0
    out = []

    while i < len(nv) and j < len(nw):
        if nv[i] == w:
            i += 1
            continue
        if nw[j] == v:
            j += 1
            continue

        if nv[i] < nw[j]:
            out.append(nv[i]); i += 1
        elif nv[i] == nw[j]:
            out.append(nv[i]); i += 1; j += 1
        else:
            out.append(nw[j]); j += 1

    while i < len(nv):
        out.append(nv[i]); i += 1
    while j < len(nw):
        out.append(nw[j]); j += 1
    return out

def _rule1(neighbours: List[List[int]], true_set: Set[int], del_set: Set[int]) -> None:
    n = len(neighbours)
    belong = [set(neighbours[v]) for v in range(n)]

    for v in range(n):
        if v in del_set:
            continue

        Nv = neighbours[v]

        N1 = []
        for u in Nv:
            dominated = False
            for t in neighbours[u]:
                if t != v and (t not in belong[v]):
                    dominated = True
                    break
            if dominated:
                N1.append(u)
        N1_set = set(N1)

        N2 = []
        for u in Nv:
            if u in N1_set:
                continue
            for t in neighbours[u]:
                if t in N1_set:
                    N2.append(u)
                    break
        N2_set = set(N2)

        if len(Nv) == (len(N1) + len(N2)):
            continue
        
        N3 = [u for u in Nv if (u not in N1_set and u not in N2_set)]

        if N3 and (v not in del_set):
            true_set.add(v)
            for u in N3:
                del_set.add(u)
            for u in N2:
                del_set.add(u)

    del_set.difference_update(true_set)

def _rule2(
    neighbours: List[List[int]],
    at_least_pairs: Set[Tuple[int, int]],
    true_set: Set[int],
    del_set: Set[int],
) -> None:
    n = len(neighbours)
    belong = [set(neighbours[v]) for v in range(n)]

    for v in range(n):
        if v in del_set:
            continue

        within3 = _bfs_within_3(neighbours, v)
        for w in sorted(x for x in within3 if x > v):
            if w in del_set:
                continue

            N_vw = _reduction_merge(neighbours, v, w)

            N1 = []
            N_vw_set = set(N_vw)
            for u in N_vw:
                dominated = False
                for t in neighbours[u]:
                    if t != v and t != w and (t not in N_vw_set):
                        dominated = True
                        break
                if dominated:
                    N1.append(u)
            N1_set = set(N1)

            N2 = []
            for u in N_vw:
                if u in N1_set:
                    continue
                for t in neighbours[u]:
                    if t in N1_set:
                        N2.append(u)
                        break
            N2_set = set(N2)

            if len(N_vw) == (len(N1) + len(N2)):
                continue

            N3 = [u for u in N_vw if (u not in N1_set and u not in N2_set)]

            if len(N3) <= 1:
                continue

            is_dominated = False
            for u in N3:
                if _dominates(belong, u, N3):
                    is_dominated = True
                    break
            if not is_dominated:
                for u in N2:
                    if _dominates(belong, u, N3):
                        is_dominated = True
                        break
            if is_dominated:
                continue

            dominate_v = _dominates(belong, v, N3)
            dominate_w = _dominates(belong, w, N3)

            if dominate_v or dominate_w:
                if dominate_v and dominate_w:
                    if (v not in true_set) and (w not in true_set):
                        a, b = (v, w) if v < w else (w, v)
                        at_least_pairs.add((a, b))
                    for u in N3:
                        del_set.add(u)
                    for u in N2:
                        if (u in belong[v]) and (u in belong[w]):
                            del_set.add(u)
                elif dominate_v and (not dominate_w):
                    true_set.add(v)
                    for u in N3:
                        del_set.add(u)
                    for u in N2:
                        if u in belong[v]:
                            del_set.add(u)
                elif dominate_w and (not dominate_v):
                    true_set.add(w)
                    for u in N3:
                        del_set.add(u)
                    for u in N2:
                        if u in belong[w]:
                            del_set.add(u)
            else:
                true_set.add(v)
                true_set.add(w)
                for u in N3:
                    if u not in true_set:
                        del_set.add(u)
                for u in N2:
                    if u not in true_set:
                        del_set.add(u)

    del_set.difference_update(true_set)

def _handle_white_nodes(
    neighbours: List[List[int]],
    is_white: List[bool],
    new_neigh: List[List[int]],
    true_set: Set[int],
    del_set: Set[int],
) -> None:
    n = len(neighbours)

    for i in range(n):
        if not is_white[i]:
            continue
        new_neigh[i] = [v for v in new_neigh[i] if not is_white[v]]

    for i in range(n):
        if not is_white[i]:
            continue
        if i in true_set or i in del_set:
            continue
        deg = len(new_neigh[i])
        if deg == 0:
            del_set.add(i)
        elif deg == 1:
            node = new_neigh[i][0]
            if i in new_neigh[node]:
                new_neigh[node].remove(i)
            del_set.add(i)
            new_neigh[i].clear()

def _addition_rule(
    neighbours: List[List[int]],
    at_least_pairs: Set[Tuple[int, int]],
    true_set: Set[int],
    del_set: Set[int],
) -> None:
    n = len(neighbours)
    is_white = [False] * n
    new_neigh = [[] for _ in range(n)]

    for i in range(n):
        if i in true_set:
            for node in neighbours[i]:
                is_white[node] = True

        if (i not in true_set) and (i not in del_set):
            new_neigh[i] = [node for node in neighbours[i] if (node not in true_set and node not in del_set)]

    _handle_white_nodes(neighbours, is_white, new_neigh, true_set, del_set)

    find = False
    for i in range(n):
        if is_white[i]:
            continue
        if i in true_set or i in del_set:
            continue
        if len(new_neigh[i]) == 1:
            node = new_neigh[i][0]
            true_set.add(node)
            del_set.add(i)

            for neigh in list(new_neigh[node]):
                if neigh != i:
                    find = True
                is_white[neigh] = True
                if node in new_neigh[neigh]:
                    new_neigh[neigh].remove(node)
            new_neigh[node].clear()

    if find:
        _handle_white_nodes(neighbours, is_white, new_neigh, true_set, del_set)

    for i in range(n):
        if is_white[i]:
            continue
        if i in true_set or i in del_set:
            continue
        if len(new_neigh[i]) == 0:
            true_set.add(i)

    at_least_pairs.difference_update({(a, b) for (a, b) in at_least_pairs if (a in true_set or b in true_set)})

    del_set.difference_update(true_set)

def compute_reduction(dist_matrix: List[List[float]], n: int, radius: float, csv_path: str = "reduction.csv"):
    neighbours = _build_neighbours(dist_matrix, n, radius)

    Nc, Nd = set(), set()
    at_least_pairs: Set[Tuple[int, int]] = set()

    n_total = n

    t0_Nc, t0_Nd = len(Nc), len(Nd)
    _rule1(neighbours, Nc, Nd)
    t1_Nc, t1_Nd = len(Nc), len(Nd)
    
    r1_fixed = t1_Nc - t0_Nc
    r1_removed = t1_Nd - t0_Nd
    
    _rule2(neighbours, at_least_pairs, Nc, Nd)
    t2_Nc, t2_Nd = len(Nc), len(Nd)
    
    r2_fixed = t2_Nc - t1_Nc
    r2_removed = t2_Nd - t1_Nd
    r2_pairs = len(at_least_pairs)

    _addition_rule(neighbours, at_least_pairs, Nc, Nd)
    t3_Nc, t3_Nd = len(Nc), len(Nd)
        
    add_fixed = t3_Nc - t2_Nc
    add_removed = t3_Nd - t2_Nd

    total_fixed = len(Nc)
    total_removed = len(Nd)
    remaining = n_total - total_fixed - total_removed
    reduction_rate = ((total_fixed + total_removed) / n_total) * 100 if n_total > 0 else 0

    print(f"[RADIUS {radius}] Rem: {remaining} | Rate: {reduction_rate:.2f}% (R1:-{r1_removed}, R2:-{r2_removed}, Add:-{add_removed})")

    if csv_path:
        file_exists = os.path.isfile(csv_path)
        
        # Định nghĩa các cột
        fieldnames = [
            'radius', 
            'total_nodes', 
            'remaining_nodes',
            'reduction_rate_percent',
            'total_fixed', 
            'total_removed',
            'r1_fixed', 'r1_removed',
            'r2_fixed', 'r2_removed', 'r2_pairs',
            'add_fixed', 'add_removed'
        ]
        
        # Dữ liệu của dòng hiện tại
        row_data = {
            'radius': radius,
            'total_nodes': n_total,
            'remaining_nodes': remaining,
            'reduction_rate_percent': round(reduction_rate, 2),
            'total_fixed': total_fixed,
            'total_removed': total_removed,
            'r1_fixed': r1_fixed,
            'r1_removed': r1_removed,
            'r2_fixed': r2_fixed,
            'r2_removed': r2_removed,
            'r2_pairs': r2_pairs,
            'add_fixed': add_fixed,
            'add_removed': add_removed
        }

        try:
            with open(csv_path, mode='a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                
                if not file_exists:
                    writer.writeheader()
                
                writer.writerow(row_data)
        except Exception as e:
            print(f"[ERROR] Could not write to CSV: {e}")

    enabled_centers = set(range(n))
    demands = list(range(n))

    return Nc, Nd, enabled_centers, demands, at_least_pairs
