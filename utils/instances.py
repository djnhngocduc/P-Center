import json
import math
import time
import heapq
from typing import List, Tuple

from encoder import PCenterSAT

INF = 10 ** 12


def load_instance_data(inst_file: str):
    with open(inst_file, "r", encoding="utf-8", errors="ignore") as f:
        return json.load(f)


def load_instance(inst_desc, use_seed_radius: bool = False):
    if "tsplib" in inst_desc:
        coords = load_tsplib_coords(inst_desc["tsplib"])
        inst = PCenterSAT.from_coordinates(coords, inst_desc["p"])

        inst.radii = sorted(set(inst.radii), reverse=True)

        if use_seed_radius:
            sr = math.floor(100.0 * float(inst_desc["seed_radius"]) + 0.5) / 100.0
            seed_idx = next((i for i, r in enumerate(inst.radii) if r <= sr), None)
        else:
            seed_idx = 0

    elif "orlib" in inst_desc:
        t0 = time.time()
        n_graph, p, edges = load_orlib_edge_list(inst_desc["orlib"])
        print(f"[LOAD] {inst_desc['name']}: n={n_graph}, m={len(edges)}, p={p}", flush=True)

        t1 = time.time()
        D = compute_apsp_dijkstra(n_graph, edges)
        print(
            f"[APSP] {inst_desc['name']}: done in {time.time() - t1:.2f}s "
            f"(total load {time.time() - t0:.2f}s)",
            flush=True
        )
        inst = PCenterSAT.from_distance_matrix(D, p)

        radii_set = set()
        for i in range(n_graph):
            for j in range(i + 1, n_graph):
                dij = D[i][j]
                if dij != INF and dij > 0:
                    radii_set.add(dij)

        inst.radii = sorted(radii_set, reverse=True)

        if use_seed_radius:
            sr = int(inst_desc["seed_radius"])
            seed_idx = next((i for i, r in enumerate(inst.radii) if r <= sr), None)
        else:
            seed_idx = 0
    else:
        raise ValueError("Instance description must contain 'tsplib' or 'orlib'")

    if inst.radii:
        seed_R = inst.radii[seed_idx]
        print(
            f"[SEED] {inst_desc['name']}: "
            f"seed_idx={seed_idx}, seed_radius={seed_R}, "
            f"radii_len={len(inst.radii)}, p={inst.p}, n={inst.n}",
            flush=True
        )
    else:
        print(f"[SEED] {inst_desc['name']}: radii empty?!", flush=True)

    return inst, seed_idx


def load_tsplib_coords(path: str):
    coords = []
    in_section = False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            u = s.upper()
            if u.startswith("NODE_COORD_SECTION"):
                in_section = True
                continue
            if u.startswith("EOF"):
                break
            if not in_section:
                continue
            parts = s.split()
            if len(parts) >= 3:
                x = float(parts[1])
                y = float(parts[2])
                coords.append((x, y))
    return coords


def load_orlib_edge_list(path: str):
    edges = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().split()

        n = int(header[0])
        m = int(header[1])
        p = int(header[2])

        for _ in range(m):
            line = f.readline()
            if not line:
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            u, v, w = int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2])
            edges.append((u, v, w))

    return n, p, edges


def compute_apsp_dijkstra(n: int, edges: List[Tuple[int, int, float]]) -> List[List[int]]:
    c = [[INF] * n for _ in range(n)]
    for i in range(n):
        c[i][i] = 0

    for u, v, w in edges:
        ww = int(w)
        if u == v:
            c[u][v] = 0
        else:
            c[u][v] = ww
            c[v][u] = ww

    adj = [[] for _ in range(n)]
    for u in range(n):
        row = c[u]
        for v in range(n):
            w = row[v]
            if u != v and w < INF:
                adj[u].append((v, w))

    D = [[INF] * n for _ in range(n)]
    for s in range(n):
        dist = [INF] * n
        dist[s] = 0
        pq = [(0, s)]
        while pq:
            d, u = heapq.heappop(pq)
            if d > dist[u]:
                continue
            for v, w in adj[u]:
                nd = d + w
                if nd < dist[v]:
                    dist[v] = nd
                    heapq.heappush(pq, (nd, v))
        D[s] = [int(x) if x < INF else INF for x in dist]

    return D