# check_symmetry.py
from collections import defaultdict
from encoder import PCenterSAT

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

def find_coverage_orbits(
    inst: PCenterSAT,
    radius: float,
):
    """
    Phát hiện orbit logic theo tập phủ,
    dùng CHÍNH reduction pipeline của encoder.py
    """

    # === DÙNG LẠI reduction gốc ===
    Nc, Nd, enabled_centers, demands, _ = inst.compute_reduction(radius)

    candidates = sorted((enabled_centers - Nc) - Nd)
    dist = inst.dist
    rad = radius + 1e-12

    # -----------------------
    # Active demands (y hệt sb.py)
    # -----------------------
    active_demands = []
    if Nc:
        for u in demands:
            for c in Nc:
                if dist[c][u] <= rad:
                    break
            else:
                active_demands.append(u)
    else:
        active_demands = list(demands)

    # -----------------------
    # Build coverage mask
    # -----------------------
    mask2cands = defaultdict(list)

    for c in candidates:
        mask = 0
        row = dist[c]
        for bit, u in enumerate(active_demands):
            if row[u] <= rad:
                mask |= (1 << bit)
        mask2cands[mask].append(c)

    # -----------------------
    # Chỉ giữ orbit size >= 2
    # -----------------------
    orbits = {m: g for m, g in mask2cands.items() if len(g) >= 2}

    return {
        "radius": radius,
        "Nc": len(Nc),
        "Nd": len(Nd),
        "candidates": len(candidates),
        "active_demands": len(active_demands),
        "num_orbits": len(orbits),
        "nodes_in_orbits": sum(len(g) for g in orbits.values()),
        "max_orbit_size": max((len(g) for g in orbits.values()), default=0),
    }

def analyze_tsplib(tsplib_path, p):
    coords = load_tsplib_coords(tsplib_path)
    inst = PCenterSAT.from_coordinates(coords, p)

    print(f"\n=== Instance {tsplib_path} | n={inst.n} p={p} ===")

    # duyệt 1 số bán kính đại diện
    radii = inst.radii[::max(1, len(inst.radii)//2000)]
    print(len(radii))
    for R in radii:
        info = find_coverage_orbits(inst, R)

        if info["num_orbits"] > 0:
            print(
                f"R={R:.2f} | "
                f"active={info['active_demands']} | "
                f"orbits={info['num_orbits']} | "
                f"nodes={info['nodes_in_orbits']} | "
                f"max_orbit={info['max_orbit_size']}"
            )


if __name__ == "__main__":
    analyze_tsplib("dataset/TSPLIB/u1060.tsp", p=150)