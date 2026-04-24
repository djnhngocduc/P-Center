import math 

def remaining_iters_from_bounds(lo: int, hi: int) -> int:
    width = hi - lo + 1
    if width <= 1:
        return 0
    return int(math.ceil(math.log2(width)))


def solver_wall_cost(status: str, wall_time: float, time_limit: float) -> float:
    if status in ("sat", "unsat"):
        return float(wall_time)
    if time_limit and time_limit > 0:
        return float(time_limit)
    return float(wall_time)