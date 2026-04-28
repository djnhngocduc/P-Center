import os
import sys
import tempfile
import threading
import traceback
import shutil
import atexit

from pysat.solvers import Solver

from solvers.backend import (
    EXTERNAL_SOLVERS, 
    run_external_solver, 
    solve_radius_cplex, 
    solve_radius_gurobi, 
    solve_radius_cpo,
)

_CANCEL_SHARED = None
_INST_SHARED = None
_LOCAL_TMPDIR = None


def _pool_initializer(cancel_proxy, inst):
    global _CANCEL_SHARED, _LOCAL_TMPDIR, _INST_SHARED

    _CANCEL_SHARED = cancel_proxy
    _INST_SHARED = inst

    base_tmp = "/dev/shm" if os.path.isdir("/dev/shm") else None
    _LOCAL_TMPDIR = tempfile.mkdtemp(prefix="pcsat_worker_", dir=base_tmp)

    print(
        f"[POOL] worker PID={os.getpid()} tmpdir={_LOCAL_TMPDIR}",
        file=sys.stderr,
        flush=True,
    )

    def _cleanup_tmp():
        try:
            if _LOCAL_TMPDIR and os.path.isdir(_LOCAL_TMPDIR):
                shutil.rmtree(_LOCAL_TMPDIR, ignore_errors=True)
        except Exception:
            pass

    atexit.register(_cleanup_tmp)


def _solve_radius_worker_proc(idx, encoding, solver_name, radius, time_limit):
    pid = os.getpid()
    print(
        f"[WORKER-START] pid={pid} idx={idx} R={radius} "
        f"enc={encoding} solver={solver_name} limit={time_limit}s",
        flush=True
    )

    try:
        global _INST_SHARED

        if solver_name in ("cplex_mip", "cplex_cp", "gurobi_mip"):
            data = _INST_SHARED._build_setcover_data(radius)
            print(
                f"[ENCODE-{solver_name.upper()}] idx={idx} R={radius} "
                f"candidates={0 if data is None else len(data.get('candidates', []))} "
                f"bound={None if data is None else data.get('p')} "
                f"rows={0 if data is None else len(data.get('cover_rows', []))}",
                flush=True
            )
            if solver_name == "cplex_mip":
                return solve_radius_cplex(
                    inst=_INST_SHARED,
                    idx=idx,
                    radius=radius,
                    time_limit=time_limit,
                    data=data,
                )
            if solver_name == "gurobi_mip":
                return solve_radius_gurobi(
                    inst=_INST_SHARED,
                    idx=idx,
                    radius=radius,
                    time_limit=time_limit,
                    data=data,
                )
            return solve_radius_cpo(
                inst=_INST_SHARED,
                idx=idx,
                radius=radius,
                time_limit=time_limit,
                data=data,
            )

        cnf, varmap = _INST_SHARED._encode_cnf(radius, encoding)

        print(
            f"[ENCODE] idx={idx} R={radius} |Nc|={len(varmap.get('Nc', []))} "
            f"|Nd|={len(varmap.get('Nd', []))} "
            f"candidates={len(varmap.get('candidates', []))} "
            f"bound={varmap.get('bound')} "
            f"clauses={0 if cnf is None else len(cnf.clauses)} "
            f"vars={0 if cnf is None else cnf.nv}",
            flush=True
        )

        if cnf is None:
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> UNSAT(by reduction)",
                flush=True
            )
            return idx, radius, "unsat", None, None, None

        if solver_name in EXTERNAL_SOLVERS:
            cancel_ev = None
            try:
                if _CANCEL_SHARED is not None:
                    cancel_ev = _CANCEL_SHARED.get(idx)
            except Exception:
                cancel_ev = None

            status, model = run_external_solver(
                solver_name=solver_name,
                cnf=cnf,
                time_limit=time_limit,
                cancel_ev=cancel_ev,
                tmpdir=_LOCAL_TMPDIR,
            )

            nvars = cnf.nv
            nclauses = len(cnf.clauses)

            if status in ("timeout", "error"):
                print(
                    f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                    f"-> {status.upper()}",
                    flush=True
                )
                return idx, radius, status, nvars, nclauses, None

            if status == "unsat":
                print(
                    f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                    f"-> UNSAT",
                    flush=True
                )
                return idx, radius, "unsat", nvars, nclauses, None

            model_set = set(model or [])
            y_vars = varmap.get("y", [])
            Nc = set(varmap.get("Nc", []))
            candidates = set(varmap.get("candidates", []))
            chosen = {
                j
                for j, v in enumerate(y_vars)
                if (v in model_set) and (j in candidates)
            }
            centers = sorted(chosen | Nc)
            print(
                f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                f"-> SAT centers={centers}",
                flush=True
            )
            return idx, radius, "sat", nvars, nclauses, centers

        with Solver(name=solver_name, bootstrap_with=cnf.clauses) as solver:
            cancel_ev = None
            try:
                if _CANCEL_SHARED is not None:
                    cancel_ev = _CANCEL_SHARED.get(idx)
            except Exception:
                cancel_ev = None

            watcher = None
            if cancel_ev is not None and hasattr(solver, "interrupt"):
                def _watch():
                    cancel_ev.wait()
                    try:
                        solver.interrupt()
                    except Exception:
                        pass
                watcher = threading.Thread(target=_watch, daemon=True)
                watcher.start()

            timer = None
            try:
                if time_limit and time_limit > 0 and hasattr(solver, "interrupt"):
                    timer = threading.Timer(time_limit, solver.interrupt)
                    timer.start()
                    try:
                        sat = solver.solve_limited(expect_interrupt=True)
                    except NotImplementedError:
                        sat = solver.solve()

                else:
                    sat = solver.solve()

                nvars = cnf.nv
                nclauses = len(cnf.clauses)

                if sat is None:
                    print(
                        f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                        f"-> TIMEOUT",
                        flush=True
                    )
                    return idx, radius, "timeout", nvars, nclauses, None
                if not sat:
                    print(
                        f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                        f"-> UNSAT",
                        flush=True
                    )
                    return idx, radius, "unsat", nvars, nclauses, None

                model = solver.get_model() or []
                if not model:
                    print(
                        f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                        f"-> TIMEOUT(no model)",
                        flush=True
                    )
                    return idx, radius, "timeout", nvars, nclauses, None

                model_set = set(model)
                y_vars = varmap.get("y", [])
                Nc = set(varmap.get("Nc", []))
                candidates = set(varmap.get("candidates", []))
                chosen = {j for j, v in enumerate(y_vars) if (v in model_set) and (j in candidates)}
                centers = sorted(chosen | Nc)
                print(
                    f"[WORKER-END] pid={pid} idx={idx} R={radius} "
                    f"-> SAT centers={centers}",
                    flush=True
                )
                return idx, radius, "sat", nvars, nclauses, centers
            finally:
                if timer:
                    timer.cancel()

    except Exception:
        tb = traceback.format_exc()
        print(
            f"[WORKER-END] pid={pid} idx={idx} R={radius}\n{tb}",
            file=sys.stderr,
            flush=True
        )
        return idx, radius, "error", None, None, None