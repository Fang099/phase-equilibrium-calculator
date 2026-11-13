"""
precompute_surfaces.py -- solve every 3D surface the website offers, once,
and save it as a small file that the calculator page loads instead of
computing it in the browser.

    python scripts/precompute_surfaces.py            # only what is missing or outdated
    python scripts/precompute_surfaces.py --force    # everything again
    python scripts/precompute_surfaces.py --edition=full

For every ordered pair of the edition's compounds, every offered model that
has database parameters for that pair, and both fixed variables (T for Pxy
surfaces, P for Txy surfaces), it runs web_api.surface() with exactly the
request the page would send with the database parameters (no slider moved):
15 diagrams of 31 points across the slider's range. Results go to
docs/surfaces/<a>~<b>~<model>~<T|P>.json, with an index the page reads.

Each file records a fingerprint of the Python modules and the database. A
file whose fingerprint no longer matches is recomputed, so after changing
the code or the data, running this script (or letting the GitHub workflow
.github/workflows/surfaces.yml run it) brings every surface up to date.
Surfaces with a moved slider are still computed in the browser, and the
browser remembers them for the next visit.
"""

import hashlib
import json
import math
import os
import sys
import time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import web_api                                                     # noqa: E402

OUT = os.path.join(ROOT, "docs", "surfaces")
LAYERS, LAYER_PTS = 15, 31                    # must match docs/assets/app.js
MODULES = ["input_checker.py", "thermo_models.py", "equilibrium_solver.py", "numerics.py",
           "compound_db.py", "data/compounds.json"]
# from web_api.py only what a surface depends on, so changes to the rest of the
# web interface (new operations, messages) do not make every surface outdated
WEB_API_PARTS = ["EOS_MODELS", "ACTIVITY_MODELS", "pair_defaults", "build_input", "_solve",
                 "_columns", "_clean", "surface"]
# columns each plotted quantity needs (the same names the page uses)
Y_COLUMNS = {"T": ["T"], "P": ["P"], "y1": ["y1"], "gamma": ["gamma1", "gamma2"],
             "lngamma": ["gamma1", "gamma2"], "K": ["K1", "K2"], "f": ["f1", "f2"],
             "phiL": ["phiL1", "phiL2"], "phiV": ["phiV1", "phiV2"], "alpha": ["alpha"],
             "GE": ["GE_RT"], "Gmix": ["Gmix_RT"], "stab": ["stability"]}
ALL_MODELS = ["PR", "SRK", "vdW", "ideal", "margules", "vanlaar", "nrtl"]


def edition():
    with open(os.path.join(ROOT, "site", "edition.json"), encoding="utf-8") as fh:
        cfg = json.load(fh)
    name = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--edition=")),
                cfg["active"])
    return name, cfg["editions"][name]


def fingerprint():
    import ast
    h = hashlib.sha256()
    for f in MODULES:
        with open(os.path.join(ROOT, f), "rb") as fh:
            h.update(fh.read())
    with open(os.path.join(ROOT, "web_api.py"), encoding="utf-8") as fh:
        src = fh.read()
    for node in ast.parse(src).body:
        names = ([node.name] if isinstance(node, (ast.FunctionDef, ast.ClassDef)) else
                 [t.id for t in getattr(node, "targets", []) if isinstance(t, ast.Name)])
        if any(n in WEB_API_PARTS for n in names):
            h.update(ast.get_source_segment(src, node).encode())
    h.update(f"{LAYERS}/{LAYER_PTS}".encode())
    return h.hexdigest()[:16]


def keep_columns(ed):
    cols = ["x1", "y1", "T", "P"]
    for y in ed.get("y_axes", Y_COLUMNS):
        cols += [c for c in Y_COLUMNS.get(y, []) if c not in cols]
    return cols


def layer_values(pd, hold):
    """The slider's range, as the page spaces it (linear in T, logarithmic in P)."""
    if hold == "T":
        lo, hi = pd["ranges"]["T"]
        return [round(lo + (hi - lo) * k / (LAYERS - 1), 2) for k in range(LAYERS)]
    lo, hi = math.log10(pd["ranges"]["P"][0]), math.log10(pd["ranges"]["P"][1])
    return [float(f"{10 ** (lo + (hi - lo) * k / (LAYERS - 1)):.5g}") for k in range(LAYERS)]


def default_params(pd, model):
    """What the page sends when no slider has been moved."""
    if model == "PR":
        return {"kij": pd["kij"]["value"]}
    if model in ("SRK", "vdW"):
        return {"kij": 0}
    if model == "vanlaar":
        return {"A12": pd["vanlaar"]["A12"], "A21": pd["vanlaar"]["A21"], "tdep": 1}
    if model == "margules":
        return {"A12": pd["A"]["A12"], "A21": pd["A"]["A21"]}
    if model == "nrtl":
        n = pd["nrtl"]
        return {"dg12": n["dg12"], "dg21": n["dg21"], "alpha": n["alpha"]}
    return {}


def jobs(ed):
    cat = web_api.catalog()
    ids = [c["id"] for c in cat["compounds"] if c["id"] in ed.get("compounds", [c["id"]])]
    models = [m for m in ALL_MODELS if m in ed.get("models", ALL_MODELS)]
    keep = keep_columns(ed)
    out = []
    for a in ids:
        for b in ids:
            if a == b:
                continue
            pd = web_api.pair_defaults(a, b)
            for m in models:
                if (m == "vanlaar" and not pd["vanlaar"]) or (m == "nrtl" and not pd["nrtl"]):
                    continue                      # the page does not offer these without data
                for hold in ("T", "P"):
                    req = {"components": [a, b], "model": m, "params": default_params(pd, m),
                           "hold": hold, "npts": LAYER_PTS, "values": layer_values(pd, hold),
                           "keep": keep}
                    out.append((f"{a}~{b}~{m}~{hold}", req))
    return out


def run_one(args):
    key, req, fp = args
    t0 = time.time()
    res = web_api.surface(req)
    with open(os.path.join(OUT, key + ".json"), "w", encoding="utf-8") as fh:
        json.dump({"version": fp, "request": req, "result": res}, fh, separators=(",", ":"))
    return key, res.get("ok", False), time.time() - t0


def main():
    name, ed = edition()
    os.makedirs(OUT, exist_ok=True)
    fp = fingerprint()
    todo, have = [], []
    for key, req in jobs(ed):
        path = os.path.join(OUT, key + ".json")
        fresh = False
        if os.path.exists(path) and "--force" not in sys.argv:
            try:
                with open(path, encoding="utf-8") as fh:
                    old = json.load(fh)
                fresh = old.get("version") == fp and old.get("request") == req
            except (OSError, ValueError):
                pass
        (have if fresh else todo).append((key, req, fp))
    print(f"{name} edition: {len(have) + len(todo)} surfaces, {len(have)} up to date, "
          f"{len(todo)} to compute")
    t0 = time.time()
    with Pool(max(1, os.cpu_count() or 1)) as pool:
        for k, (key, ok, dt) in enumerate(pool.imap_unordered(run_one, todo), 1):
            print(f"  [{k}/{len(todo)}] {key}: {'ok' if ok else 'FAILED'} ({dt:.1f} s)", flush=True)
    keys = sorted(k for k, _, _ in have + todo)
    wanted = {k + ".json" for k in keys} | {"index.json"}
    for f in os.listdir(OUT):                   # surfaces no longer offered
        if f.endswith(".json") and f not in wanted:
            os.remove(os.path.join(OUT, f))
    with open(os.path.join(OUT, "index.json"), "w", encoding="utf-8") as fh:
        json.dump({"version": fp, "edition": name, "layers": LAYERS, "points": LAYER_PTS,
                   "keys": keys}, fh, separators=(",", ":"))
    print(f"done in {time.time() - t0:.0f} s; {len(keys)} surfaces in docs/surfaces/")


if __name__ == "__main__":
    main()
