"""Geometry and fit checks for one rebuild. Read-only."""

import bmesh
import numpy as np

from . import nearest
from .build import inside
from .sampling import mm_per_local

BUDGETS = {'NUT': 150, 'FITTING': 400, 'ELBOW': 800, 'HOUSING': 2000}


def topology(ob):
    bm = bmesh.new()
    bm.from_mesh(ob.data)
    V, E, F = len(bm.verts), len(bm.edges), len(bm.faces)
    q = sum(1 for f in bm.faces if len(f.verts) == 4)
    t = sum(1 for f in bm.faces if len(f.verts) == 3)
    nonman = sum(1 for e in bm.edges if len(e.link_faces) != 2)
    zero = sum(1 for f in bm.faces if f.calc_area() <= 1e-14)
    tris = sum(len(f.verts) - 2 for f in bm.faces)

    seen, comps = set(), 0
    for v in bm.verts:
        if v in seen:
            continue
        comps += 1
        stack = [v]
        seen.add(v)
        while stack:
            w = stack.pop()
            for e in w.link_edges:
                o = e.other_vert(w)
                if o not in seen:
                    seen.add(o)
                    stack.append(o)
    bm.free()
    F1 = max(F, 1)
    return dict(verts=V, faces=F, tris=tris, quad_pct=100.0 * q / F1,
                tri_pct=100.0 * t / F1, ngon_pct=100.0 * (F - q - t) / F1,
                non_manifold=nonman, zero_area=zero, components=comps,
                euler=V - E + F, genus_ok=(comps == 1 and V - E + F == 2))


def _stats(d):
    if len(d) == 0:
        return dict(mean=0.0, p95=0.0, p99=0.0, max=0.0)
    return dict(mean=float(d.mean()), p95=float(np.percentile(d, 95)),
                p99=float(np.percentile(d, 99)), max=float(d.max()))


def _bbox(P):
    return P.max(axis=0) - P.min(axis=0)


def _tris(ob, depsgraph):
    """Evaluated vertices (local) and the triangles that are drawn. No writes."""
    ev = ob.evaluated_get(depsgraph)
    me = ev.to_mesh()
    try:
        V = np.empty(len(me.vertices) * 3)
        me.vertices.foreach_get("co", V)
        T = np.empty(len(me.loop_triangles) * 3, np.int64)
        me.loop_triangles.foreach_get("vertices", T)
    finally:
        ev.to_mesh_clear()
    return V.reshape(-1, 3), T.reshape(-1, 3)


def fit(ref, rb, depsgraph, scene, max_points=20000, on_surface_mm=0.01, exclude_inside=True):
    """Surface deviation BOTH ways with exact nearest-point distances.

    Ray-casting lies on hollow, open CAD shells, so this uses nearest points.
    Reference points INSIDE the rebuild are excluded (unless exclude_inside
    is off): for a closed-shell rebuild they are the bores, hole walls and
    interior deliberately dropped.
    """
    mm = mm_per_local(ref, scene)
    Vr, Tr = _tris(ref, depsgraph)
    ref_surf = nearest.Surface(Vr, Tr, exact_above=on_surface_mm / mm)

    # the rebuild in the reference's space, measured on the triangles that
    # are drawn and exported, not on the polygons: a curved quad split the
    # other way is a different surface
    Vb, Tb = _tris(rb, depsgraph)
    M = np.array(ref.matrix_world.inverted() @ rb.matrix_world)
    Pb = Vb @ M[:3, :3].T + M[:3, 3]
    rb_surf = nearest.Surface(Pb, Tb, exact_above=on_surface_mm / mm)
    # edge midpoints and face centres too: a fill over a hole shows only there
    X = Pb[Tb]
    Pb_extra = np.vstack([X.mean(1), (X + np.roll(X, -1, axis=1)).reshape(-1, 3) / 2.0])

    rng = np.random.default_rng(0)
    sample = Vr if len(Vr) <= max_points else Vr[rng.choice(len(Vr), max_points, False)]

    d_ref, excluded = [], 0
    for p, dist in zip(sample, rb_surf(sample)):
        if not np.isfinite(dist):
            continue
        dmm = dist * mm
        if exclude_inside and dmm > on_surface_mm and inside(rb_surf.tree, p):
            excluded += 1
            continue
        d_ref.append(dmm)

    d_rb = ref_surf(np.vstack([Pb, Pb_extra]))
    d_rb = d_rb[np.isfinite(d_rb)] * mm

    dims_ref = _bbox(Vr) * mm
    dims_rb = _bbox(Pb) * mm
    mdiff = max(abs(ref.matrix_world[i][j] - rb.matrix_world[i][j])
                for i in range(4) for j in range(4))
    return dict(ref_to_rb=_stats(np.array(d_ref)), rb_to_ref=_stats(np.array(d_rb)),
                excluded_inside=excluded, dims_ref=dims_ref.tolist(),
                dims_rb=dims_rb.tolist(),
                dims_err=float(np.abs(dims_ref - dims_rb).max()),
                matrix_diff=float(mdiff), ref_verts=len(Vr))


def budget_of(rb):
    return int(rb.get("xr_budget", 0)) or 2000


def run(ref, rb, depsgraph, scene):
    topo = topology(rb)
    res = dict(topo=topo, budget=budget_of(rb))
    if ref is not None:
        res["fit"] = fit(ref, rb, depsgraph, scene)
        ref_faces = len(ref.data.polygons)
        res["reduction"] = ref_faces / max(topo["faces"], 1)
    return res


def store(rb, res):
    t = res["topo"]
    rb["xr_faces"] = t["faces"]
    if rb.get("xr_mode") == "CAD_XR":
        # A faithful copy, not a closed shell: judged by the rebuild's own
        # two-way measure, taken where it knew which CAD faces were hidden and
        # removed. A plain distance check here would count each one as an error.
        tol = float(rb.get("xr_tol_mm", 0.2))
        rb["xr_genus_ok"] = bool(float(rb.get("xr_dev_max_mm", 0.0)) <= tol + 1e-6)
        rb["xr_over_budget"] = False
        return
    rb["xr_genus_ok"] = bool(t["genus_ok"] and t["non_manifold"] == 0
                             and t["zero_area"] == 0)
    rb["xr_over_budget"] = t["faces"] > res["budget"]
    if "fit" in res:
        rb["xr_dev_max_mm"] = res["fit"]["ref_to_rb"]["max"]
        rb["xr_dev_mean_mm"] = res["fit"]["ref_to_rb"]["mean"]


def lines(name, res):
    t = res["topo"]
    out = [f"== {name}",
           f"  verts {t['verts']}  faces {t['faces']}  (budget {res['budget']}"
           f"{' OVER' if t['faces'] > res['budget'] else ''})  tris {t['tris']}",
           f"  quad {t['quad_pct']:.0f}%  tri {t['tri_pct']:.0f}%  ngon {t['ngon_pct']:.0f}%",
           f"  non-manifold {t['non_manifold']}  zero-area {t['zero_area']}  "
           f"components {t['components']}  V-E+F {t['euler']}  "
           f"genus0 {'YES' if t['genus_ok'] else 'NO'}"]
    if "reduction" in res:
        out.append(f"  reduction x{res['reduction']:.1f} vs reference")
    if "fit" in res:
        f = res["fit"]
        a, b = f["ref_to_rb"], f["rb_to_ref"]
        out += [f"  dims ref {_fmt(f['dims_ref'])}  rebuilt {_fmt(f['dims_rb'])}"
                f"  (max err {f['dims_err']:.3f} mm)",
                f"  matrix max diff {f['matrix_diff']:.2e}",
                f"  ref->rebuild mm: mean {a['mean']:.3f} p95 {a['p95']:.3f} "
                f"p99 {a['p99']:.3f} max {a['max']:.3f}  "
                f"(excluded {f['excluded_inside']} interior/bore points)",
                f"  rebuild->ref mm: mean {b['mean']:.3f} p95 {b['p95']:.3f} "
                f"p99 {b['p99']:.3f} max {b['max']:.3f}"]
    return out


def _fmt(v):
    return "x".join(f"{x:.2f}" for x in v)
