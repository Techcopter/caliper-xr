"""CAD -> XR: the shape-preserving pipeline for one object.

    load -> (hide what no one can see) -> recover CAD faces -> height-field
    pieces -> simplify the curve network -> rebuild every face lean ->
    assemble: n-gons on planes, quads on strips, CAD normals, materials.
"""

import math
import time

import bmesh
import bpy
import numpy as np

from . import build, cadmesh, nearest, param, retess, segment, visibility


def vertex_normals(cm, patch):
    """Smooth normal of each (vertex, CAD face): the source's own normals,
    area-weighted, averaged inside the face only."""
    np1 = int(patch.max()) + 1
    v = cm.T.reshape(-1)
    p = np.repeat(patch, 3)
    w = np.repeat(cm.A, 3)
    n = cm.CN.reshape(-1, 3) * w[:, None]
    key = v.astype(np.int64) * np1 + p
    uk, inv = np.unique(key, return_inverse=True)
    acc = np.zeros((len(uk), 3))
    np.add.at(acc, inv.reshape(-1), n)
    L = np.linalg.norm(acc, axis=1)
    L[L == 0] = 1.0
    acc /= L[:, None]
    return {(int(k // np1), int(k % np1)): tuple(acc[i]) for i, k in enumerate(uk)}


def vertex_uvs(cm, patch):
    """UV of each (vertex, CAD face). Seams are patch borders, so one UV per key."""
    if getattr(cm, "UV", None) is None:
        return None
    np1 = int(patch.max()) + 1
    key = cm.T.reshape(-1).astype(np.int64) * np1 + np.repeat(patch, 3)
    uk, first = np.unique(key, return_index=True)
    uv = cm.UV.reshape(-1, 2)[first]
    return {(int(k // np1), int(k % np1)): (float(a), float(b)) for k, (a, b) in zip(uk, uv)}


def left_out(cm, patch, ob, occ=None, shell=None):
    """CAD faces to leave out, as two masks per face: never seen by the
    camera test (`occ`), and not on the outer shell (`shell`)."""
    npatch = int(patch.max()) + 1
    hidden, inner = np.zeros(npatch, bool), np.zeros(npatch, bool)
    if occ is not None:
        method, data = occ
        if method == "gpu":
            vis = data.get(ob.name)
            if vis is not None:
                hidden = visibility.hidden_from_faces(patch, vis[cm.src])
        else:
            tree, diag = data
            hidden = visibility.hidden_patches(cm, patch, ob.matrix_world, tree, diag)
    if hidden.any():
        # A hidden face whose neighbours are all seen would leave a hole in a
        # visible surface: the camera test samples 128 directions, and a
        # narrow bore can hide from every one. Whole hidden bodies stay out.
        two = cm.ef[:, 0] >= 0
        a, b = patch[cm.ef[two, 0]], patch[cm.ef[two, 1]]
        a, b = a[a != b], b[a != b]
        nb_hidden, nb = np.zeros(npatch, bool), np.zeros(npatch, bool)
        np.logical_or.at(nb_hidden, a, hidden[b])
        np.logical_or.at(nb_hidden, b, hidden[a])
        nb[a] = nb[b] = True
        hidden &= ~(nb & ~nb_hidden)
    if shell is not None and shell.get(ob.name) is not None:
        inner = visibility.hidden_from_faces(patch, shell[ob.name][cm.src]) & ~hidden
    return hidden, inner


def run(ob, depsgraph, s, lpm, occ=None, progress=None, shell=None):
    """Returns (bmesh, loop_normals, info, stats); bmesh is None if hidden."""
    t0 = time.time()
    info = []
    cm = cadmesh.load(ob, depsgraph)
    tol = s.xr_tol_mm * lpm
    hard = segment.hard_edges(cm, s.xr_crease_deg)
    patch = segment.patches(cm, hard)
    n_src = len(cm.T)
    info.append(f"  source {n_src} triangles; {patch.max() + 1} CAD faces recovered")

    hidden, inner = left_out(cm, patch, ob, occ, shell)
    drop = hidden | inner
    if drop.all():
        info.append("  nothing of it is on the outside: nothing to export")
        return None, None, info, {"src_tris": len(cm.T), "hidden": True}
    if drop.any():
        keepf = ~drop[patch]
        info.append(f"  left out {int(drop.sum())} CAD faces ({int((~keepf).sum())} triangles): "
                    f"{int(hidden.sum())} hidden, {int(inner.sum())} inside the outer shell")
        cm = cadmesh.subset(cm, keepf)
        hard = segment.hard_edges(cm, s.xr_crease_deg)
        patch = segment.patches(cm, hard)

    patch = segment.extract_planes(cm, patch, hard, max(tol * 0.05, cm.diag * 1e-6), (2.0 * lpm) ** 2)
    flat = segment.planar_patches(cm, patch)
    UV, uv_good = param.unwrap(cm, patch, ~flat)
    sub, pieces = segment.make_pieces(cm, patch, hard, flat, uv_good)
    bnd = segment.boundary_edges(cm, sub, hard)
    vn = np.zeros_like(cm.V)
    np.add.at(vn, cm.T.reshape(-1), np.repeat(cm.N * cm.A[:, None], 3, axis=0))
    vn /= np.maximum(np.linalg.norm(vn, axis=1), 1e-30)[:, None]
    near = nearest.Surface(cm.V, cm.T, exact_above=tol)       # distance to the source
    keep = segment.simplify_chains(cm, sub, bnd, tol, hard, vn, surf=near)
    vnorm = vertex_normals(cm, patch)
    res, state = retess.rebuild(cm, sub, pieces, UV, bnd, keep, tol, patch, vnorm, progress,
                                near=near)
    # Global guarantee, both ways: every source point within tol of the
    # result, and every result face within tol of the source. First put back
    # dropped outline vertices near what is off; a piece still off after that
    # goes back to its exact source triangles, its neighbours rebuilt against it.
    for it in range(11):
        m = _measure(cm, res, near, tol)
        if m["worst"] <= tol or it == 10:
            break
        redo = _restore(cm, sub, bnd, keep, tol, m) if it < 4 else set()
        if not redo:
            redo = retess.force_source(cm, sub, bnd, keep, state, _over(cm, res, sub, tol, m))
        res, state = retess.rebuild(cm, sub, pieces, UV, bnd, keep, tol, patch, vnorm,
                                    progress, state=state, only=redo, near=near)
    # the passes above only needed to know what exceeds tol; the reported
    # figure is the exact largest distance each way
    if "to_res" in m:
        m["fwd_max"] = m["to_res"].exact_max(m["S"], m["fwd"])
        m["rev_max"] = near.exact_max(m["R"], m["d"])
        m["worst"] = max(m["fwd_max"], m["rev_max"])
    info.append(f"  largest deviation {m['worst'] / lpm:.3f} mm (tolerance {tol / lpm:.3f}): "
                f"CAD to XR {m['fwd_max'] / lpm:.3f}, XR to CAD {m['rev_max'] / lpm:.3f}")
    info_worst = m["worst"] / lpm
    planar = np.array([pc.mode == 'plane' for pc in pieces], bool)
    nuv = sum(1 for pc in pieces if pc.mode == 'uv')
    nproj = sum(1 for pc in pieces if pc.mode == 'proj')
    info.append(f"  {int(flat.sum())} flat faces, {nuv} curved faces unwrapped whole, "
                f"{nproj} projection pieces")
    vuv = vertex_uvs(cm, patch)
    bm, loopn = assemble(cm, res, patch, sub, planar, vnorm, s.xr_crease_deg, vuv)
    stats = {"src_tris": n_src, "worst_mm": info_worst, "faces": len(set(patch.tolist()))}
    info.append(f"  rebuilt {len(pieces)} pieces"
                + (f" ({res.fallback} kept as in the source)" if res.fallback else "")
                + f" in {time.time() - t0:.1f}s")
    return bm, loopn, info, stats


def _verts(cm, res):
    return np.vstack([cm.V, np.array(res.pts, float).reshape(-1, 3)])


def _measure(cm, res, near, tol):
    """Deviation both ways at sample points: every source vertex and face
    centre to the result, every result face centre and edge midpoint to the
    source (a fill over a hole shows only in that direction)."""
    S = np.vstack([cm.V[np.unique(cm.T)], cm.V[cm.T].mean(1)])
    T = np.asarray(res.tris, np.int64).reshape(-1, 3)
    if not len(T):
        return {"S": S, "fwd": np.full(len(S), np.inf), "rev": np.zeros(0),
                "rp": np.zeros((0, 3)), "fwd_max": float("inf"), "rev_max": 0.0,
                "worst": float("inf")}
    VA = _verts(cm, res)
    to_res = nearest.Surface(VA, T, exact_above=tol)
    fwd = to_res(S)
    X = VA[T]
    R = np.stack([X.mean(1), (X[:, 0] + X[:, 1]) / 2, (X[:, 1] + X[:, 2]) / 2,
                  (X[:, 2] + X[:, 0]) / 2], 1)
    d = near(R.reshape(-1, 3)).reshape(-1, 4)
    k = d.argmax(1)
    i = np.arange(len(d))
    fmax, rmax = float(fwd.max()), float(d.max())
    return {"S": S, "fwd": fwd, "rev": d[i, k], "rp": R[i, k], "worst": max(fmax, rmax),
            "fwd_max": fmax, "rev_max": rmax,
            "to_res": to_res, "R": R.reshape(-1, 3), "d": d.reshape(-1)}


def _over(cm, res, sub, tol, m):
    """Pieces with a sample farther than tol, either way."""
    nvs = len(m["S"]) - len(cm.T)
    bad = np.flatnonzero(m["fwd"] > tol)
    over = set(int(x) for x in np.unique(sub[bad[bad >= nvs] - nvs]))
    vids = np.unique(cm.T)[bad[bad < nvs]]
    if len(vids):
        fs = np.flatnonzero(np.isin(cm.T, vids).any(axis=1))
        over.update(int(x) for x in np.unique(sub[fs]))
    over.update(int(x) for x in np.unique(np.asarray(res.tsub)[m["rev"] > tol]))
    return over


def _restore(cm, sub, bnd, keep, tol, m):
    """Put back the dropped outline vertex nearest each spot that is off;
    returns the pieces to rebuild (none when there is nothing to put back)."""
    from mathutils import Vector
    from mathutils.kdtree import KDTree
    bad = np.vstack([m["S"][m["fwd"] > tol], m["rp"][m["rev"] > tol]])
    bv = np.unique(cm.E[bnd].ravel())
    dropped = bv[~keep[bv]]
    if not len(bad) or not len(dropped):
        return set()
    kd = KDTree(len(dropped))
    for i, v in enumerate(dropped):
        kd.insert(Vector(cm.V[v]), i)
    kd.balance()
    newly = np.array(sorted({int(dropped[kd.find(Vector(p))[1]]) for p in bad}))
    keep[newly] = True
    touched = np.flatnonzero(np.isin(cm.T, newly).any(axis=1))
    return set(int(x) for x in np.unique(sub[touched]))


def assemble(cm, res, patch, sub, planar, vnorm, crease_deg, vuv=None):
    nv = len(cm.V)
    nsub = len(planar)
    sub_patch = np.zeros(nsub, np.int64)
    sub_patch[sub] = patch
    sub_mat = np.zeros(nsub, np.int64)
    sub_mat[sub] = cm.mat
    flat_patch = np.zeros(int(patch.max()) + 1, bool)
    flat_patch[sub_patch[planar]] = True

    bm = bmesh.new()
    gl = bm.verts.layers.int.new("xr_gid")
    pl = bm.faces.layers.int.new("xr_patch")
    tris = np.array(res.tris, np.int64).reshape(-1, 3)
    vmap = {}
    for g in np.unique(tris):
        g = int(g)
        co = cm.V[g] if g < nv else res.pts[g - nv]
        v = bm.verts.new(tuple(co))
        v[gl] = g
        vmap[g] = v
    for (a, b, c), sp in zip(tris, res.tsub):
        if a == b or b == c or a == c:
            continue
        try:
            f = bm.faces.new((vmap[int(a)], vmap[int(b)], vmap[int(c)]))
        except ValueError:
            continue                                    # duplicate face
        f[pl] = int(sub_patch[sp])
        f.material_index = int(sub_mat[sp])
        f.smooth = True
    bm.normal_update()

    # Seams on CAD-face borders keep quads and n-gons inside one face.
    for e in bm.edges:
        lf = e.link_faces
        e.seam = len(lf) != 2 or lf[0][pl] != lf[1][pl]
    curved = [f for f in bm.faces if not flat_patch[f[pl]]]
    if curved:
        bmesh.ops.join_triangles(bm, faces=curved, cmp_seam=True, cmp_materials=True,
                                 angle_face_threshold=math.radians(10.0),
                                 angle_shape_threshold=math.radians(40.0))
        # A quad is drawn and exported split along loops 0-2. Make that the
        # diagonal the rebuild chose and measured: the other one can sit
        # farther from the CAD than the tolerance on a curved face.
        own = {(min(a, b), max(a, b)) for t in res.tris for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
        for f in [f for f in bm.faces if len(f.verts) == 4 and not flat_patch[f[pl]]]:
            vs = list(f.verts)
            if (min(vs[0][gl], vs[2][gl]), max(vs[0][gl], vs[2][gl])) in own:
                continue
            mi, sm, pv = f.material_index, f.smooth, f[pl]
            bm.faces.remove(f)
            nf = bm.faces.new(vs[1:] + vs[:1])
            nf.material_index, nf.smooth, nf[pl] = mi, sm, pv
        bm.normal_update()
    inner = [e for e in bm.edges
             if not e.seam and len(e.link_faces) == 2 and flat_patch[e.link_faces[0][pl]]]
    if inner:
        bmesh.ops.dissolve_limit(bm, angle_limit=math.radians(1.0), verts=[], edges=inner,
                                 delimit={'SEAM', 'MATERIAL'})

    lim = math.cos(math.radians(crease_deg))
    for e in bm.edges:
        lf = e.link_faces
        e.smooth = not (len(lf) == 2 and e.seam and lf[0].normal.dot(lf[1].normal) < lim)
        e.seam = False
    uvl = bm.loops.layers.uv.new("UVMap") if vuv is not None else None
    loopn = []
    for f in bm.faces:
        p = f[pl]
        if uvl is not None:
            for lp in f.loops:
                g = lp.vert[gl]
                if g < nv:
                    lp[uvl].uv = vuv.get((g, p), (0.0, 0.0))
                else:
                    u = res.puv[g - nv]
                    lp[uvl].uv = (float(u[0]), float(u[1])) if u is not None else (0.0, 0.0)
        fn = tuple(f.normal)
        for lp in f.loops:
            g = lp.vert[gl]
            loopn.append(vnorm.get((g, p), fn) if g < nv else tuple(res.pnrm[g - nv]))
    return bm, loopn


def make_object(bm, loopn, name, ref, col, keep_materials, meta):
    """Write the result as a new object in REBUILT; the source is untouched."""
    old = bpy.data.objects.get(name)
    if old is not None and build.is_rebuild(old):
        me_old = old.data
        bpy.data.objects.remove(old, do_unlink=True)
        if me_old is not None and me_old.users == 0:
            bpy.data.meshes.remove(me_old)
    me = bpy.data.meshes.new(name)
    bm.to_mesh(me)
    bm.free()
    for a in ("xr_gid", "xr_patch"):
        if me.attributes.get(a) is not None:
            me.attributes.remove(me.attributes[a])
    if keep_materials:
        # the effective material of each slot: parts from CAD often link
        # materials to the object, not the mesh, and the colour lives there
        for slot in ref.material_slots:
            me.materials.append(slot.material)
    if len(loopn) == len(me.loops):
        me.normals_split_custom_set(loopn)
    ob = bpy.data.objects.new(name, me)
    col.objects.link(ob)
    ob.rotation_mode = ref.rotation_mode
    ob.matrix_world = ref.matrix_world.copy()
    for k, v in meta.items():
        ob[k] = v
    return ob
