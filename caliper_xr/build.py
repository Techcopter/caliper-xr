"""Build the lean outer shell and clean it up. Writes only NEW datablocks."""

import math

import bmesh
import bpy
import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree

from . import geom

REBUILT = "REBUILT"


# --------------------------------------------------------------------------
# Collection + naming
# --------------------------------------------------------------------------

def rebuilt_collection(scene):
    col = bpy.data.collections.get(REBUILT)
    if col is None:
        col = bpy.data.collections.new(REBUILT)
    if col.name not in scene.collection.children:
        scene.collection.children.link(col)
    # An excluded layer collection has no evaluated mesh, and verification
    # then fails with "has no evaluated mesh data".
    lc = next((c for c in bpy.context.view_layer.layer_collection.children
               if c.name == REBUILT), None)
    if lc is not None and lc.exclude:
        lc.exclude = False
    return col


def is_rebuild(ob):
    return any(c.name == REBUILT for c in ob.users_collection)


def rebuild_name(ref):
    return ref.name + "_REBUILT"


# --------------------------------------------------------------------------
# Loft
# --------------------------------------------------------------------------

def loft(grid, station_idx, angle_idx):
    """Rings at the kept stations, sharing one angle set; capped both ends."""
    C, a, u, v = grid["frame"]
    lpm = grid["lpm"]
    thetas = grid["thetas"][angle_idx]
    cs, sn = np.cos(thetas), np.sin(thetas)
    pole = 1e-4                                   # mm

    bm = bmesh.new()
    rings = []
    last = len(station_idx) - 1
    for n, k in enumerate(station_idx):
        z = grid["z"][k] * lpm
        rr = grid["rows"][k][angle_idx]
        base = C + a * z
        if rr.max() <= pole and n in (0, last):
            rings.append([bm.verts.new(tuple(base))])
            continue
        rr = np.maximum(rr, pole) * lpm
        pts = base + np.outer(rr * cs, u) + np.outer(rr * sn, v)
        rings.append([bm.verts.new(tuple(p)) for p in pts])

    for r0, r1 in zip(rings, rings[1:]):
        if len(r0) == 1 and len(r1) == 1:
            continue
        if len(r0) == 1:
            m = len(r1)
            for j in range(m):
                bm.faces.new((r0[0], r1[j], r1[(j + 1) % m]))
        elif len(r1) == 1:
            m = len(r0)
            for j in range(m):
                bm.faces.new((r0[j], r1[0], r0[(j + 1) % m]))
        else:
            m = len(r0)
            for j in range(m):
                bm.faces.new((r0[j], r1[j], r1[(j + 1) % m], r0[(j + 1) % m]))

    if len(rings[0]) > 2:
        bm.faces.new(list(reversed(rings[0])))
    if len(rings[-1]) > 2:
        bm.faces.new(rings[-1])
    return bm


def hull(V):
    """Convex hull fallback: genus-0 by construction, coarse by nature."""
    bm = bmesh.new()
    for p in V:
        bm.verts.new(tuple(p))
    bmesh.ops.convex_hull(bm, input=bm.verts, use_existing_faces=False)
    loose = [v for v in bm.verts if not v.link_faces]
    if loose:
        bmesh.ops.delete(bm, geom=loose, context='VERTS')
    return bm


# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------

def orient(bm):
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])
    if bm.calc_volume(signed=True) < 0:
        bmesh.ops.reverse_faces(bm, faces=bm.faces[:])


def _health(bm):
    return (len(bm.verts) - len(bm.edges) + len(bm.faces),
            sum(1 for v in bm.verts if not v.is_manifold),
            sum(1 for e in bm.edges if not e.is_manifold))


def _guarded(bm, step):
    """Run a cleanup step only if it keeps the shell as sound as it was.

    Welding at 0.005 mm can pinch two boolean slivers into a bow-tie vertex,
    turning a genus-0 union into V-E+F=0. The step is tried on a copy first
    and skipped if it adds non-manifold geometry or changes V-E+F.
    """
    before = _health(bm)
    trial = bm.copy()
    try:
        step(trial)
        after = _health(trial)
    finally:
        trial.free()
    if after[0] == before[0] and after[1] <= before[1] and after[2] <= before[2]:
        step(bm)


def cleanup(bm, lpm, sharp_deg=24.0):
    """Weld, dissolve coplanar/collinear, quad-up, drop loose, orient, shade.

    Dissolving coplanar faces is not decimation: it removes no silhouette,
    only edges that describe nothing.
    """
    def weld(b):
        bmesh.ops.remove_doubles(b, verts=b.verts[:], dist=0.005 * lpm)

    def degenerate(b):
        bmesh.ops.dissolve_degenerate(b, edges=b.edges[:], dist=0.0001 * lpm)

    def limit(b):
        bmesh.ops.dissolve_limit(b, angle_limit=math.radians(0.5),
                                 verts=b.verts[:], edges=b.edges[:],
                                 delimit={'NORMAL'})

    def quads(b):
        tris = [f for f in b.faces if len(f.verts) == 3]
        if tris:
            bmesh.ops.join_triangles(b, faces=tris, angle_face_threshold=0.01,
                                     angle_shape_threshold=3.14)

    for step in (weld, degenerate, limit, quads):
        _guarded(bm, step)
    wire = [e for e in bm.edges if not e.link_faces]
    if wire:
        bmesh.ops.delete(bm, geom=wire, context='EDGES')
    lone = [v for v in bm.verts if not v.link_edges]
    if lone:
        bmesh.ops.delete(bm, geom=lone, context='VERTS')
    orient(bm)
    lim = math.radians(sharp_deg)
    for e in bm.edges:
        e.smooth = not (len(e.link_faces) == 2 and e.calc_face_angle(0.0) > lim)
    for f in bm.faces:
        f.smooth = True


def to_object(bm, name, ref, col, meta):
    """Write a new mesh object; copy the reference transform exactly."""
    old = bpy.data.objects.get(name)
    if old is not None and is_rebuild(old):
        me_old = old.data
        bpy.data.objects.remove(old, do_unlink=True)
        if me_old is not None and me_old.users == 0:
            bpy.data.meshes.remove(me_old)
    me = bpy.data.meshes.new(name)
    bm.to_mesh(me)
    bm.free()
    me.materials.clear()
    ob = bpy.data.objects.new(name, me)
    col.objects.link(ob)
    if ref is not None:
        ob.rotation_mode = ref.rotation_mode
        if ref.rotation_mode == 'QUATERNION':
            ob.rotation_quaternion = ref.rotation_quaternion.copy()
        ob.matrix_world = ref.matrix_world.copy()
    for k, val in meta.items():
        ob[k] = val
    return ob


# --------------------------------------------------------------------------
# Union (elbows, tees: legs built separately, then merged)
# --------------------------------------------------------------------------

def _append(dst, src):
    vmap = {v: dst.verts.new(v.co) for v in src.verts}
    for f in src.faces:
        dst.faces.new([vmap[v] for v in f.verts])


def _components(bm):
    seen, comps = set(), []
    for f in bm.faces:
        if f in seen:
            continue
        stack, comp = [f], []
        seen.add(f)
        while stack:
            g = stack.pop()
            comp.append(g)
            for e in g.edges:
                for h in e.link_faces:
                    if h not in seen:
                        seen.add(h)
                        stack.append(h)
        comps.append(comp)
    return comps


def inside(bvh, p, direction=(0.5773, 0.5774, 0.5775), limit=64):
    """Parity test: odd number of crossings along a ray means inside."""
    d = Vector(direction).normalized()
    o = Vector(p)
    hits = 0
    for _ in range(limit):
        loc, _n, _i, dist = bvh.ray_cast(o, d)
        if loc is None:
            break
        hits += 1
        o = loc + d * max(dist * 1e-6, 1e-7)
    return hits % 2 == 1


def drop_interior(bm):
    """Delete any shell lying wholly inside another: that is interior."""
    comps = _components(bm)
    if len(comps) < 2:
        return 0
    doomed = []
    for i, comp in enumerate(comps):
        others = [f for j, c in enumerate(comps) if j != i for f in c]
        verts = list({v for f in others for v in f.verts})
        index = {v: n for n, v in enumerate(verts)}
        tree = BVHTree.FromPolygons([v.co for v in verts],
                                    [[index[v] for v in f.verts] for f in others])
        if inside(tree, comp[0].verts[0].co):
            doomed += comp
    if doomed:
        bmesh.ops.delete(bm, geom=doomed, context='FACES')
    return len(doomed)


def union(objs, target, col, lpm):
    """EXACT boolean union of rebuild shells into the target's local space."""
    inv = target.matrix_world.inverted()
    joined = bmesh.new()
    for ob in objs:
        b = bmesh.new()
        b.from_mesh(ob.data)
        b.transform(inv @ ob.matrix_world)
        orient(b)                      # every shell oriented on its own
        _append(joined, b)
        b.free()

    me = bpy.data.meshes.new("_XR_UNION")
    joined.to_mesh(me)
    joined.free()
    tmp = bpy.data.objects.new("_XR_UNION", me)
    tmp["xr_temp"] = True
    col.objects.link(tmp)
    tmp.matrix_world = target.matrix_world.copy()

    vl = bpy.context.view_layer
    prev = vl.objects.active
    for o in vl.objects:
        o.select_set(False)
    vl.objects.active = tmp
    tmp.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.intersect_boolean(operation='UNION', use_self=True, solver='EXACT')
    bpy.ops.object.mode_set(mode='OBJECT')

    bm = bmesh.new()
    bm.from_mesh(tmp.data)
    dropped = drop_interior(bm)
    cleanup(bm, lpm)
    bpy.data.objects.remove(tmp, do_unlink=True)
    bpy.data.meshes.remove(me)
    if prev is not None and prev.name in vl.objects:
        vl.objects.active = prev
    return bm, dropped


# --------------------------------------------------------------------------
# Solve: grid -> stations -> angles -> bmesh
# --------------------------------------------------------------------------

def solve(grid, s, tol_scale=1.0):
    st_tol = s.station_tol_mm * tol_scale
    an_tol = s.angle_tol_mm * tol_scale
    stations = geom.dp_stations(grid["z"], grid["rows"], st_tol, grid["forced"])
    rows = [grid["rows"][k] for k in stations]
    segs = geom.segments_for(grid["rmax_mm"], s.seg_small, s.seg_medium, s.seg_large)
    angles = geom.choose_angles(rows, grid["thetas"], an_tol, segs)
    return stations, angles
