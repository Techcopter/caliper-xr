"""Flatten large curved CAD faces whole, with Blender's angle-based unwrap.

Rebuilding a curved face needs a 2D domain in which it is one-to-one. A
whole-face unwrap gives that without cutting the face into artificial
pieces (whose zig-zag borders cost vertices and look messy). Seams go on
every CAD-face border, so each face is its own island. A face whose island
folds over itself is reported as unusable; it falls back to projection
pieces.
"""

import bmesh
import bpy
import numpy as np


def unwrap(cm, patch, want):
    """UV per triangle corner, (m,3,2), NaN where not unwrapped; and a bool
    per patch saying whether its island is fold-free."""
    m = len(cm.T)
    UV = np.full((m, 3, 2), np.nan)
    npatch = int(patch.max()) + 1
    good = np.zeros(npatch, bool)
    faces = np.flatnonzero(want[patch])
    if len(faces) == 0:
        return UV, good
    T = cm.T[faces]
    used, tloc = np.unique(T, return_inverse=True)
    tloc = tloc.reshape(-1, 3)
    me = bpy.data.meshes.new("_xr_uv")
    me.from_pydata(cm.V[used].tolist(), [], tloc.tolist())
    me.uv_layers.new(name="xr")
    bm = bmesh.new()
    bm.from_mesh(me)
    bm.faces.ensure_lookup_table()
    fp = patch[faces]
    for e in bm.edges:
        lf = e.link_faces
        e.seam = bool(len(lf) != 2 or fp[lf[0].index] != fp[lf[1].index])
    bm.to_mesh(me)
    bm.free()
    ob = bpy.data.objects.new("_xr_uv", me)
    ob["xr_temp"] = True
    bpy.context.scene.collection.objects.link(ob)
    vl = bpy.context.view_layer
    prev = vl.objects.active
    prev_mode = bpy.context.mode
    try:
        if prev_mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        for o in vl.objects:
            o.select_set(False)
        vl.objects.active = ob
        ob.select_set(True)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.uv.unwrap(method='ANGLE_BASED', fill_holes=True, margin=0.0)
        bpy.ops.object.mode_set(mode='OBJECT')
        uv = np.zeros(len(me.loops) * 2)
        me.uv_layers["xr"].data.foreach_get("uv", uv)
        uv = uv.reshape(-1, 3, 2)
    finally:
        bpy.data.objects.remove(ob, do_unlink=True)
        bpy.data.meshes.remove(me)
        if prev is not None and prev.name in vl.objects:
            vl.objects.active = prev
    UV[faces] = uv
    a = ((uv[:, 1, 0] - uv[:, 0, 0]) * (uv[:, 2, 1] - uv[:, 0, 1])
         - (uv[:, 1, 1] - uv[:, 0, 1]) * (uv[:, 2, 0] - uv[:, 0, 0]))
    for p in np.unique(fp):
        s = a[fp == p]
        pos, neg = int((s > 0).sum()), int((s < 0).sum())
        tiny = int((np.abs(s) <= np.abs(s).max() * 1e-12).sum()) if len(s) else 0
        # usable: consistently oriented (either way; a mirrored island is fine)
        good[p] = min(pos, neg) == 0 and tiny == 0
    return UV, good
