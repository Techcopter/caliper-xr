"""Outer shell: faces that border outside air stay, internal faces go.

A sealed box with a part inside, a box with a 1 mm slit, two blocks
touching face to face, and a block with a deep 5 x 40 mm drilled hole.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bmesh  # noqa: E402
import bpy  # noqa: E402
import numpy as np  # noqa: E402
from harness import MM, _cut, _cutter, check, done, reset  # noqa: E402

from caliper_xr import shell  # noqa: E402
from caliper_xr.visibility import world_tris  # noqa: E402


def box(bm, size, loc=(0, 0, 0), inward=False):
    vs = bmesh.ops.create_cube(bm, size=1.0)["verts"]
    bmesh.ops.scale(bm, vec=[s * MM for s in size], verts=vs)
    bmesh.ops.translate(bm, vec=[c * MM for c in loc], verts=vs)
    if inward:
        bmesh.ops.reverse_faces(bm, faces=list({f for v in vs for f in v.link_faces}))


def obj(name, *boxes):
    bm = bmesh.new()
    for b in boxes:
        box(bm, *b)
    me = bpy.data.meshes.new(name)
    bm.to_mesh(me)
    bm.free()
    ob = bpy.data.objects.new(name, me)
    bpy.context.scene.collection.objects.link(ob)
    return ob


reset()
A = obj("sealed_box", ((100, 100, 100),), ((98, 98, 98), (0, 0, 0), True))
A_in = obj("inside_a", ((30, 30, 30),))
B = obj("slit_box", ((100, 100, 100), (300, 0, 0)), ((98, 98, 98), (300, 0, 0), True))
_cut(B, _cutter('BOX', (1.0, 60, 10), (350, 0, 0)))
B_in = obj("inside_b", ((30, 30, 30), (300, 0, 0)))
C = obj("touching_pair", ((40, 40, 40), (600, 0, 0)), ((40, 40, 40), (640, 0, 0)))
D = obj("drilled_block", ((60, 60, 60), (900, 0, 0)))
_cut(D, _cutter('CYL', (5, 5, 80), (900, 0, 30), verts=48))
objs = [A, A_in, B, B_in, C, D]
bpy.context.view_layer.update()
dg = bpy.context.evaluated_depsgraph_get()
cen = dict(zip([o.name for o in objs], [x.mean(1) / MM for x in world_tris(objs, dg)]))


def frac(f, m):
    return float(f[m].mean()) if m.any() else 1.0


for gap, each in ((2.0, False), (2.0, True)):
    flags, hs = shell.analyse(objs, dg, gap * MM, each)
    tag = f"[{'each part' if each else 'selection'}, voxels {max(hs) / MM:.2g} mm]"
    f, c = flags["sealed_box"], cen["sealed_box"]
    outer = np.abs(c).max(1) > 49.5
    check(frac(f, outer) == 1.0, f"{tag} sealed box: outside faces kept")
    check(frac(f, ~outer) == 0.0, f"{tag} sealed box: inner walls left out")
    f, c = flags["touching_pair"], cen["touching_pair"]
    mid = np.abs(c[:, 0] - 620) < 0.01
    check(frac(f, mid) == 0.0 and frac(f, ~mid) == 1.0, f"{tag} touching faces left out, the rest kept")
    f, c = flags["drilled_block"], cen["drilled_block"]
    hole = np.hypot(c[:, 0] - 900, c[:, 1]) < 3.0
    need = 1.0 if each else 0.9
    check(frac(f, hole) >= need, f"{tag} deep 5 mm hole: {frac(f, hole):.0%} of its faces kept")
    if each:
        check(frac(flags["inside_a"], np.ones(len(flags["inside_a"]), bool)) == 1.0,
              f"{tag} a part judged on its own keeps its skin")
    else:
        check(frac(flags["inside_a"], np.ones(len(flags["inside_a"]), bool)) == 0.0,
              f"{tag} part inside the sealed box left out")
        check(frac(flags["inside_b"], np.ones(len(flags["inside_b"]), bool)) == 0.0,
              f"{tag} part behind a 1 mm slit left out (slit sealed)")
done("shell")
