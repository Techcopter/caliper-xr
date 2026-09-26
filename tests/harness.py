"""Shared set-up for the headless tests.

Imports the add-on from this checkout, builds CAD-like parts procedurally
(fillets, drilled holes, pockets - dense, fan-triangulated, like a CAD
tessellation), and fails loudly: every check raises, and the runner starts
Blender with --python-exit-code so a failure is a non-zero exit.
"""

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import bpy  # noqa: E402

import caliper_xr  # noqa: E402,F401

MM = 0.001
RESULTS = []


def reset():
    """An empty scene in metres, whatever the startup file holds."""
    bpy.ops.wm.read_factory_settings(use_empty=True)
    bpy.context.scene.unit_settings.scale_length = 1.0


def check(ok, what):
    RESULTS.append((bool(ok), what))
    print(("  PASS  " if ok else "  FAIL  ") + what)
    if not ok:
        raise AssertionError(what)


def _apply_all(ob):
    bpy.context.view_layer.objects.active = ob
    for m in list(ob.modifiers):
        bpy.ops.object.modifier_apply(modifier=m.name)


def _cutter(kind, size, loc, rot=(0.0, 0.0, 0.0), verts=64):
    if kind == 'CYL':
        bpy.ops.mesh.primitive_cylinder_add(radius=size[0] * MM / 2, depth=size[2] * MM,
                                            vertices=verts, location=[v * MM for v in loc],
                                            rotation=rot)
    else:
        bpy.ops.mesh.primitive_cube_add(size=1.0, location=[v * MM for v in loc], rotation=rot)
        bpy.context.object.scale = [v * MM for v in size]
    return bpy.context.object


def _cut(ob, cutter):
    m = ob.modifiers.new("cut", 'BOOLEAN')
    m.operation, m.object, m.solver = 'DIFFERENCE', cutter, 'EXACT'
    _apply_all(ob)
    bpy.data.objects.remove(cutter, do_unlink=True)


def _finish(ob, name, fillet_mm, segments):
    if fillet_mm:
        b = ob.modifiers.new("fillet", 'BEVEL')
        b.width, b.segments, b.limit_method = fillet_mm * MM, segments, 'ANGLE'
        _apply_all(ob)
    t = ob.modifiers.new("tess", 'TRIANGULATE')
    t.quad_method, t.ngon_method = 'FIXED', 'CLIP'
    _apply_all(ob)
    ob.name = ob.data.name = name
    return ob


def bracket(loc=(0.0, 0.0, 0.0)):
    """80 x 40 x 10 mm plate, 3 mm fillets, two 8 mm through holes."""
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=[v * MM for v in loc])
    ob = bpy.context.object
    ob.scale = (80 * MM, 40 * MM, 10 * MM)
    bpy.ops.object.transform_apply(scale=True)
    _finish(ob, "bracket", 3.0, 8)
    for x in (-25.0, 25.0):
        _cut(ob, _cutter('CYL', (8, 8, 30), (loc[0] + x, loc[1], loc[2])))
    return _finish(ob, "bracket", 0, 0)


def shaft(loc=(0.0, 0.0, 0.0)):
    """24 mm shaft, 60 mm long, 1 mm chamfers, 6 mm cross hole."""
    bpy.ops.mesh.primitive_cylinder_add(radius=12 * MM, depth=60 * MM, vertices=128,
                                        location=[v * MM for v in loc])
    ob = _finish(bpy.context.object, "shaft", 1.0, 6)
    _cut(ob, _cutter('CYL', (6, 6, 40), (loc[0], loc[1], loc[2] + 10), rot=(math.pi / 2, 0, 0)))
    return _finish(ob, "shaft", 0, 0)


def housing(loc=(0.0, 0.0, 0.0)):
    """60 x 60 x 40 mm block, 4 mm fillets, 40 x 40 x 30 mm open pocket."""
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=[v * MM for v in loc])
    ob = bpy.context.object
    ob.scale = (60 * MM, 60 * MM, 40 * MM)
    bpy.ops.object.transform_apply(scale=True)
    _finish(ob, "housing", 4.0, 10)
    _cut(ob, _cutter('BOX', (40, 40, 30), (loc[0], loc[1], loc[2] + 10)))
    return _finish(ob, "housing", 0, 0)


def triangles(ob):
    return sum(len(p.vertices) - 2 for p in ob.data.polygons)


def done(name):
    n = sum(ok for ok, _ in RESULTS)
    print(f"== {name}: {n}/{len(RESULTS)} checks passed")
