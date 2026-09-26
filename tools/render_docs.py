"""Render the README images from procedural parts (no customer geometry).

    blender -b --factory-startup --python tools/render_docs.py

docs/img/before_after.png  CAD tessellation (top) and Caliper XR output (bottom)
docs/img/shell_check.png   a sealed enclosure seen through, internal parts in red
"""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests"))

import bmesh  # noqa: E402
import bpy  # noqa: E402
import numpy as np  # noqa: E402
from harness import MM, ROOT, bracket, housing, reset, shaft  # noqa: E402
from mathutils import Vector  # noqa: E402

import caliper_xr  # noqa: E402
from caliper_xr import build, cadxr, sampling  # noqa: E402

OUT = os.path.join(ROOT, "docs", "img")
STEEL, YELLOW, RED = (0.46, 0.49, 0.53, 1), (0.95, 0.72, 0.02, 1), (0.85, 0.16, 0.10, 1)


class Settings:
    xr_tol_mm = 0.2
    xr_crease_deg = 30.0


def stage():
    sc = bpy.context.scene
    sc.render.engine = 'BLENDER_WORKBENCH'
    sc.world = sc.world or bpy.data.worlds.new("w")
    sc.world.color = (0.012, 0.014, 0.017)
    sh = sc.display.shading
    sh.light, sh.color_type, sh.show_cavity = 'STUDIO', 'OBJECT', True
    sh.cavity_type, sh.show_specular_highlight = 'BOTH', True
    sc.render.film_transparent = False
    return sc


def wire(ob, width_mm=0.22):
    """A yellow wireframe riding on the part: the edges the XR device draws."""
    w = bpy.data.objects.new(ob.name + "_wire", ob.data.copy())
    bpy.context.scene.collection.objects.link(w)
    w.matrix_world = ob.matrix_world.copy()
    m = w.modifiers.new("w", 'WIREFRAME')
    # no even offset: it spikes at the long thin triangles CAD is made of
    m.thickness, m.use_replace, m.use_even_offset = width_mm * MM, True, False
    w.color = YELLOW
    ob.color = STEEL
    return w


def camera(objs, az, el, res):
    sc = bpy.context.scene
    P = np.array([o.matrix_world @ Vector(c) for o in objs for c in o.bound_box])
    lo, hi = P.min(0), P.max(0)
    ctr, rad = Vector(((lo + hi) / 2).tolist()), float(np.linalg.norm(hi - lo)) / 2
    cam = bpy.data.objects.get("cam") or bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
    if cam.name not in sc.collection.objects:
        sc.collection.objects.link(cam)
    cam.data.type, cam.data.ortho_scale = 'ORTHO', rad * 2.45
    cam.data.clip_start, cam.data.clip_end = rad * 0.01, rad * 20
    a, e = math.radians(az), math.radians(el)
    d = Vector((math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)))
    cam.location = ctr + d * rad * 4
    cam.rotation_euler = (-d).to_track_quat('-Z', 'Y').to_euler()
    sc.camera = cam
    sc.render.resolution_x, sc.render.resolution_y = res


def render(path, show):
    for o in bpy.context.scene.objects:
        o.hide_render = o.type == 'MESH' and o not in show
    bpy.context.scene.render.filepath = path
    bpy.ops.render.render(write_still=True)


def before_after():
    reset()
    stage()
    parts = [bracket((0, 0, 0)), shaft((95, 0, 0)), housing((190, 0, 0))]
    col = build.rebuilt_collection(bpy.context.scene)
    outs = []
    for ob in parts:
        bpy.context.view_layer.update()
        dg = bpy.context.evaluated_depsgraph_get()
        bm, loopn, _info, _st = cadxr.run(ob, dg, Settings, 1.0 / sampling.mm_per_local(ob, bpy.context.scene))
        outs.append(cadxr.make_object(bm, loopn, ob.name + "_XR", ob, col, False, {}))
    src_w = [wire(o) for o in parts]
    out_w = [wire(o) for o in outs]
    camera(parts, -58, 28, (1600, 720))
    tmp = [os.path.join(OUT, "_a.png"), os.path.join(OUT, "_b.png")]
    render(tmp[0], parts + src_w)
    render(tmp[1], outs + out_w)
    stack(tmp, os.path.join(OUT, "before_after.png"))
    n_src = sum(len(p.data.polygons) for p in parts)
    n_out = sum(sum(len(f.vertices) - 2 for f in o.data.polygons) for o in outs)
    print(f"before_after: {n_src} -> {n_out} triangles")


def shell_check():
    reset()
    stage()
    caliper_xr.register()
    try:
        bm = bmesh.new()
        vs = bmesh.ops.create_cube(bm, size=1.0)["verts"]               # a sealed skin
        bmesh.ops.scale(bm, vec=(160 * MM, 110 * MM, 90 * MM), verts=vs)
        me = bpy.data.meshes.new("enclosure")
        bm.to_mesh(me)
        bm.free()
        box = bpy.data.objects.new("enclosure", me)
        bpy.context.scene.collection.objects.link(box)
        inner = [bracket((-30, -18, -25)), shaft((40, 20, -10)), housing((-35, 25, 5))]
        for o in inner:
            o.scale = (0.8, 0.8, 0.8)
        bpy.context.view_layer.update()
        for o in bpy.context.view_layer.objects:
            o.select_set(o in inner + [box])
        s = bpy.context.scene.caliper
        s.xr_shell_only, s.xr_remove_hidden = True, False
        bpy.ops.caliper.check_shell()
        prev = bpy.data.objects["Caliper_Shell_Check"]
        prev.color, box.color = RED, STEEL
        sh = bpy.context.scene.display.shading
        sh.show_xray, sh.xray_alpha = True, 0.35
        camera([box], -52, 30, (1200, 760))
        render(os.path.join(OUT, "shell_check.png"), [box, prev])
    finally:
        caliper_xr.unregister()


def stack(paths, out):
    imgs = []
    for p in paths:
        im = bpy.data.images.load(p, check_existing=False)
        imgs.append(np.array(im.pixels[:]).reshape(im.size[1], im.size[0], 4))
        bpy.data.images.remove(im)
        os.remove(p)
    band = np.zeros((6, imgs[0].shape[1], 4))
    band[..., :] = (0.95, 0.72, 0.02, 1)
    grid = np.concatenate([imgs[1], band, imgs[0]], 0)        # image rows run bottom-up
    im = bpy.data.images.new("stack", grid.shape[1], grid.shape[0])
    im.pixels = grid.ravel().tolist()
    im.filepath_raw, im.file_format = out, 'PNG'
    im.save()


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    before_after()
    shell_check()
