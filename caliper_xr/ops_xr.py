"""CAD -> XR operators: shape-preserving retessellation of the selection, and
a check of what the outer-shell and hidden-geometry tests would leave out."""

import time

import bpy
import numpy as np

from . import build, cadmesh, cadxr, sampling, segment, shell, visibility
from .ops import REPORT, report_text

QUALITY = {'HIGH': 0.1, 'BALANCED': 0.2, 'LIGHT': 0.5}
CHECK = "CALIPER_CHECK"
PREVIEW = "Caliper_Shell_Check"


def tolerance_mm(s):
    return s.xr_tol_mm if s.xr_quality == 'CUSTOM' else QUALITY[s.xr_quality]


class _S:
    """The few settings the pipeline reads, resolved from the scene."""

    def __init__(self, s):
        self.xr_tol_mm = tolerance_mm(s)
        self.xr_crease_deg = s.xr_crease_deg


def _refs(context):
    """Selected, visible CAD parts: not our rebuilds, not our previews."""
    return [o for o in context.selected_objects
            if o.type == 'MESH' and not build.is_rebuild(o) and not o.get("xr_temp")
            and o.visible_get()]


def _tests(context, refs, dg, out, shell_on):
    """Run the camera test (always on the whole selection) and the outer-shell
    test (on the selection or each part) the settings ask for. Returns
    (occ, shell flags) for cadxr, and says in `out` what ran."""
    s = context.scene.caliper
    each = s.xr_shell_scope == 'EACH'
    occ = flags = None
    if s.xr_remove_hidden:
        t0 = time.time()
        vis = visibility.visible_gpu(refs, dg)
        if vis is not None:
            occ = ("gpu", vis)
            out.append(f"  visibility: GPU, 128 views of {len(refs)} parts, {time.time() - t0:.1f}s")
        else:
            occ = ("cpu", visibility.occluder(refs, dg))
            out.append("  visibility: ray sampling (no GPU here; less certain)")
    if shell_on:
        t0 = time.time()
        wpm = 1.0 / (context.scene.unit_settings.scale_length * 1000.0)     # world units per mm
        flags, hs = shell.analyse(refs, dg, s.xr_shell_gap_mm * wpm, each)
        if hs:
            lo, hi = min(hs) / wpm, max(hs) / wpm
            size = f"{lo:.2g} mm" if hi - lo < 1e-6 else f"{lo:.2g}-{hi:.2g} mm"
            out.append(f"  outer shell: {'each part on its own' if each else 'selection as one'}, "
                       f"voxels {size}, openings under about {2 * hi:.2g} mm sealed, "
                       f"{time.time() - t0:.1f}s")
    return occ, flags


class CALIPER_OT_retessellate(bpy.types.Operator):
    """Rebuild the selected CAD parts lean for VR/AR: every CAD face is re-triangulated
from its simplified outline and stays within the tolerance of the original surface.
Flat faces become single polygons, crisp edges stay where the CAD has them, shading
is copied from the CAD. Parts are copied into REBUILT; the originals are never touched"""
    bl_idname = "caliper.retessellate"
    bl_label = "Retessellate for XR"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        s = context.scene.caliper
        refs = _refs(context)
        if not refs:
            self.report({'WARNING'}, "Select one or more visible CAD parts")
            return {'CANCELLED'}
        if context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        col = build.rebuilt_collection(context.scene)
        dg = context.evaluated_depsgraph_get()
        settings = _S(s)
        t0 = time.time()
        out = [f"# CAD -> XR  {time.strftime('%Y-%m-%d %H:%M')}  tolerance "
               f"{settings.xr_tol_mm:.3f} mm, {len(refs)} part(s)"]
        occ, flags = _tests(context, refs, dg, out, s.xr_shell_only)
        wm = context.window_manager
        wm.progress_begin(0, len(refs))
        src_total = new_total = hidden_parts = failed = 0
        made = []
        for i, ref in enumerate(refs):
            wm.progress_update(i)
            lpm = 1.0 / sampling.mm_per_local(ref, context.scene)
            try:
                bm, loopn, info, st = cadxr.run(ref, dg, settings, lpm, occ=occ, shell=flags)
            except Exception as e:                      # one bad part must not stop an assembly
                failed += 1
                out.append(f"== {ref.name}: FAILED - {e}")
                continue
            src_total += st.get("src_tris", 0)
            if bm is None:
                hidden_parts += 1
                out.append(f"== {ref.name}: not on the outside, skipped")
                out += info
                continue
            ob = cadxr.make_object(bm, loopn, build.rebuild_name(ref), ref, col, s.xr_keep_materials,
                                   {"xr_source": ref.name, "xr_mode": "CAD_XR",
                                    "xr_tol_mm": settings.xr_tol_mm,
                                    "xr_dev_max_mm": st.get("worst_mm", 0.0),
                                    "xr_built": time.strftime("%Y-%m-%d %H:%M")})
            tris = sum(len(p.vertices) - 2 for p in ob.data.polygons)
            new_total += tris
            ok = st.get("worst_mm", 0.0) <= settings.xr_tol_mm + 1e-6
            ob["xr_genus_ok"] = bool(ok)                    # export gate: shape held
            ob["xr_over_budget"] = False
            ob["xr_faces"] = len(ob.data.polygons)
            made.append(ob)
            out.append(f"== {ref.name} -> {ob.name}: {st.get('src_tris', 0)} -> {tris} triangles"
                       + ("" if ok else "  (deviation above tolerance: check it)"))
            out += info
        wm.progress_end()
        ratio = src_total / max(new_total, 1)
        out.append(f"  TOTAL {src_total} -> {new_total} triangles (x{ratio:.1f}), "
                   f"{hidden_parts} hidden part(s) skipped, {failed} failed, "
                   f"{time.time() - t0:.0f}s")
        report_text(out)
        for o in context.view_layer.objects:
            o.select_set(False)
        for o in made:
            o.select_set(True)
        if made:
            context.view_layer.objects.active = made[0]
        self.report({'INFO' if not failed else 'WARNING'},
                    f"{len(made)} part(s): {src_total} -> {new_total} triangles (x{ratio:.1f}); "
                    f"see '{REPORT}'")
        return {'FINISHED'}


class CALIPER_OT_toggle_source(bpy.types.Operator):
    """Hide or show the CAD originals of the selected rebuilds, to compare"""
    bl_idname = "caliper.toggle_source"
    bl_label = "Show / Hide Originals"
    bl_options = {'REGISTER'}

    def execute(self, context):
        col = bpy.data.collections.get(build.REBUILT)
        srcs = []
        for o in (col.all_objects if col else []):
            src = bpy.data.objects.get(o.get("xr_source", ""))
            if src is not None and src.name in context.view_layer.objects:
                srcs.append(src)
        if not srcs:
            self.report({'WARNING'}, "No rebuilds with originals in this view layer")
            return {'CANCELLED'}
        hide = not srcs[0].hide_get()
        for o in srcs:
            o.hide_set(hide)
        self.report({'INFO'}, f"{'Hid' if hide else 'Showed'} {len(srcs)} original(s)")
        return {'FINISHED'}


class CALIPER_OT_check_shell(bpy.types.Operator):
    """Check what the rebuild would leave out of the selected parts - faces inside
their outer shell, and faces no one can see when Remove hidden geometry is on -
without building anything. The report lists them per part, and a red preview
object holds them (hide the parts or isolate it to look); Housekeeping deletes it"""
    bl_idname = "caliper.check_shell"
    bl_label = "Check Shell"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        s = context.scene.caliper
        refs = _refs(context)
        if not refs:
            self.report({'WARNING'}, "Select one or more visible CAD parts")
            return {'CANCELLED'}
        if context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        dg = context.evaluated_depsgraph_get()
        t0 = time.time()
        scope = "each part on its own" if s.xr_shell_scope == 'EACH' else "selection as one"
        out = [f"# Shell check  {time.strftime('%Y-%m-%d %H:%M')}  seal gaps under "
               f"{s.xr_shell_gap_mm:g} mm, {scope}, {len(refs)} part(s)"]
        occ, flags = _tests(context, refs, dg, out, True)
        total = hid_t = inn_t = gone = 0
        rows, verts, tris = [], [], []
        for ref in refs:
            cm = cadmesh.load(ref, dg)
            if not len(cm.T):
                continue
            patch = segment.patches(cm, segment.hard_edges(cm, s.xr_crease_deg))
            hidden, inner = cadxr.left_out(cm, patch, ref, occ, flags)
            h, i = int(hidden[patch].sum()), int(inner[patch].sum())
            total, hid_t, inn_t = total + len(cm.T), hid_t + h, inn_t + i
            if not h + i:
                continue
            drop = hidden | inner
            gone += int(drop.all())
            rows.append((h + i, f"  {ref.name}: {h + i} of {len(cm.T)} triangles "
                                f"({int(drop.sum())} of {len(drop)} CAD faces) - {h} hidden, "
                                f"{i} inside the shell" + ("; nothing left" if drop.all() else "")))
            M = np.array(ref.matrix_world)
            T = cm.T[drop[patch]]
            used, inv = np.unique(T, return_inverse=True)
            tris.append(inv.reshape(-1, 3) + sum(len(v) for v in verts))
            verts.append(cm.V[used] @ M[:3, :3].T + M[:3, 3])
        out += [r for _n, r in sorted(rows, key=lambda r: -r[0])]
        out.append(f"  TOTAL {total} triangles: {hid_t + inn_t} left out ({hid_t} hidden, {inn_t} "
                   f"inside the shell; {gone} part(s) entirely), {total - hid_t - inn_t} kept, "
                   f"{time.time() - t0:.0f}s")
        ob = _preview(context, np.concatenate(verts) if verts else np.zeros((0, 3)),
                      np.concatenate(tris) if tris else np.zeros((0, 3), np.int64))
        if ob is not None:
            out.append(f"  preview: '{ob.name}' holds what is left out, in red - it sits inside the "
                       f"parts, so hide them or isolate it (numpad /) to look")
        report_text(out)
        self.report({'INFO'}, f"{hid_t + inn_t} of {total} triangles would be left out "
                              f"({inn_t} inside the shell); see '{REPORT}'")
        return {'FINISHED'}


def _preview(context, V, T):
    """Replace the red preview of left-out faces (a temp object, never exported)."""
    old = bpy.data.objects.get(PREVIEW)
    if old is not None and old.get("xr_temp"):
        me_old = old.data
        bpy.data.objects.remove(old, do_unlink=True)
        if me_old is not None and me_old.users == 0:
            bpy.data.meshes.remove(me_old)
    if not len(T):
        return None
    me = bpy.data.meshes.new(PREVIEW)
    me.from_pydata(V.tolist(), [], T.tolist())
    mat = bpy.data.materials.get(PREVIEW) or bpy.data.materials.new(PREVIEW)
    mat.diffuse_color = (0.9, 0.12, 0.08, 1.0)
    me.materials.append(mat)
    col = bpy.data.collections.get(CHECK)
    if col is None:
        col = bpy.data.collections.new(CHECK)
    if col.name not in context.scene.collection.children:
        context.scene.collection.children.link(col)
    ob = bpy.data.objects.new(PREVIEW, me)
    col.objects.link(ob)
    ob["xr_temp"] = True
    ob.hide_render = True
    return ob


CLASSES = (CALIPER_OT_retessellate, CALIPER_OT_toggle_source, CALIPER_OT_check_shell)
