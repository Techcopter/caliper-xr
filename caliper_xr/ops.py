"""Operators. None of them writes to an object outside the REBUILT collection.

The hard rule: the original CAD is never deleted, edited,
decimated, modified or moved. References are only ever *read* (evaluated
copy + private BVH). Every write goes to a new object in REBUILT.
"""

import time

import bmesh
import bpy
from mathutils import Matrix

from . import build, geom, sampling, verify

REPORT = "Caliper_Report"


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def report_text(lines_, reset=False):
    txt = bpy.data.texts.get(REPORT) or bpy.data.texts.new(REPORT)
    if reset:
        txt.clear()
    txt.write("\n".join(lines_) + "\n")
    for ln in lines_:
        print(ln)


def references(context):
    """Selected meshes that are CAD references (i.e. not our own rebuilds)."""
    obs = context.selected_objects or ([context.active_object] if context.active_object else [])
    return [o for o in obs if o and o.type == 'MESH' and not build.is_rebuild(o)]


def rebuilds(context):
    out = []
    for o in context.selected_objects:
        if o.type != 'MESH':
            continue
        if build.is_rebuild(o):
            out.append(o)
        else:
            rb = bpy.data.objects.get(build.rebuild_name(o))
            if rb is not None and build.is_rebuild(rb):
                out.append(rb)
    return list(dict.fromkeys(out))


def cursor_frame(context, ob):
    """3D cursor as a local-space axis: its Z is the axis, its origin a point on it."""
    inv = ob.matrix_world.inverted()
    cur = context.scene.cursor.matrix
    C = inv @ cur.translation
    a = (inv.to_3x3() @ (cur.to_3x3().col[2])).normalized()
    return (tuple(C), tuple(a))


def is_partial(s):
    return s.axis == 'CURSOR' or s.range_from > 0.0 or s.range_to < 1.0


def target_name(ref, s):
    """Whole part -> <name>_REBUILT (replaced on re-run). A leg or span ->
    <name>_REBUILT_partN, so building leg B never overwrites leg A."""
    base = build.rebuild_name(ref)
    if not is_partial(s):
        return base
    n = 1
    while bpy.data.objects.get(f"{base}_part{n}") is not None:
        n += 1
    return f"{base}_part{n}"


def budget(s, rmax_mm):
    if s.budget_class == 'AUTO':
        return geom.budget_for(rmax_mm)
    return verify.BUDGETS[s.budget_class]


def rebuild_one(context, ref, s, col):
    """Build ref's shell. Returns (object, info lines)."""
    dg = context.evaluated_depsgraph_get()
    lpm = 1.0 / sampling.mm_per_local(ref, context.scene)   # local per mm
    V, ref_polys, bvh = sampling.read_reference(ref, dg)
    if len(V) < 4:
        raise ValueError("reference has no geometry")

    mode = s.mode
    info = []
    grid = None
    if mode in {'AUTO', 'SECTION'}:
        cur = cursor_frame(context, ref) if s.axis == 'CURSOR' else None
        fr = sampling.frame(V, s.axis, cur)
        grid = sampling.sample_grid(bvh, V, fr, lpm, s)
        info.append(f"  rays {grid['rays']}  holes capped {grid['holes']} "
                    f"({grid['hole_cells']} samples)  misses {grid['miss'] * 100:.0f}%")
        if grid.get("extended"):
            info.append(f"  leg extended straight through the junction over "
                        f"{grid['extended']} station(s)")
        if mode == 'AUTO' and grid['miss'] > s.hull_miss_fraction:
            info.append("  section sampling missed too often -> convex hull fallback")
            mode = 'HULL'
        else:
            mode = 'SECTION'

    if mode == 'HULL':
        bm = build.hull(V)
        build.cleanup(bm, lpm)
        rmax = max(float(((V.max(0) - V.min(0)) ** 2).sum() ** 0.5) * 0.5 / lpm, 1.0)
        target = budget(s, rmax)
        used = "hull"
    else:
        target = budget(s, grid["rmax_mm"])
        bm, used = None, ""
        scales = (1.0, 1.5, 2.25, 3.4) if s.fit_budget else (1.0,)
        for sc in scales:
            if bm is not None:
                bm.free()
            stations, angles = build.solve(grid, s, sc)
            bm = build.loft(grid, stations, angles)
            build.cleanup(bm, lpm)
            used = f"{len(stations)} rings x {len(angles)} angles (tol x{sc:g})"
            if len(bm.faces) <= target:
                break

    meta = {"xr_source": ref.name, "xr_budget": target, "xr_mode": mode,
            "xr_built": time.strftime("%Y-%m-%d %H:%M"), "xr_partial": is_partial(s)}
    ob = build.to_object(bm, target_name(ref, s), ref, col, meta)
    info.insert(0, f"  {mode.lower()}: {used}; ref faces {ref_polys}")
    return ob, info


# --------------------------------------------------------------------------
# Operators
# --------------------------------------------------------------------------

class CALIPER_OT_analyse(bpy.types.Operator):
    """Dump size, components, axis, roundness and the holes that will be capped"""
    bl_idname = "caliper.analyse"
    bl_label = "Analyse"
    bl_options = {'REGISTER'}

    def execute(self, context):
        s = context.scene.caliper
        refs = references(context)
        if not refs:
            self.report({'WARNING'}, "Select one or more CAD parts")
            return {'CANCELLED'}
        dg = context.evaluated_depsgraph_get()
        out = [f"# Analyse  {time.strftime('%Y-%m-%d %H:%M')}"]
        for ref in refs:
            mm = sampling.mm_per_local(ref, context.scene)
            V, polys, bvh = sampling.read_reference(ref, dg)
            dims = (V.max(0) - V.min(0)) * mm
            fr = sampling.frame(V, 'AUTO')
            line = [f"== {ref.name}",
                    f"  verts {len(V)}  faces {polys}  components {sampling.component_count(ref)}",
                    f"  dims mm {'x'.join(f'{d:.2f}' for d in dims)}  scale "
                    f"{tuple(round(x, 7) for x in ref.matrix_world.to_scale())}",
                    f"  auto axis (local) {tuple(round(float(x), 3) for x in fr[1])}"]
            try:
                g = sampling.sample_grid(bvh, V, fr, 1.0 / mm, s)
                rows = g["rows"]
                rnd = sum(geom.is_round(r, s.angle_tol_mm) for r in rows) / len(rows)
                line.append(f"  max radius {g['rmax_mm']:.2f} mm  round sections {rnd * 100:.0f}%"
                            f"  budget {budget(s, g['rmax_mm'])}")
                line.append(f"  holes/openings to cap: {g['holes']} region(s)"
                            f"  ray misses {g['miss'] * 100:.0f}%")
            except ValueError as e:
                line.append(f"  sampling failed: {e}")
            out += line
        report_text(out)
        self.report({'INFO'}, f"Analysed {len(refs)} part(s) - see text '{REPORT}'")
        return {'FINISHED'}


class CALIPER_OT_rebuild(bpy.types.Operator):
    """Build one closed, hole-free, lean outer shell per selected part into REBUILT.
The original is only read, never modified"""
    bl_idname = "caliper.rebuild"
    bl_label = "Rebuild Selected"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        s = context.scene.caliper
        refs = references(context)
        if not refs:
            self.report({'WARNING'}, "Select one or more CAD parts (not rebuilds)")
            return {'CANCELLED'}
        if context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        col = build.rebuilt_collection(context.scene)
        wm = context.window_manager
        wm.progress_begin(0, len(refs))
        out = [f"# Rebuild  {time.strftime('%Y-%m-%d %H:%M')}"]
        ok = failed = 0
        dg = context.evaluated_depsgraph_get()
        for i, ref in enumerate(refs):
            wm.progress_update(i)
            try:
                ob, info = rebuild_one(context, ref, s, col)
                out.append(f"== {ref.name} -> {ob.name}")
                out += info
                if s.verify_after:
                    context.view_layer.update()
                    # A leg is only part of the reference: fit is meaningless
                    # until the union, so only its topology is checked here.
                    res = verify.run(None if ob.get("xr_partial") else ref,
                                     ob, dg, context.scene)
                    verify.store(ob, res)
                    out += verify.lines(ob.name, res)[1:]
                ok += 1
            except Exception as e:          # one bad part must not stop 1300
                failed += 1
                out.append(f"== {ref.name}: FAILED - {e}")
        wm.progress_end()
        report_text(out)
        self.report({'INFO' if not failed else 'WARNING'},
                    f"Rebuilt {ok}, failed {failed} - see text '{REPORT}'")
        return {'FINISHED'}


class CALIPER_OT_union(bpy.types.Operator):
    """EXACT boolean union of selected rebuild shells (elbow/tee legs) into the
active one; interior shells are deleted. Works only on REBUILT objects"""
    bl_idname = "caliper.union"
    bl_label = "Union Selected Rebuilds"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        objs = [o for o in context.selected_objects if o.type == 'MESH']
        target = context.active_object
        if len(objs) < 2 or target not in objs:
            self.report({'WARNING'}, "Select 2+ rebuilds, the active one names the result")
            return {'CANCELLED'}
        if not all(build.is_rebuild(o) for o in objs):
            self.report({'ERROR'}, "Union only runs on REBUILT objects - never on the CAD")
            return {'CANCELLED'}
        if context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        col = build.rebuilt_collection(context.scene)
        src = target.get("xr_source", "")
        ref = bpy.data.objects.get(src) if src else None
        lpm = 1.0 / sampling.mm_per_local(ref or target, context.scene)
        bm, dropped = build.union(objs, target, col, lpm)
        sources = {o.get("xr_source", "") for o in objs}
        name = build.rebuild_name(ref) if ref is not None and len(sources) == 1 else target.name
        meta = {k: target[k] for k in target.keys() if k.startswith("xr_")}
        meta["xr_mode"] = "UNION"
        meta["xr_partial"] = False
        mat = target.matrix_world.copy()
        for o in objs:
            me = o.data
            bpy.data.objects.remove(o, do_unlink=True)
            if me.users == 0:
                bpy.data.meshes.remove(me)
        ob = build.to_object(bm, name, ref, col, meta)
        if ref is None:
            ob.matrix_world = mat
        ob.select_set(True)
        context.view_layer.objects.active = ob
        self.report({'INFO'}, f"Union -> {ob.name}, {len(ob.data.polygons)} faces, "
                              f"{dropped} interior faces removed")
        return {'FINISHED'}


class CALIPER_OT_corner_sphere(bpy.types.Operator):
    """Add a sphere shell centred on the 3D cursor for an elbow's outer corner.
Offset elbows are a sharp polyline swept by a radius: the corner is a sphere
on the axis intersection, not a torus. Union it with the two legs"""
    bl_idname = "caliper.corner_sphere"
    bl_label = "Add Corner Sphere"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        s = context.scene.caliper
        act = context.active_object
        if act is None or act.type != 'MESH':
            self.report({'WARNING'}, "Make the CAD part (or one of its legs) active")
            return {'CANCELLED'}
        ref = bpy.data.objects.get(act.get("xr_source", "")) if build.is_rebuild(act) else act
        if ref is None:
            self.report({'WARNING'}, "Could not find the reference part")
            return {'CANCELLED'}
        lpm = 1.0 / sampling.mm_per_local(ref, context.scene)
        # Scaled a hair (x1.002) so nothing lands near-coincident with a leg.
        r = s.corner_radius_mm * 1.002
        segs = geom.segments_for(r, s.seg_small, s.seg_medium, s.seg_large)
        bm = bmesh.new()
        centre = ref.matrix_world.inverted() @ context.scene.cursor.location
        bmesh.ops.create_uvsphere(bm, u_segments=segs, v_segments=max(4, segs // 2),
                                  radius=r * lpm,
                                  matrix=Matrix.Translation(centre))
        build.cleanup(bm, lpm)
        col = build.rebuilt_collection(context.scene)
        base = build.rebuild_name(ref)
        n = 1
        while bpy.data.objects.get(f"{base}_part{n}") is not None:
            n += 1
        meta = {"xr_source": ref.name, "xr_mode": "CORNER",
                "xr_budget": verify.BUDGETS['ELBOW']}
        ob = build.to_object(bm, f"{base}_part{n}", ref, col, meta)
        self.report({'INFO'}, f"Added {ob.name} ({len(ob.data.polygons)} faces)")
        return {'FINISHED'}


class CALIPER_OT_verify(bpy.types.Operator):
    """Report topology, genus, budget and surface deviation both ways"""
    bl_idname = "caliper.verify"
    bl_label = "Verify"
    bl_options = {'REGISTER'}

    def execute(self, context):
        rbs = rebuilds(context)
        if not rbs:
            self.report({'WARNING'}, "Select rebuilds (or their CAD parts)")
            return {'CANCELLED'}
        build.rebuilt_collection(context.scene)
        context.view_layer.update()
        dg = context.evaluated_depsgraph_get()
        out = [f"# Verify  {time.strftime('%Y-%m-%d %H:%M')}"]
        bad = 0
        for rb in rbs:
            ref = bpy.data.objects.get(rb.get("xr_source", ""))
            res = verify.run(ref, rb, dg, context.scene)
            verify.store(rb, res)
            bad += int(not rb["xr_genus_ok"] or rb["xr_over_budget"])
            out += verify.lines(rb.name, res)
        report_text(out)
        self.report({'INFO' if not bad else 'WARNING'},
                    f"Verified {len(rbs)}, {bad} need attention - see text '{REPORT}'")
        return {'FINISHED'}


class CALIPER_OT_audit(bpy.types.Operator):
    """Audit every object in REBUILT against the shell checks (genus 0, one
component, face budget) and list the ones to redo. Changes nothing"""
    bl_idname = "caliper.audit"
    bl_label = "Audit REBUILT"
    bl_options = {'REGISTER'}

    def execute(self, context):
        col = bpy.data.collections.get(build.REBUILT)
        if col is None:
            self.report({'WARNING'}, "No REBUILT collection")
            return {'CANCELLED'}
        out = [f"# Audit  {time.strftime('%Y-%m-%d %H:%M')}"]
        redo = []
        for ob in col.all_objects:
            if ob.type != 'MESH':
                continue
            if ob.get("xr_mode") == "CAD_XR":
                tol = float(ob.get("xr_tol_mm", 0.2))
                dev = float(ob.get("xr_dev_max_mm", 0.0))
                if dev > tol + 1e-6:
                    redo.append(ob.name)
                    out.append(f"  REDO {ob.name}: deviation {dev:.3f} mm > {tol:.3f} mm")
                continue
            t = verify.topology(ob)
            b = verify.budget_of(ob)
            why = []
            if t["components"] != 1:
                why.append(f"{t['components']} components")
            if t["euler"] != 2:
                why.append(f"V-E+F={t['euler']} (hole/tunnel)")
            if t["non_manifold"]:
                why.append(f"{t['non_manifold']} non-manifold edges")
            if t["faces"] > b:
                why.append(f"{t['faces']} faces > {b}")
            if why:
                redo.append(ob.name)
                out.append(f"  REDO {ob.name}: {', '.join(why)}")
        out.append(f"  {len(redo)} to redo of {len(col.all_objects)}")
        report_text(out)
        self.report({'INFO'}, f"{len(redo)} rebuild(s) fail the shell checks - see '{REPORT}'")
        return {'FINISHED'}


class CALIPER_OT_housekeeping(bpy.types.Operator):
    """Delete this add-on's temp objects, purge orphan meshes, un-exclude REBUILT"""
    bl_idname = "caliper.housekeeping"
    bl_label = "Housekeeping"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        # Only objects this add-on tagged as temporary - never a user's '_' object.
        temps = [o for o in bpy.data.objects if o.get("xr_temp")]
        for o in temps:
            bpy.data.objects.remove(o, do_unlink=True)
        orphans = [m for m in bpy.data.meshes if m.users == 0]
        for m in orphans:
            bpy.data.meshes.remove(m)
        if bpy.data.collections.get(build.REBUILT):
            build.rebuilt_collection(context.scene)
        self.report({'INFO'}, f"Removed {len(temps)} temp objects, {len(orphans)} orphan meshes")
        return {'FINISHED'}


CLASSES = (CALIPER_OT_analyse, CALIPER_OT_rebuild, CALIPER_OT_union, CALIPER_OT_corner_sphere,
           CALIPER_OT_verify,
           CALIPER_OT_audit, CALIPER_OT_housekeeping)
