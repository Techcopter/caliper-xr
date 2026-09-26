"""GLB export for VR / AR / WebXR, with a pre-flight report."""

import os
import time

import bmesh
import bpy
from bpy.props import BoolProperty, EnumProperty, IntProperty, StringProperty
from bpy_extras.io_utils import ExportHelper

from . import build
from .ops import REPORT, report_text


class CALIPER_OT_export(bpy.types.Operator, ExportHelper):
    """Export rebuilds as a GLB for VR / AR / WebXR, with a pre-flight report"""
    bl_idname = "caliper.export_glb"
    bl_label = "Export XR GLB"
    bl_options = {'REGISTER'}

    filename_ext = ".glb"
    filter_glob: StringProperty(default="*.glb", options={'HIDDEN'})

    scope: EnumProperty(name="Objects", items=[
        ('REBUILT', "All REBUILT", "Every visible object in the REBUILT collection"),
        ('SELECTED', "Selected rebuilds", "Only selected objects in REBUILT")],
        default='REBUILT')
    draco: BoolProperty(name="Draco compression", default=True,
                        description="Smaller file; the viewer must ship a Draco decoder")
    draco_level: IntProperty(name="Draco level", default=6, min=0, max=10)
    merge: BoolProperty(name="Merge into one mesh", default=False,
                        description="One draw call, but part names are lost: parts can "
                                    "no longer be highlighted or selected individually")
    y_up: BoolProperty(name="+Y up", default=True)
    require_verified: BoolProperty(name="Block unverified / failing parts", default=True)
    write_report: BoolProperty(name="Write report beside GLB", default=True)

    def targets(self, context):
        col = bpy.data.collections.get(build.REBUILT)
        if col is None:
            return []
        obs = [o for o in col.all_objects if o.type == 'MESH' and o.visible_get()]
        if self.scope == 'SELECTED':
            obs = [o for o in obs if o.select_get()]
        return obs

    def execute(self, context):
        s = context.scene.caliper
        obs = self.targets(context)
        if not obs:
            self.report({'WARNING'}, "Nothing to export in REBUILT")
            return {'CANCELLED'}
        if context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')

        tris = sum(sum(len(p.vertices) - 2 for p in o.data.polygons) for o in obs)
        failing = [o.name for o in obs if not o.get("xr_genus_ok", False)]
        over = [o.name for o in obs if o.get("xr_over_budget", False)]
        lines = [f"# XR export  {time.strftime('%Y-%m-%d %H:%M')}  -> {self.filepath}",
                 f"  objects {len(obs)}  triangles {tris}  (scene budget {s.scene_tri_budget})",
                 f"  unverified/failing {len(failing)}  over face budget {len(over)}"]
        if tris > s.scene_tri_budget:
            lines.append("  WARNING: over the scene triangle budget for mobile XR")
        lines += [f"    not verified/failing: {n}" for n in failing[:50]]
        if failing and self.require_verified:
            report_text(lines + ["  BLOCKED: run Verify, fix the failing parts, or untick the block"])
            self.report({'ERROR'}, f"{len(failing)} part(s) unverified or failing - see '{REPORT}'")
            return {'CANCELLED'}

        vl = context.view_layer
        prev_sel = [o for o in vl.objects if o.select_get()]
        prev_act = vl.objects.active
        tmp = None
        try:
            for o in vl.objects:
                o.select_set(False)
            if self.merge:
                tmp = self._merged(context, obs)
                tmp.select_set(True)
                vl.objects.active = tmp
            else:
                for o in obs:
                    o.select_set(True)
                vl.objects.active = obs[0]
            self._gltf(self.filepath)
        finally:
            if tmp is not None:
                me = tmp.data
                bpy.data.objects.remove(tmp, do_unlink=True)
                bpy.data.meshes.remove(me)
            for o in vl.objects:
                o.select_set(False)
            for o in prev_sel:
                if o.name in vl.objects:
                    o.select_set(True)
            if prev_act is not None and prev_act.name in vl.objects:
                vl.objects.active = prev_act

        size = os.path.getsize(self.filepath) if os.path.exists(self.filepath) else 0
        lines.append(f"  wrote {size / 1024:.0f} KB")
        report_text(lines)
        if self.write_report:
            with open(os.path.splitext(self.filepath)[0] + ".report.txt", "w") as fh:
                fh.write("\n".join(lines) + "\n")
        self.report({'INFO'}, f"Exported {len(obs)} objects, {tris} tris, {size / 1024:.0f} KB")
        return {'FINISHED'}

    def _merged(self, context, obs):
        bm = bmesh.new()
        for o in obs:
            b = bmesh.new()
            b.from_mesh(o.data)
            b.transform(o.matrix_world)
            vmap = {v: bm.verts.new(v.co) for v in b.verts}
            for f in b.faces:
                nf = bm.faces.new([vmap[v] for v in f.verts])
                nf.smooth = f.smooth
            for e in b.edges:
                if not e.smooth:
                    ne = bm.edges.get([vmap[v] for v in e.verts])
                    if ne:
                        ne.smooth = False
            b.free()
        me = bpy.data.meshes.new("_XR_MERGED")
        bm.to_mesh(me)
        bm.free()
        ob = bpy.data.objects.new("XR_Merged", me)
        ob["xr_temp"] = True
        build.rebuilt_collection(context.scene).objects.link(ob)
        return ob

    def _gltf(self, path):
        want = dict(filepath=path, export_format='GLB', use_selection=True,
                    export_apply=True, export_yup=self.y_up, export_normals=True,
                    export_texcoords=False, export_materials='EXPORT',
                    export_draco_mesh_compression_enable=self.draco,
                    export_draco_mesh_compression_level=self.draco_level,
                    export_draco_position_quantization=14,
                    export_draco_normal_quantization=10)
        # Exporter options change between Blender versions: pass only the
        # ones this build actually has.
        have = {p.identifier for p in bpy.ops.export_scene.gltf.get_rna_type().properties}
        bpy.ops.export_scene.gltf(**{k: v for k, v in want.items() if k in have})


CLASSES = (CALIPER_OT_export,)
