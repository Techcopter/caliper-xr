"""Caliper XR - heavy CAD in, lean XR out, measured both ways.

Main method, "Retessellate for XR": every CAD face is recovered from the
tessellated mesh and rebuilt lean from its simplified outline, never farther
than the tolerance from the original surface. Flat faces become single
polygons, crisp edges stay where the CAD has them, shading comes from the
CAD's own normals, materials are kept, and geometry nobody can see from
outside the selection is removed - or, with "Outer shell only", everything
inside the outer skin, which "Check Shell" shows before anything is built.

Legacy method (collapsed in the panel): the outer-shell rebuild of the
section-loft method, for revolved parts.

Everything is written to the REBUILT collection; the original CAD is only
ever read. Sidebar: View3D > N panel > Caliper.
"""

import bpy
from bpy.props import (BoolProperty, EnumProperty, FloatProperty,
                       IntProperty, PointerProperty)

from . import export, ops, ops_xr


class CALIPER_Settings(bpy.types.PropertyGroup):
    # --- CAD -> XR -------------------------------------------------------
    xr_quality: EnumProperty(name="Quality", default='BALANCED', items=[
        ('HIGH', "High (0.1 mm)", "Close-up inspection; heaviest"),
        ('BALANCED', "Balanced (0.2 mm)", "Indistinguishable at arm's length"),
        ('LIGHT', "Light (0.5 mm)", "Room-scale scenes, mobile headsets"),
        ('CUSTOM', "Custom", "Set the tolerance yourself")])
    xr_tol_mm: FloatProperty(name="Tolerance (mm)", default=0.2, min=0.005, max=20.0,
                             description="Largest allowed distance from the original surface")
    xr_crease_deg: FloatProperty(name="Crease angle", default=30.0, min=5.0, max=89.0,
                                 description="Edges sharper than this are always kept as crisp "
                                             "edges, even where the CAD's normals are smooth")
    xr_remove_hidden: BoolProperty(name="Remove hidden geometry", default=True,
                                   description="Drop CAD faces no one can see from outside the "
                                               "selected parts (buried brackets, inner walls). "
                                               "A face counts as visible if any point of it can "
                                               "be seen from any direction")
    xr_shell_only: BoolProperty(name="Outer shell only", default=False,
                                description="Build only the outer skin: leave out faces inside it - "
                                            "inner bodies, walls of closed cavities, faces between "
                                            "touching bodies, and whatever sits behind a seam or "
                                            "vent narrower than the gap. Check Shell shows them")
    xr_shell_gap_mm: FloatProperty(name="Seal gaps under (mm)", default=2.0, min=0.0, max=100.0,
                                   description="Openings narrower than this are not a way in: what "
                                               "is only reachable through them is left out. Large "
                                               "assemblies are limited by the grid; the report "
                                               "says what was sealed")
    xr_shell_scope: EnumProperty(name="Shell of", default='ALL', items=[
        ('ALL', "Whole selection", "One skin around all the selected parts: faces the other "
                                   "parts cover are left out too"),
        ('EACH', "Each part", "Every part is judged on its own and keeps its complete skin - "
                              "for parts that move, open or come apart in XR. Turn off Remove "
                              "hidden geometry too: that test always judges the whole selection")])
    xr_keep_materials: BoolProperty(name="Keep materials", default=True)
    show_legacy: BoolProperty(name="Legacy: outer-shell rebuild", default=False)

    # --- Legacy: outer-shell rebuild (revolved parts) ---------------------
    mode: EnumProperty(name="Method", default='AUTO', items=[
        ('AUTO', "Auto", "Section loft; convex hull if the part is not sampleable"),
        ('SECTION', "Section loft", "Ray-sample rings along an axis (cylinders, hex, steps)"),
        ('HULL', "Convex hull", "Coarse fallback for parts no axis describes")])
    axis: EnumProperty(name="Axis", default='AUTO', items=[
        ('AUTO', "Auto", "Principal symmetry axis, snapped to a local axis within 5 deg"),
        ('X', "Local X", ""), ('Y', "Local Y", ""), ('Z', "Local Z", ""),
        ('CURSOR', "3D Cursor Z", "Axis along the cursor's Z through its origin - "
                                  "use for elbow / tee legs")])
    range_from: FloatProperty(name="From", default=0.0, min=0.0, max=1.0, subtype='FACTOR')
    range_to: FloatProperty(name="To", default=1.0, min=0.0, max=1.0, subtype='FACTOR')
    ray_radius_mm: FloatProperty(name="Ray start radius (mm)", default=0.0, min=0.0,
                                 description="0 = auto. Set just outside a leg's radius so "
                                             "rays near a junction do not hit the other leg")
    corner_radius_mm: FloatProperty(name="Corner radius (mm)", default=10.0, min=0.01)
    angular_samples: IntProperty(name="Angular samples", default=360, min=24, max=1440)
    sample_step_mm: FloatProperty(name="Station step (mm)", default=0.1, min=0.01, max=10.0)
    max_stations: IntProperty(name="Max stations", default=400, min=8, max=4000)
    station_tol_mm: FloatProperty(name="Station tol (mm)", default=0.15, min=0.001, max=5.0)
    angle_tol_mm: FloatProperty(name="Section tol (mm)", default=0.1, min=0.001, max=5.0)
    step_mm: FloatProperty(name="Hard step (mm)", default=0.3, min=0.01, max=10.0)
    seg_small: IntProperty(name="Small", default=8, min=3, max=64, description="r < 10 mm")
    seg_medium: IntProperty(name="Medium", default=16, min=3, max=128, description="r < 60 mm")
    seg_large: IntProperty(name="Large", default=32, min=3, max=256, description="r >= 60 mm")
    fill_mode: EnumProperty(name="Fill", default='PITS', items=[
        ('PITS', "Holes & pockets", "Cap bolt holes, ports, pockets; keep grooves and slots"),
        ('AGGRESSIVE', "Also slots & grooves", "Cap anything bounded in either direction")])
    max_hole_mm: FloatProperty(name="Max hole (mm)", default=25.0, min=0.5, max=1000.0)
    max_hole_deg: FloatProperty(name="Max hole (deg)", default=40.0, min=2.0, max=180.0)
    budget_class: EnumProperty(name="Budget", default='AUTO', items=[
        ('AUTO', "Auto (by size)", ""), ('NUT', "Nut / plug / cap (150)", ""),
        ('FITTING', "Fitting / adapter (400)", ""), ('ELBOW', "Elbow / tee (800)", ""),
        ('HOUSING', "Housing / motor (2000)", "")])
    fit_budget: BoolProperty(name="Loosen to fit budget", default=False)
    hull_miss_fraction: FloatProperty(name="Hull if misses >", default=0.35, min=0.0, max=1.0,
                                      subtype='FACTOR')
    verify_after: BoolProperty(name="Verify after rebuild", default=True)
    scene_tri_budget: IntProperty(name="Scene triangle budget", default=300000, min=1000,
                                  description="Warn above this on export (mobile headset class)")


class CALIPER_PT_panel(bpy.types.Panel):
    bl_label = "Caliper XR"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Caliper"

    def draw(self, context):
        s = context.scene.caliper
        L = self.layout

        b = L.box()
        b.label(text="CAD -> XR (shape-true)", icon='MESH_ICOSPHERE')
        b.prop(s, "xr_quality")
        if s.xr_quality == 'CUSTOM':
            b.prop(s, "xr_tol_mm")
        b.prop(s, "xr_crease_deg")
        b.prop(s, "xr_remove_hidden")
        b.prop(s, "xr_shell_only")
        col = b.column(align=True)
        col.active = s.xr_shell_only
        col.prop(s, "xr_shell_gap_mm")
        col.prop(s, "xr_shell_scope", text="")
        b.prop(s, "xr_keep_materials")
        r = b.row()
        r.scale_y = 1.4
        r.operator("caliper.retessellate", icon='MOD_DECIM')
        r = b.row(align=True)
        r.operator("caliper.check_shell", icon='VIEWZOOM')
        r.operator("caliper.toggle_source", icon='HIDE_OFF')

        b = L.box()
        b.label(text="Check & export", icon='EXPORT')
        r = b.row(align=True)
        r.operator("caliper.verify", icon='CHECKMARK')
        r.operator("caliper.audit", icon='VIEWZOOM')
        b.prop(s, "scene_tri_budget")
        b.operator("caliper.export_glb", icon='EXPORT')
        b.operator("caliper.housekeeping", icon='TRASH')

        ob = context.active_object
        if ob is not None and "xr_faces" in ob:
            b = L.box()
            b.label(text=ob.name)
            b.label(text=f"Faces {ob['xr_faces']}")
            if "xr_dev_max_mm" in ob:
                ok = ob.get("xr_genus_ok")
                b.label(text=f"Largest deviation {ob['xr_dev_max_mm']:.3f} mm",
                        icon='CHECKMARK' if ok else 'ERROR')

        b = L.box()
        b.prop(s, "show_legacy", icon='TRIA_DOWN' if s.show_legacy else 'TRIA_RIGHT', emboss=False)
        if not s.show_legacy:
            return
        b.label(text="Closed outer shells for revolved parts")
        col = b.column(align=True)
        col.operator("caliper.analyse", icon='VIEWZOOM')
        col.operator("caliper.rebuild", icon='MOD_REMESH')
        col.operator("caliper.union", icon='MOD_BOOLEAN')
        b.prop(s, "mode", text="")
        b.prop(s, "axis")
        r = b.row(align=True)
        r.prop(s, "range_from")
        r.prop(s, "range_to")
        b.prop(s, "ray_radius_mm")
        r = b.row(align=True)
        r.prop(s, "corner_radius_mm")
        r.operator("caliper.corner_sphere", text="", icon='SPHERE')
        b.prop(s, "budget_class")
        b.prop(s, "fit_budget")
        b.prop(s, "verify_after")
        b.prop(s, "fill_mode", text="")
        b.prop(s, "max_hole_mm")
        b.prop(s, "max_hole_deg")
        b.prop(s, "station_tol_mm")
        b.prop(s, "angle_tol_mm")
        b.prop(s, "step_mm")
        r = b.row(align=True)
        r.label(text="Segments")
        r.prop(s, "seg_small", text="S")
        r.prop(s, "seg_medium", text="M")
        r.prop(s, "seg_large", text="L")
        b.prop(s, "angular_samples")
        b.prop(s, "sample_step_mm")
        b.prop(s, "max_stations")
        b.prop(s, "hull_miss_fraction")


def _menu_export(self, context):
    self.layout.operator("caliper.export_glb", text="Caliper XR (.glb)")


CLASSES = (CALIPER_Settings, CALIPER_PT_panel) + ops.CLASSES + ops_xr.CLASSES + export.CLASSES


def register():
    for c in CLASSES:
        bpy.utils.register_class(c)
    bpy.types.Scene.caliper = PointerProperty(type=CALIPER_Settings)
    bpy.types.TOPBAR_MT_file_export.append(_menu_export)


def unregister():
    bpy.types.TOPBAR_MT_file_export.remove(_menu_export)
    del bpy.types.Scene.caliper
    for c in reversed(CLASSES):
        bpy.utils.unregister_class(c)
