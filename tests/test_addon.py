"""The add-on end to end, through its operators, as a user drives it:
register, Check Shell, Retessellate for XR with Outer shell only, Verify,
Audit, export a GLB, Housekeeping, unregister. The CAD is never modified.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bpy  # noqa: E402
import numpy as np  # noqa: E402
from harness import bracket, check, done, housing, reset, shaft, triangles  # noqa: E402

import caliper_xr  # noqa: E402
from caliper_xr import build, verify  # noqa: E402


def fingerprint(ob):
    co = np.empty(len(ob.data.vertices) * 3)
    ob.data.vertices.foreach_get("co", co)
    return (len(ob.data.polygons), float(np.abs(co).sum()), tuple(ob.matrix_world.col[3]))


reset()
caliper_xr.register()
try:
    parts = [bracket((0, 0, 0)), housing((150, 0, 0)), shaft((300, 0, 0))]
    before = {o.name: fingerprint(o) for o in parts}
    for o in bpy.context.view_layer.objects:
        o.select_set(o in parts)
    bpy.context.view_layer.objects.active = parts[0]
    s = bpy.context.scene.caliper
    s.xr_quality, s.xr_shell_only, s.xr_shell_gap_mm = 'BALANCED', True, 2.0

    check(bpy.ops.caliper.check_shell() == {'FINISHED'}, "Check Shell runs")
    check("Caliper_Report" in bpy.data.texts, "the report is written")
    check(bpy.ops.caliper.retessellate() == {'FINISHED'}, "Retessellate for XR runs")
    col = bpy.data.collections.get(build.REBUILT)
    rbs = [o for o in col.all_objects if o.get("xr_mode") == "CAD_XR"]
    check(len(rbs) == len(parts), f"{len(rbs)} rebuilds in {build.REBUILT}")
    for rb in rbs:
        src = bpy.data.objects[rb["xr_source"]]
        check(rb["xr_genus_ok"] and rb["xr_dev_max_mm"] <= 0.2 + 1e-6,
              f"{rb.name}: {triangles(src)} -> {triangles(rb)} triangles, "
              f"deviation {rb['xr_dev_max_mm']:.3f} mm")
    check(all(fingerprint(o) == before[o.name] for o in parts), "the CAD parts are untouched")
    # the shaft's 6 x 24 mm cross bore can hide from all 128 camera directions;
    # its faces must stay, or the closed part comes out with pinholes
    t = verify.topology(bpy.data.objects[build.rebuild_name(bpy.data.objects["shaft"])])
    check(t["non_manifold"] == 0, f"closed part stays closed ({t['non_manifold']} open edges)")

    for o in bpy.context.view_layer.objects:
        o.select_set(o in rbs)
    check(bpy.ops.caliper.verify() == {'FINISHED'}, "Verify runs")
    check(bpy.ops.caliper.audit() == {'FINISHED'}, "Audit runs")
    path = os.path.join(tempfile.mkdtemp(), "caliper_test.glb")
    check(bpy.ops.caliper.export_glb(filepath=path) == {'FINISHED'} and os.path.getsize(path) > 0,
          f"GLB exported ({os.path.getsize(path) // 1024 if os.path.exists(path) else 0} KB)")
    check(bpy.ops.caliper.housekeeping() == {'FINISHED'}, "Housekeeping runs")
    check(bpy.data.objects.get("Caliper_Shell_Check") is None, "the check preview is cleaned up")
finally:
    caliper_xr.unregister()
check(not hasattr(bpy.types.Scene, "caliper"), "unregisters cleanly")
done("addon")
