"""Shape-true retessellation of CAD-like parts.

Every part must come out lighter, and never farther than the tolerance from
the CAD - by the engine's own two-way measure, and again by Verify's
independent measure on the triangles that are drawn and exported.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bpy  # noqa: E402
from harness import bracket, check, done, housing, reset, shaft, triangles  # noqa: E402

from caliper_xr import build, cadxr, sampling, verify  # noqa: E402


class Settings:
    xr_tol_mm = 0.2
    xr_crease_deg = 30.0


reset()
parts = [bracket((0, 0, 0)), shaft((150, 0, 0)), housing((300, 0, 0))]
col = build.rebuilt_collection(bpy.context.scene)
for tol in (0.2, 0.1):
    Settings.xr_tol_mm = tol
    for ob in parts:
        bpy.context.view_layer.update()
        dg = bpy.context.evaluated_depsgraph_get()
        lpm = 1.0 / sampling.mm_per_local(ob, bpy.context.scene)
        bm, loopn, info, st = cadxr.run(ob, dg, Settings, lpm)
        rb = cadxr.make_object(bm, loopn, build.rebuild_name(ob), ob, col, True,
                               {"xr_source": ob.name, "xr_mode": "CAD_XR"})
        bpy.context.view_layer.update()
        dg = bpy.context.evaluated_depsgraph_get()
        src, new = triangles(ob), triangles(rb)
        f = verify.fit(ob, rb, dg, bpy.context.scene, exclude_inside=False)
        a, b = f["ref_to_rb"]["max"], f["rb_to_ref"]["max"]
        tag = f"{ob.name} @ {tol:g} mm"
        check(new < src, f"{tag}: {src} -> {new} triangles (x{src / new:.1f})")
        check(st["worst_mm"] <= tol + 1e-6, f"{tag}: engine, largest deviation {st['worst_mm']:.3f} mm")
        check(a <= tol + 1e-3 and b <= tol + 1e-3,
              f"{tag}: Verify, CAD->XR {a:.3f} mm, XR->CAD {b:.3f} mm")
        check(f["dims_err"] <= tol + 1e-3, f"{tag}: bounding box within {f['dims_err']:.3f} mm")
done("retess")
