# Changelog

## 1.0.0 - first public release

**Caliper XR** - heavy CAD in, lean XR out, measured both ways.

- **Retessellate for XR.** CAD faces are recovered from the tessellation and rebuilt lean on the
  original surface: n-gons on flat faces, quads on strips, crisp edges and cut-outs kept, the
  CAD's normals, UVs and materials carried over. Quality classes 0.1 / 0.2 / 0.5 mm or custom.
- **Two-way tolerance guarantee.** Every CAD point within tolerance of the result and every
  result face within tolerance of the CAD, or the part is flagged and export refuses it.
- **Exact distances.** Blender's BVH is single precision and misjudges sheet-metal slivers by up
  to 0.67 mm; Caliper finds candidates with it and measures in double precision.
- **On-surface refinement.** Where no CAD vertex is available, points are placed exactly on the
  CAD surface, with interpolated normals and UVs.
- **Shortcuts stay on the surface.** Outline simplification never cuts an arc across the hole it
  runs around.
- **Remove hidden geometry.** GPU visibility from 128 directions; a hidden face surrounded by
  visible ones is kept, so a narrow bore never ends up with holes.
- **Outer shell only.** Voxel reachability with gap sealing: inner bodies, closed cavities,
  faces between touching bodies and parts behind narrow seams are left out; per selection or
  per part.
- **Check Shell.** Shows what would be left out - report and red preview - without building.
- **Pinned quad diagonals.** Quads are split on export along the diagonal that was measured.
- **Verify, Audit, GLB export** with Draco, pre-flight report and export gate.
- **Headless.** GPU visibility runs from the command line in Blender 5.
- 4 headless test suites, 55 checks, on procedural CAD-like parts.
