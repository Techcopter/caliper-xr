# CALIPER XR · OPERATOR'S MANUAL

```
DOC  CXR-OM-1.0          UNIT  mm          TOLERANCE  two-way          SOURCE CAD  read-only
```

| § | Section |
|---|---|
| 1 | [Principles](#1--principles) |
| 2 | [Installation](#2--installation) |
| 3 | [Preparing the scene](#3--preparing-the-scene) |
| 4 | [Retessellate for XR](#4--retessellate-for-xr) |
| 5 | [Remove hidden geometry](#5--remove-hidden-geometry) |
| 6 | [Outer shell only and Check Shell](#6--outer-shell-only-and-check-shell) |
| 7 | [Reading the report](#7--reading-the-report) |
| 8 | [Verify, Audit, originals](#8--verify-audit-originals) |
| 9 | [Export](#9--export) |
| 10 | [Batch and command line](#10--batch-and-command-line) |
| 11 | [Legacy: outer-shell section loft](#11--legacy-outer-shell-section-loft) |
| 12 | [Troubleshooting](#12--troubleshooting) |

---

## 1 · Principles

Caliper XR turns tessellated CAD into lean meshes for VR, AR and XR under four rules:

1. **The source CAD is read-only.** Every result is a new object in the `REBUILT` collection,
   named `<part>_REBUILT`, with the source's transform. The original is never edited, moved,
   decimated or deleted.
2. **Every vertex lies on the CAD surface.** Faces are rebuilt, not decimated: a vertex is either
   one of the CAD's own or a point placed exactly on a CAD triangle.
3. **The tolerance is a two-way hard limit.** Every CAD vertex and face centre lies within the
   tolerance of the result, and every result face centre and edge midpoint lies within the
   tolerance of the CAD. A part that cannot meet it is flagged, never passed quietly.
4. **Nothing is removed out of sight.** *Check Shell* shows what the hidden-geometry and
   outer-shell tests leave out, before a build.

## 2 · Installation

| Step | Action |
|---|---|
| 1 | Download `caliper_xr-<version>.zip` from the repository's **Releases** page |
| 2 | Blender ▸ **Edit ▸ Preferences ▸ Get Extensions** ▸ **⌄** (top right) ▸ **Install from Disk…** ▸ the zip |
| 3 | 3D Viewport ▸ `N` ▸ **Caliper** tab |

**Update:** install the newer zip the same way; it replaces the old version.
**Remove:** Preferences ▸ Get Extensions ▸ Caliper XR ▸ **⌄** ▸ **Uninstall**.
**Requirements:** Blender 4.2 LTS or newer. Developed and tested on Blender 5.2, macOS.

## 3 · Preparing the scene

- **Scale.** Tolerances are in millimetres. Caliper converts with the scene's
  *Unit Scale* (`Scene ▸ Units ▸ Unit Scale`) and each object's scale: 1 unit = 1 m at scale 1.0
  means 0.2 mm is 0.0002 units. Import CAD at its real size; if it arrives 1000× too large,
  fix the import scale rather than the tolerance.
- **One object per part** gives the best result - an assembly imported part by part. A whole
  device merged into one mesh works too; it is treated as one part.
- **Modifiers** on the CAD are evaluated, not applied: Caliper reads the evaluated mesh.
- **Normals and UVs.** The CAD's custom normals, UV seams, sharp edges and material boundaries
  are read as clues to where one CAD face ends and the next begins. Keep them from the import.
- **Visibility.** Only *visible, selected* mesh objects are processed. Objects in `REBUILT` and
  Caliper's own previews are skipped.

## 4 · Retessellate for XR

Select the parts and press **Retessellate for XR**.

| Control | Default | Meaning |
|---|---|---|
| **Quality** | Balanced | *High* 0.1 mm - close inspection. *Balanced* 0.2 mm - indistinguishable at arm's length. *Light* 0.5 mm - room-scale scenes, mobile headsets. *Custom* - set **Tolerance (mm)** yourself |
| **Crease angle** | 30° | Edges sharper than this always stay crisp edges, even where the CAD's normals are smooth |
| **Keep materials** | on | The rebuild takes each part's effective materials, including ones linked to the object |

What happens, per part:

1. **CAD faces are recovered** from the tessellation - by exporter splits, normal creases, UV
   seams, material changes and dihedral angle - and flat regions are extracted as planes.
2. **Each face is flattened** in 2D: flat faces by their plane, curved faces by a whole-face
   unwrap, faces that will not unwrap by projection pieces.
3. **Outlines are simplified once, globally**, so two faces sharing a border share every vertex
   on it - no cracks. A shortcut is accepted only if it stays within the tolerance of the CAD
   surface, so an arc around a hole is never cut straight across it.
4. **Each face is rebuilt** as a constrained Delaunay triangulation of its outline and refined
   only where the result strays: with CAD vertices, or points placed on the CAD surface. A face
   never comes out heavier than its source; if it cannot be rebuilt, its source triangles are
   kept and its neighbours are rebuilt against it.
5. **Both directions are measured** in double precision. Anything over the tolerance gets its
   dropped outline vertices back; a face still over goes back to its exact source triangles.
6. **The result is assembled**: n-gons on planes, quads on strips, the CAD's normals and UVs,
   every quad split along the diagonal that was measured - so what is exported is what was
   checked.

The rebuild's panel line shows **Faces** and **Largest deviation**, with a check mark when it is
within the tolerance. Custom properties on the object keep the record: `xr_source`, `xr_mode`,
`xr_tol_mm`, `xr_dev_max_mm`, `xr_genus_ok`, `xr_built`.

## 5 · Remove hidden geometry

**On by default.** The selection is drawn on the GPU from 128 directions around it, every
triangle in its own colour, both sides, no culling. A CAD face that never shows up in any view is
left out.

- The test sees through vents, gaps and grazing angles down to about a pixel.
- It always judges the **whole selection**: a face covered by another part counts as hidden.
- A hidden face whose neighbouring faces are all visible is **kept** - 128 directions can miss
  the inside of a narrow bore, and dropping it would open a hole in a visible surface.
- A part with nothing visible is skipped entirely: *hidden part skipped* in the report.
- Without a GPU (rare), ray sampling runs instead and the report says it is less certain.
  From the command line, Blender 5 initialises the GPU automatically.

Turn it off when parts will move, open or come apart in XR: a face hidden in the assembled state
may be exposed later.

## 6 · Outer shell only and Check Shell

**Outer shell only** builds just the outer skin. The geometry is sampled onto a voxel grid,
openings narrower than the gap are sealed, and the air reachable from outside is flood-filled.
A face stays when it borders that outside air; internal faces go:

- bodies inside other bodies,
- the inner walls of closed boxes and cavities,
- faces between touching bodies,
- anything reachable only through a seam or vent narrower than the gap.

| Control | Default | Meaning |
|---|---|---|
| **Seal gaps under (mm)** | 2 | Openings narrower than this are not a way in |
| **Shell of** | Whole selection | *Whole selection* - one skin around all selected parts. *Each part* - every part judged on its own, keeping its complete skin |

**Choosing the gap.** Reachability is all-or-nothing per cavity: one opening wider than the gap
makes a whole interior count as outside. With a small gap, a device with vents keeps its
interior; a gap wider than its openings removes the interior entirely and leaves the skin. On a
sheet-metal cabinet of 820 parts, a 2 mm gap changed little, while a 25 mm gap took it from
54,979 to about 35,000 triangles with no visible change from outside.

**Resolution.** The grid is capped near 640 voxels along the longest side and 24 million
in total. On a large assembly the smallest sealed opening is set by that cap, not by the gap
setting - the report states the real figure (*openings under about 5.6 mm sealed*).

**Combined with the camera test.** A face is left out when *either* test says so. For *Each
part*, turn off *Remove hidden geometry* too - the camera test always judges the whole selection.

**Check Shell** runs both tests without building anything:

- the report lists every part that would lose faces - how many, and whether they are hidden or
  inside the shell - and the total kept;
- a red object, `Caliper_Shell_Check`, in the `CALIPER_CHECK` collection, holds exactly the
  faces that would go. They sit inside the parts, so hide the parts (`H`) or isolate the preview
  (numpad `/`) to look. It never renders and never exports; **Housekeeping** deletes it.

Check Shell always runs the shell test, whether or not *Outer shell only* is ticked, so you can
see what it would do before switching it on.

## 7 · Reading the report

Every run appends to the text block **`Caliper_Report`** (Text Editor) and prints to the console.

| Line | Meaning |
|---|---|
| `# CAD -> XR <date> tolerance 0.200 mm, N part(s)` | Run header |
| `visibility: GPU, 128 views of N parts, 2.1s` | The camera test ran on the GPU |
| `outer shell: selection as one, voxels 2.8 mm, openings under about 5.6 mm sealed` | The shell test, its grid and what it really sealed |
| `== part -> part_REBUILT: 1492 -> 442 triangles` | Source and result triangles |
| `source 1492 triangles; 3 CAD faces recovered` | What segmentation found |
| `left out 302 CAD faces (1085 triangles): 299 hidden, 3 inside the outer shell` | What the two tests removed |
| `largest deviation 0.180 mm (tolerance 0.200): CAD to XR 0.180, XR to CAD 0.172` | The two-way measure - exact |
| `6 flat faces, 1 curved faces unwrapped whole, 8 projection pieces` | How faces were flattened |
| `rebuilt 15 pieces (2 kept as in the source) in 0.2s` | Pieces that could not beat their source keep it |
| `(deviation above tolerance: check it)` | The part missed the tolerance: it is flagged and blocked from export |
| `== part: not on the outside, skipped` | Nothing of the part is visible or on the shell |
| `TOTAL 6964 -> 1254 triangles (x5.6), 0 hidden part(s) skipped, 0 failed, 4s` | Run summary |

## 8 · Verify, Audit, originals

- **Verify** re-measures the selected rebuilds (or the rebuilds of selected parts) against their
  CAD: topology, bounding box, matrix, and deviation both ways on the triangles that are drawn and
  exported. A CAD → XR rebuild is judged by its own two-way record, which knows which faces were
  hidden and removed.
- **Audit** lists every rebuild in `REBUILT` that needs redoing, and changes nothing.
- **Show / Hide Originals** toggles the CAD parts of the rebuilds in view, to compare.
- **Housekeeping** deletes Caliper's temporary objects (the Check Shell preview) and orphan
  meshes.

## 9 · Export

**Export XR GLB** (also *File ▸ Export ▸ Caliper XR (.glb)*) writes the rebuilds as one GLB.

| Option | Default | Meaning |
|---|---|---|
| **Objects** | All REBUILT | Or only the selected rebuilds |
| **Draco compression** | on (level 6) | Much smaller files; the viewer needs a Draco decoder |
| **Merge into one mesh** | off | One draw call, but parts can no longer be picked or highlighted |
| **+Y up** | on | glTF convention |
| **Block unverified / failing parts** | on | Refuses to export a part that is flagged |
| **Write report beside GLB** | on | `<name>.report.txt` next to the file |

The pre-flight check warns above the **Scene triangle budget** (default 300,000, a mobile-headset
class figure).

## 10 · Batch and command line

Everything the panel does is an operator, so it runs headless:

```python
# batch.py - blender -b assembly.blend --python batch.py
import bpy

bpy.ops.preferences.addon_enable(module="bl_ext.user_default.caliper_xr")
for ob in bpy.context.view_layer.objects:
    ob.select_set(ob.type == 'MESH' and ob.visible_get())

s = bpy.context.scene.caliper
s.xr_quality = 'BALANCED'          # HIGH / BALANCED / LIGHT / CUSTOM (+ s.xr_tol_mm)
s.xr_shell_only = False            # True, with s.xr_shell_gap_mm and s.xr_shell_scope

bpy.ops.caliper.check_shell()      # optional: report what would be left out
bpy.ops.caliper.retessellate()
bpy.ops.caliper.export_glb(filepath="/tmp/assembly.glb")
print(bpy.data.texts["Caliper_Report"].as_string())
```

Operators: `caliper.retessellate`, `caliper.check_shell`, `caliper.toggle_source`,
`caliper.verify`, `caliper.audit`, `caliper.export_glb`, `caliper.housekeeping`.

## 11 · Legacy: outer-shell section loft

Collapsed at the bottom of the panel. It rebuilds **revolved parts** - fittings, nuts, elbows -
as closed, single-shell meshes by ray-sampling rings along an axis and lofting them, capping
holes and pockets up to a set size, within a face budget per part class. Use it when a closed,
genus-0 shell matters more than matching every CAD face; for everything else, Retessellate for XR
is the method.

## 12 · Troubleshooting

| Symptom | Cause and remedy |
|---|---|
| *Select one or more visible CAD parts* | Nothing selected, the parts are hidden, or they are rebuilds |
| Deviations look 1000× too large or too small | Wrong scale: check *Unit Scale* and the object scale (§3) |
| A part is flagged *deviation above tolerance* | Rare - the report says by how much. Rerun at a slightly larger tolerance, or keep the part as is |
| A result is heavier than expected | Curved faces that do not unwrap become projection pieces; faces that cannot beat their source keep it. Try a larger tolerance or crease angle |
| Something visible was removed | Turn off *Remove hidden geometry*, or run *Check Shell* to see what it drops |
| *Outer shell only* removed too little | The interior is reachable through an opening wider than the gap - raise the gap |
| *Outer shell only* removed too much | Lower the gap, or use *Each part* |
| Export is blocked | A part is flagged or unverified: run *Verify*, fix or rebuild it, or untick the block |
| It is slow | Time goes to the shell grid and the two-way measure. A 190k-triangle single mesh takes about 2 minutes; an 820-part assembly about 25 s |

---

<sub>CALIPER XR · GPL-3.0-or-later · maintained by Techcopter</sub>
