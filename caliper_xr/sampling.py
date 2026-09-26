"""Read-only access to the reference CAD part: frame, rays, radius grid.

Nothing here writes to the reference object. Geometry is read from the
evaluated depsgraph copy and ray-cast against a private BVH in the object's
LOCAL space, so neighbouring parts can never be hit and the rebuild can take
the reference's matrix_world unchanged.
"""

import math

import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree

from . import geom


def mm_per_local(ob, scene):
    """How many millimetres one local unit of this object is."""
    s = ob.matrix_world.to_scale()
    mean = (abs(s.x) + abs(s.y) + abs(s.z)) / 3.0
    return scene.unit_settings.scale_length * 1000.0 * mean


def read_reference(ob, depsgraph):
    """Evaluated verts (N,3 local), polygon count and a BVH. No writes."""
    ev = ob.evaluated_get(depsgraph)
    me = ev.to_mesh()
    try:
        n = len(me.vertices)
        co = np.empty(n * 3, dtype=np.float64)
        me.vertices.foreach_get("co", co)
        polys = len(me.polygons)
    finally:
        ev.to_mesh_clear()
    bvh = BVHTree.FromObject(ob, depsgraph)
    return co.reshape(-1, 3), polys, bvh


def component_count(ob):
    """Connected components of the reference mesh (walks link edges)."""
    me = ob.data
    n = len(me.vertices)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for e in me.edges:
        a, b = find(e.vertices[0]), find(e.vertices[1])
        if a != b:
            parent[a] = b
    return len({find(i) for i in range(n)})


# --------------------------------------------------------------------------
# Axis frame
# --------------------------------------------------------------------------

_LOCAL = np.eye(3)


def _snap(a, deg=5.0):
    c = math.cos(math.radians(deg))
    for ax in _LOCAL:
        d = float(np.dot(a, ax))
        if abs(d) >= c:
            return ax * np.sign(d)
    return a / np.linalg.norm(a)


def frame(V, mode, cursor=None):
    """Return (C, a, u, v): a point on the axis, the axis, and a basis.

    AUTO picks the principal direction whose other two variances are most
    alike - the symmetry axis of a revolved part, whether long (shaft) or
    flat (nut) - and snaps it to a local axis within 5 deg, since CAD parts
    usually import axis-aligned.
    """
    if mode == 'CURSOR' and cursor is not None:
        C, a = cursor
        a = np.asarray(a, float)
        a = a / np.linalg.norm(a)
        C = np.asarray(C, float)
    elif mode in {'X', 'Y', 'Z'}:
        a = _LOCAL['XYZ'.index(mode)].copy()
        C = None
    else:
        X = V - V.mean(axis=0)
        w, E = np.linalg.eigh(np.cov(X.T) if len(V) > 3 else np.eye(3))
        best, score = 2, float("inf")
        for i in range(3):
            j, k = [x for x in range(3) if x != i]
            s = abs(w[j] - w[k]) / max(w[j], w[k], 1e-18)
            if s < score:
                best, score = i, s
        a = _snap(E[:, best])
        C = None

    ref = _LOCAL[int(np.argmin(np.abs(_LOCAL @ a)))]
    u = ref - a * np.dot(ref, a)
    u /= np.linalg.norm(u)
    v = np.cross(a, u)

    if C is None:
        pu, pv = V @ u, V @ v
        C = u * (pu.min() + pu.max()) * 0.5 + v * (pv.min() + pv.max()) * 0.5
    else:
        C = C - a * np.dot(C, a)
    return C, a, u, v


# --------------------------------------------------------------------------
# Rays
# --------------------------------------------------------------------------

class RingSampler:
    """Casts rays inward toward the axis and records the first hit radius."""

    def __init__(self, bvh, C, a, u, v, RM, thetas):
        self.bvh = bvh
        self.C, self.a = Vector(C), Vector(a)
        self.RM = RM
        self.dirs = [Vector(math.cos(t) * u + math.sin(t) * v) for t in thetas]
        self.rays = 0

    def ring(self, z):
        base = self.C + self.a * z
        RM = self.RM
        out = np.full(len(self.dirs), np.nan)
        cast = self.bvh.ray_cast
        for j, d in enumerate(self.dirs):
            # No normal filter: CAD shells often import with inverted normals,
            # and a facing test silently rejects real geometry.
            loc, _n, _i, dist = cast(base + d * RM, -d, RM * 2.0)
            if loc is not None:
                r = RM - dist
                if r > 0:
                    out[j] = r
        self.rays += len(self.dirs)
        return out


def _extend_leg(R, raw, RM_mm, s):
    """On a leg span, carry the clean profile straight through the junction.

    Near the junction the rays see the other leg: they start inside its
    material, pass down its bore, or miss. Those rows describe the other leg,
    not this one, so they are replaced (in place) by the last clean row -
    the leg is extended past the junction as a capped solid, and the corner
    sphere plus the other leg cover what it does not. Returns rows replaced.
    """
    if not (s.axis == 'CURSOR' or s.range_from > 0.0 or s.range_to < 1.0):
        return 0
    S = len(R)
    if s.range_from > 0.0:
        order = range(S)                    # junction at the low end
    elif s.range_to < 1.0:
        order = range(S - 1, -1, -1)
    else:
        return 0
    free = R[-max(1, S // 10):] if s.range_from > 0.0 else R[:max(1, S // 10)]
    spread_ref = float(np.median(free.max(axis=1) - free.min(axis=1)))
    limit = 3.0 * spread_ref + 4.0 * s.angle_tol_mm + s.step_mm
    bad = []
    for k in order:
        spread = float(R[k].max() - R[k].min())
        dirty = (np.isnan(raw[k]).any() or spread > limit
                 or np.nanmax(raw[k]) >= 0.97 * RM_mm)
        if not dirty:
            break
        bad.append(k)
    if bad and len(bad) < S:
        clean = bad[-1] + (1 if s.range_from > 0.0 else -1)
        R[bad] = R[clean]
    return len(bad)


def sample_grid(bvh, V, frame_, lpm, s):
    """Sample the outer surface as a (station x angle) radius grid in mm.

    Returns a dict with z (mm along the axis), rows (hole-filled radii, mm),
    raw rows, dz, thetas, forced step pairs and the filled-cell count.
    """
    C, a, u, v = frame_
    mm = 1.0 / lpm                       # mm per local unit
    zv = (V - C) @ a
    z0, z1 = float(zv.min()), float(zv.max())
    L = z1 - z0
    lo = z0 + L * s.range_from
    hi = z0 + L * s.range_to
    if hi - lo <= 1e-9:
        raise ValueError("Axis range is empty")

    # Centre and ray start come from the geometry inside the sampled span
    # only - on an elbow leg, the other leg must not drag the axis sideways.
    inside = (zv >= lo - 1e-9) & (zv <= hi + 1e-9)
    Vs = V[inside] if inside.sum() >= 3 else V
    if s.axis != 'CURSOR':
        # A leg's free end is a clean section of that leg alone; its
        # junction end is not. Centre on the free end's 10 % band.
        band = 0.1 * (hi - lo)
        if s.range_from > 0.0 and s.range_to >= 1.0:
            sel = zv >= hi - band
        elif s.range_to < 1.0 and s.range_from <= 0.0:
            sel = zv <= lo + band
        else:
            sel = inside
        Vc = V[sel] if sel.sum() >= 3 else Vs
        pu, pv = Vc @ u, Vc @ v
        C = u * (pu.min() + pu.max()) * 0.5 + v * (pv.min() + pv.max()) * 0.5
        frame_ = (C, a, u, v)
    zs_in = (Vs - C) @ a
    rv = np.linalg.norm((Vs - C) - np.outer(zs_in, a), axis=1)
    if s.ray_radius_mm > 0:
        RM = s.ray_radius_mm * lpm            # user-capped: just outside the leg
    else:
        RM = float(rv.max()) * 1.05 + 1.0 * lpm   # just outside the part

    A = int(s.angular_samples)
    thetas = np.arange(A) * (2.0 * math.pi / A)
    rs = RingSampler(bvh, C, a, u, v, RM, thetas)

    step = max(s.sample_step_mm * lpm, (hi - lo) / float(s.max_stations))
    n = max(2, int(math.ceil((hi - lo) / step)) + 1)
    eps = min(step * 0.25, 0.01 * lpm)
    zs = np.linspace(lo + eps, hi - eps, n)
    raw = np.array([rs.ring(z) for z in zs]) * mm

    dz_mm = (zs[1] - zs[0]) * mm if n > 1 else 1.0
    miss = float(np.mean(~np.isfinite(raw)))
    R, hole_mask = geom.fill_holes(raw, thetas[1] - thetas[0], dz_mm,
                                s.max_hole_mm, s.max_hole_deg, s.fill_mode)

    if not np.isfinite(R).any():
        raise ValueError("No ray hit the part on this axis")
    R = np.where(np.isfinite(R), R, 0.0)
    extended = _extend_leg(R, raw, RM * mm, s)

    # Hard steps: locate by bisection, then place two rings 0.01 mm apart,
    # each carrying the filled values of the plateau on its side.
    zs_out = zs.copy()
    forced = []
    jump = np.abs(np.diff(R, axis=0)).max(axis=1) if n > 1 else np.array([])
    for k in np.flatnonzero(jump > s.step_mm):
        za, zb = zs[k], zs[k + 1]
        ra, rb = raw[k], raw[k + 1]
        for _ in range(8):
            zm = 0.5 * (za + zb)
            rm = rs.ring(zm) * mm
            da = np.nanmedian(np.abs(rm - ra)) if np.isfinite(rm - ra).any() else np.inf
            db = np.nanmedian(np.abs(rm - rb)) if np.isfinite(rm - rb).any() else np.inf
            if da <= db:
                za = zm
            else:
                zb = zm
        zstep = 0.5 * (za + zb)
        half = min(0.01 * lpm, (zs[k + 1] - zs[k]) * 0.25)
        zs_out[k] = max(zstep - half, zs_out[k])
        zs_out[k + 1] = zstep + half
        forced += [k, k + 1]

    # First and last rings sit exactly on the true extremes.
    zs_out[0], zs_out[-1] = lo, hi
    return dict(z=zs_out * mm, rows=R, raw=raw, thetas=thetas, extended=extended,
                forced=forced, holes=geom.count_regions(hole_mask),
                hole_cells=int(hole_mask.sum()), miss=miss, rays=rs.rays,
                frame=frame_, lpm=lpm, rmax_mm=float(np.nanmax(R)))
