"""Read a CAD mesh as welded triangle arrays, with the signals that mark
where one CAD face ends and the next begins. Read-only on the source.

A tessellated CAD part carries its B-rep structure implicitly:
  * exporters often split vertices along face boundaries,
  * custom (split) normals are continuous inside a face and break across
    a face boundary,
  * material changes follow faces,
  * real creases show as large dihedral angles.
Every one of those is recorded per edge so segmentation can use them.
"""

import numpy as np


class CadMesh:
    """Welded triangle soup plus adjacency. Coordinates are object-local."""

    def __init__(self):
        self.V = None          # (n,3) welded vertex positions
        self.T = None          # (m,3) triangles, welded vertex ids
        self.N = None          # (m,3) unit face normals
        self.A = None          # (m,)  face areas
        self.mat = None        # (m,)  material index per triangle
        self.CN = None         # (m,3,3) corner normals per triangle corner
        self.E = None          # (k,2) unique edges (sorted welded ids)
        self.TE = None         # (m,3) edge id of tri edge (v0v1, v1v2, v2v0)
        self.ecount = None     # (k,) faces per edge
        self.ef = None         # (k,2) the two faces of a manifold edge, -1 else
        self.split = None      # (k,) edge was split in the source (B-rep hint)
        self.crease = None     # (k,) custom normals break across the edge
        self.UV = None         # (m,3,2) source UVs per corner, or None
        self.uvseam = None     # (k,) UVs break across the edge
        self.diag = 1.0


def load(ob, depsgraph):
    """Evaluated mesh of ob as a CadMesh. Never writes to ob."""
    ev = ob.evaluated_get(depsgraph)
    me = ev.to_mesh()
    try:
        nv = len(me.vertices)
        co = np.empty(nv * 3)
        me.vertices.foreach_get("co", co)
        V = co.reshape(-1, 3)
        lt = me.loop_triangles
        nt = len(lt)
        T = np.empty(nt * 3, np.int64)
        lt.foreach_get("vertices", T)
        T = T.reshape(-1, 3)
        L = np.empty(nt * 3, np.int64)
        lt.foreach_get("loops", L)
        L = L.reshape(-1, 3)
        P = np.empty(nt, np.int64)
        lt.foreach_get("polygon_index", P)
        nl = len(me.loops)
        cn = np.empty(nl * 3)
        me.corner_normals.foreach_get("vector", cn)
        CN = cn.reshape(-1, 3)
        mi = np.zeros(len(me.polygons), np.int64)
        me.polygons.foreach_get("material_index", mi)
        UVc = None
        uvl = me.uv_layers.active
        if uvl is not None and len(uvl.data):
            uv = np.empty(nl * 2)
            uvl.data.foreach_get("uv", uv)
            UVc = uv.reshape(-1, 2)
        custom = bool(me.has_custom_normals)
        sharp_attr = me.attributes.get("sharp_edge")
        has_sharp = False
        if sharp_attr is not None and sharp_attr.domain == 'EDGE':
            sh = np.zeros(len(me.edges), bool)
            sharp_attr.data.foreach_get("value", sh)
            has_sharp = bool(sh.any())
        sf = me.attributes.get("sharp_face")
        all_flat = False
        if sf is not None and sf.domain == 'FACE':
            fl = np.zeros(len(me.polygons), bool)
            sf.data.foreach_get("value", fl)
            all_flat = bool(fl.all())
    finally:
        ev.to_mesh_clear()
    cm = build(V, T, CN[L], mi[P], use_normal_breaks=(custom or has_sharp) and not all_flat,
               TUV=None if UVc is None else UVc[L])
    return cm


def build(V, T, TCN, tmat, use_normal_breaks=True, TUV=None):
    cm = CadMesh()
    if len(T) == 0:
        raise ValueError("mesh has no faces")
    diag = float(np.linalg.norm(V.max(0) - V.min(0))) or 1.0
    cm.diag = diag

    # Raw half-edge keys, before welding: an exporter that split vertices
    # along a CAD face boundary leaves those faces unconnected here.
    raw = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    raw = np.sort(raw, axis=1)

    # Weld exact duplicates (CAD exporters duplicate bit-for-bit).
    eps = max(diag * 1e-8, 1e-12)
    q = np.round(V / eps).astype(np.int64)
    _u, first, inv = np.unique(q, axis=0, return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    Vw = V[first]
    Tw = inv[T]

    # Drop triangles that collapsed or have no area.
    ok = (Tw[:, 0] != Tw[:, 1]) & (Tw[:, 1] != Tw[:, 2]) & (Tw[:, 2] != Tw[:, 0])
    a, b, c = Vw[Tw[:, 0]], Vw[Tw[:, 1]], Vw[Tw[:, 2]]
    cr = np.cross(b - a, c - a)
    area2 = np.linalg.norm(cr, axis=1)
    ok &= area2 > (diag * 1e-9) ** 2
    m0 = len(Tw)
    keep_rows = np.concatenate([ok, ok, ok])
    Tw, TCN, tmat, cr, area2 = Tw[ok], TCN[ok], tmat[ok], cr[ok], area2[ok]
    TUV = TUV[ok] if TUV is not None else None
    src_index = np.flatnonzero(ok)
    raw = raw[keep_rows]
    m = len(Tw)

    cm.V, cm.T = Vw, Tw
    cm.src = src_index             # source loop-triangle index of each triangle
    cm.A = 0.5 * area2
    cm.N = cr / area2[:, None]
    cm.mat = tmat
    cm.CN = TCN
    cm.UV = TUV

    he = np.concatenate([Tw[:, [0, 1]], Tw[:, [1, 2]], Tw[:, [2, 0]]])
    hs = np.sort(he, axis=1)
    E, einv = np.unique(hs, axis=0, return_inverse=True)
    einv = einv.reshape(-1)
    cm.E = E
    cm.TE = einv.reshape(3, m).T
    cm.ecount = np.bincount(einv, minlength=len(E))

    # The two faces of each manifold edge, and which tri-edge slot they use.
    order = np.argsort(einv, kind="stable")
    starts = np.searchsorted(einv[order], np.arange(len(E)))
    ef = np.full((len(E), 2), -1, np.int64)
    es = np.full((len(E), 2), -1, np.int64)
    two = np.flatnonzero(cm.ecount == 2)
    r0 = order[starts[two]]
    r1 = order[starts[two] + 1]
    ef[two, 0], ef[two, 1] = r0 % m, r1 % m
    es[two, 0], es[two, 1] = r0 // m, r1 // m
    cm.ef = ef

    split = np.zeros(len(E), bool)
    split[two] = np.any(raw[r0] != raw[r1], axis=1)
    cm.split = split

    crease = np.zeros(len(E), bool)
    if use_normal_breaks and len(two):
        # corner normals of both faces at both edge ends
        f0, f1 = ef[two, 0], ef[two, 1]
        s0, s1 = es[two, 0], es[two, 1]
        va = Tw[f0, s0]
        vb = Tw[f0, (s0 + 1) % 3]
        n0a, n0b = TCN[f0, s0], TCN[f0, (s0 + 1) % 3]
        # in f1 the same vertices sit at s1 / s1+1, in either order
        w0, w1 = Tw[f1, s1], Tw[f1, (s1 + 1) % 3]
        m0a = np.where((w0 == va)[:, None], TCN[f1, s1], TCN[f1, (s1 + 1) % 3])
        m0b = np.where((w1 == vb)[:, None], TCN[f1, (s1 + 1) % 3], TCN[f1, s1])
        ca = np.einsum("ij,ij->i", _unit(n0a), _unit(m0a))
        cb = np.einsum("ij,ij->i", _unit(n0b), _unit(m0b))
        crease[two] = np.minimum(ca, cb) < np.cos(np.radians(2.0))
    cm.crease = crease
    uvseam = np.zeros(len(E), bool)
    if TUV is not None and len(two):
        f0, f1 = ef[two, 0], ef[two, 1]
        s0, s1 = es[two, 0], es[two, 1]
        va, vb = Tw[f0, s0], Tw[f0, (s0 + 1) % 3]
        u0a, u0b = TUV[f0, s0], TUV[f0, (s0 + 1) % 3]
        w0 = Tw[f1, s1]
        u1a = np.where((w0 == va)[:, None], TUV[f1, s1], TUV[f1, (s1 + 1) % 3])
        u1b = np.where((w0 == va)[:, None], TUV[f1, (s1 + 1) % 3], TUV[f1, s1])
        d = np.maximum(np.abs(u0a - u1a).max(1), np.abs(u0b - u1b).max(1))
        uvseam[two] = d > 1e-5
        if uvseam[two].mean() > 0.3:
            # per-face or auto-generated UVs: seams everywhere mean nothing,
            # and honouring them would stop every face from simplifying
            uvseam[:] = False
            cm.UV = None
    cm.uvseam = uvseam
    cm.dropped = m0 - m
    return cm


def _unit(x):
    n = np.linalg.norm(x, axis=1)
    n[n == 0] = 1.0
    return x / n[:, None]


def dihedral_cos(cm):
    """cos of the angle between the two face normals of each manifold edge."""
    c = np.ones(len(cm.E))
    two = cm.ecount == 2
    f0, f1 = cm.ef[two, 0], cm.ef[two, 1]
    c[two] = np.einsum("ij,ij->i", cm.N[f0], cm.N[f1])
    return c


def subset(cm, keepf):
    """A CadMesh of only the kept triangles, keeping the source's split and
    crease signals (they are properties of edges, found again by key)."""
    key_old = cm.E[:, 0] * (len(cm.V) + 1) + cm.E[:, 1]
    T = cm.T[keepf]
    out = CadMesh()
    out.V, out.T, out.diag = cm.V, T, cm.diag
    out.src = cm.src[keepf]
    out.N, out.A, out.mat, out.CN = cm.N[keepf], cm.A[keepf], cm.mat[keepf], cm.CN[keepf]
    out.UV = cm.UV[keepf] if getattr(cm, "UV", None) is not None else None
    m = len(T)
    he = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    E, einv = np.unique(np.sort(he, axis=1), axis=0, return_inverse=True)
    einv = einv.reshape(-1)
    out.E, out.TE = E, einv.reshape(3, m).T
    out.ecount = np.bincount(einv, minlength=len(E))
    order = np.argsort(einv, kind="stable")
    starts = np.searchsorted(einv[order], np.arange(len(E)))
    ef = np.full((len(E), 2), -1, np.int64)
    two = np.flatnonzero(out.ecount == 2)
    ef[two, 0] = order[starts[two]] % m
    ef[two, 1] = order[starts[two] + 1] % m
    out.ef = ef
    key_new = E[:, 0] * (len(cm.V) + 1) + E[:, 1]
    srt = np.argsort(key_old)
    pos = np.searchsorted(key_old[srt], key_new)
    pos = np.clip(pos, 0, len(srt) - 1)
    found = key_old[srt][pos] == key_new
    idx = srt[pos]
    out.split = np.where(found, cm.split[idx], False)
    out.crease = np.where(found, cm.crease[idx], False)
    out.uvseam = np.where(found, cm.uvseam[idx], False)
    out.dropped = 0
    return out
