"""Split a CadMesh into recovered CAD faces, then into height-field pieces,
and simplify the curve network between them once, globally.

  patches      smooth regions bounded by creases, CAD face boundaries,
               material changes, open and non-manifold edges.
  sub-patches  each patch cut so all its normals lie within MAX_DEV of one
               direction: it then projects one-to-one onto a plane, which is
               where it is re-triangulated.
  chains       the edges between sub-patches, walked corner to corner and
               reduced with Douglas-Peucker. Both sides of a chain use the
               same kept vertices, so the result has no cracks.
"""

import math

import numpy as np

from .cadmesh import dihedral_cos

MAX_DEV = math.radians(70.0)


class _UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.p[b] = a


def _labels(n, pairs):
    uf = _UF(n)
    for a, b in pairs:
        uf.union(int(a), int(b))
    roots = np.array([uf.find(i) for i in range(n)])
    _u, lab = np.unique(roots, return_inverse=True)
    return lab.reshape(-1)


def hard_edges(cm, crease_deg):
    """Edges no patch may cross."""
    two = cm.ecount == 2
    hard = ~two
    same_mat = np.ones(len(cm.E), bool)
    same_mat[two] = cm.mat[cm.ef[two, 0]] == cm.mat[cm.ef[two, 1]]
    hard |= ~same_mat
    hard |= cm.split | cm.crease
    if getattr(cm, "uvseam", None) is not None:
        hard |= cm.uvseam                          # a texture must not smear across a seam
    hard |= dihedral_cos(cm) < math.cos(math.radians(crease_deg))
    return hard


def patches(cm, hard):
    soft = np.flatnonzero(~hard)
    return _labels(len(cm.T), cm.ef[soft])


def extract_planes(cm, patch, hard, dist_tol, min_area):
    """Split flat CAD faces out of smooth patches.

    Custom normals often flow smoothly across tangent fillets, so a whole
    housing skin arrives as one smooth patch. Its flat faces are recovered as
    coplanar regions that are much wider than the facet strips beside them:
    a strip of a curved surface is about as wide as its neighbouring strips,
    a real flat face is not. Returns a new patch id per face.
    """
    m = len(cm.T)
    two = cm.ecount == 2
    dc = dihedral_cos(cm)
    flatj = two & ~hard & (dc > math.cos(math.radians(0.25)))
    reg = _labels(m, cm.ef[flatj])
    nreg = int(reg.max()) + 1
    area = np.bincount(reg, weights=cm.A, minlength=nreg)
    width = np.sqrt(area)                                  # rough, for small regions
    order = np.argsort(reg, kind="stable")
    bounds = np.searchsorted(reg[order], np.arange(nreg + 1))
    coplanar = np.zeros(nreg, bool)
    for r in np.flatnonzero(area >= min_area):
        f = order[bounds[r]:bounds[r + 1]]
        n = (cm.N[f] * cm.A[f, None]).sum(0)
        n /= max(np.linalg.norm(n), 1e-30)
        P = cm.V[np.unique(cm.T[f])]
        c = P.mean(0)
        coplanar[r] = np.abs((P - c) @ n).max() <= dist_tol
        Q = P - c
        Q -= np.outer(Q @ n, n)
        if len(Q) >= 2:
            _u, _s, vt = np.linalg.svd(Q, full_matrices=False)
            length = float(np.ptp(Q @ vt[0]))
            width[r] = area[r] / max(length, 1e-30)
    curvj = two & ~hard & ~flatj
    ra, rb = reg[cm.ef[curvj, 0]], reg[cm.ef[curvj, 1]]
    minnb = np.full(nreg, np.inf)
    np.minimum.at(minnb, ra, width[rb])
    np.minimum.at(minnb, rb, width[ra])
    plane = coplanar & (area >= min_area) & (width >= 3.0 * minnb)
    plane |= coplanar & (area >= min_area) & np.isinf(minnb)
    # new patches: each accepted plane alone; the rest of each smooth patch
    # split into its connected remainders
    is_plane_face = plane[reg]
    keepj = two & ~hard & (is_plane_face[np.maximum(cm.ef[:, 0], 0)] == is_plane_face[np.maximum(cm.ef[:, 1], 0)])
    keepj &= ~(is_plane_face[np.maximum(cm.ef[:, 0], 0)] & (reg[np.maximum(cm.ef[:, 0], 0)] != reg[np.maximum(cm.ef[:, 1], 0)]))
    return _labels(m, cm.ef[keepj])


def planar_patches(cm, patch, plane_deg=0.5):
    """True per patch whose normals all agree within plane_deg."""
    npatch = int(patch.max()) + 1
    out = np.zeros(npatch, bool)
    order = np.argsort(patch, kind="stable")
    bounds = np.searchsorted(patch[order], np.arange(npatch + 1))
    c = math.cos(math.radians(plane_deg))
    for p in range(npatch):
        f = order[bounds[p]:bounds[p + 1]]
        if len(f):
            mean = (cm.N[f] * cm.A[f, None]).sum(0)
            L = np.linalg.norm(mean)
            out[p] = L > 0 and (cm.N[f] @ (mean / L)).min() >= c
    return out


def make_pieces(cm, patch, hard, flat, uv_good, max_dev=MAX_DEV):
    """Piece id per face and a Piece per id: one piece per flat face ('plane'),
    one per curved face that unwrapped cleanly ('uv'); the rest are cut into
    projection pieces ('proj') whose normals stay within max_dev."""
    from .retess import Piece
    m = len(cm.T)
    sub = np.full(m, -1, np.int64)
    pieces = []
    cosmax = math.cos(max_dev)
    npatch = int(patch.max()) + 1
    order = np.argsort(patch, kind="stable")
    bounds = np.searchsorted(patch[order], np.arange(npatch + 1))
    soft = ~hard
    for p in range(npatch):
        faces = order[bounds[p]:bounds[p + 1]]
        if len(faces) == 0:
            continue
        n, a = cm.N[faces], cm.A[faces]
        mean = (n * a[:, None]).sum(0)
        L = np.linalg.norm(mean)
        mean = mean / L if L > 0 else n[0]
        if flat[p]:
            sub[faces] = len(pieces)
            pieces.append(Piece('plane', mean))
            continue
        if uv_good[p]:
            sub[faces] = len(pieces)
            pieces.append(Piece('uv'))
            continue
        if (n @ mean).min() >= cosmax:
            sub[faces] = len(pieces)
            pieces.append(Piece('proj', mean))
            continue
        lab, cents = _cluster(n, a, cosmax)
        local = {int(f): i for i, f in enumerate(faces)}
        pairs = []
        for i, f in enumerate(faces):
            for k in range(3):
                e = cm.TE[f, k]
                if not soft[e]:
                    continue
                g = cm.ef[e, 0] if cm.ef[e, 0] != f else cm.ef[e, 1]
                j = local.get(int(g))
                if j is not None and lab[j] == lab[i]:
                    pairs.append((i, j))
        comp = _labels(len(faces), pairs)
        for c in range(comp.max() + 1):
            sel = faces[comp == c]
            d = (cm.N[sel] * cm.A[sel, None]).sum(0)
            d /= max(np.linalg.norm(d), 1e-30)
            if (cm.N[sel] @ d).min() < cosmax * 0.5:
                d = cents[lab[np.flatnonzero(comp == c)[0]]]
            sub[sel] = len(pieces)
            pieces.append(Piece('proj', d))
    return sub, pieces


def _cluster(n, a, cosmax):
    """Spherical k-means on normals until every face is within max_dev."""
    for k in range(2, 17):
        cents = [n[int(np.argmax(a))]]
        for _ in range(k - 1):                           # farthest-point init
            d = np.max(np.stack([n @ c for c in cents]), axis=0)
            cents.append(n[int(np.argmin(d))])
        C = np.array(cents)
        for _ in range(10):
            lab = np.argmax(n @ C.T, axis=1)
            for j in range(k):
                s = (n[lab == j] * a[lab == j, None]).sum(0)
                L = np.linalg.norm(s)
                if L > 0:
                    C[j] = s / L
        lab = np.argmax(n @ C.T, axis=1)
        if np.all(np.einsum("ij,ij->i", n, C[lab]) >= cosmax):
            return lab, list(C)
    return lab, list(C)


# --------------------------------------------------------------------------
# Curve network
# --------------------------------------------------------------------------

def boundary_edges(cm, sub, hard):
    """Edges that must survive: sub-patch borders, open/non-manifold, creases."""
    two = cm.ecount == 2
    b = ~two | hard
    b[two] |= sub[cm.ef[two, 0]] != sub[cm.ef[two, 1]]
    return b


def simplify_chains(cm, sub, bnd, tol, hard=None, vn=None, min_loop=3, surf=None):
    """Walk boundary edges corner to corner, DP each chain. Returns keep mask.

    A chain made only of soft edges (a border between pieces of one smooth
    surface, or a tangent border between a flat face and a fillet) is
    reduced by its deviation along the surface normal: sliding along the
    surface changes nothing you can see, so zig-zag borders collapse to a
    few points. But a shortcut must still run on the surface: `surf` gives
    the distance of points to the source, and a chord that leaves it (an arc
    cut straight across the hole it runs around) keeps its vertices. Creases,
    open edges and material borders keep the full 3D test, because where a
    crease runs is part of the shape.
    """
    V, E = cm.V, cm.E
    nv = len(V)
    keep = np.zeros(nv, bool)
    be = np.flatnonzero(bnd)
    if len(be) == 0:
        return keep
    adj = {}
    ehard = {}
    for e in be:
        a, b = int(E[e, 0]), int(E[e, 1])
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
        ehard[(min(a, b), max(a, b))] = True if hard is None else bool(hard[e])

    def dp(chain):
        soft = vn is not None and not any(ehard[(min(x, y), max(x, y))]
                                          for x, y in zip(chain, chain[1:]))
        _keep_dp(V, chain, tol, keep, vn if soft else None, surf if soft else None)
    # corners: not degree 2, or three or more sub-patches meet
    vs = np.concatenate([cm.T[:, 0], cm.T[:, 1], cm.T[:, 2]])
    ss = np.concatenate([sub, sub, sub])
    pairs = np.unique(np.stack([vs, ss], 1), axis=0)
    nsub = np.bincount(pairs[:, 0], minlength=nv)
    corner = set(v for v, nb in adj.items() if len(nb) != 2 or nsub[v] >= 3)
    for v in corner:
        keep[v] = True
    used = set()

    def walk(a, b):
        chain = [a, b]
        used.add((min(a, b), max(a, b)))
        prev, cur = a, b
        while cur not in corner:
            nxt = [w for w in adj[cur] if w != prev and (min(cur, w), max(cur, w)) not in used]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            used.add((min(prev, cur), max(prev, cur)))
            chain.append(cur)
        return chain

    for c in corner:
        for nb in adj[c]:
            if (min(c, nb), max(c, nb)) not in used:
                dp(walk(c, nb))
    # closed loops with no corner at all (a hole rim, a circular edge)
    for v in list(adj):
        for nb in adj[v]:
            if (min(v, nb), max(v, nb)) in used:
                continue
            corner.add(v)
            chain = walk(v, nb)
            corner.discard(v)
            if len(chain) > 2 and chain[-1] == chain[0]:
                chain = chain[:-1]
                j = int(np.argmax(np.linalg.norm(V[chain] - V[chain[0]], axis=1)))
                keep[chain[0]] = keep[chain[j]] = True
                dp(chain[:j + 1])
                dp(chain[j:] + [chain[0]])
                ring = [w for w in chain if keep[w]]
                if len(ring) < min_loop:                     # never a sliver loop
                    step = max(1, len(chain) // min_loop)
                    for w in chain[::step][:min_loop]:
                        keep[w] = True
            else:
                dp(chain)
    return keep


def _keep_dp(V, chain, tol, keep, vn=None, surf=None):
    keep[chain[0]] = keep[chain[-1]] = True
    if len(chain) < 3:
        return
    P = V[chain]
    Nn = vn[chain] if vn is not None else None
    st = [(0, len(chain) - 1)]
    while st:
        a, b = st.pop()
        if b - a < 2:
            continue
        d = P[b] - P[a]
        L = np.linalg.norm(d)
        seg = P[a + 1:b] - P[a]
        if L < 1e-30:
            off = seg
        else:
            u = d / L
            off = seg - np.outer(seg @ u, u)          # from the chord to the point
        if Nn is None:
            err = np.linalg.norm(off, axis=1)
        else:
            err = np.abs(np.einsum("ij,ij->i", off, Nn[a + 1:b]))
            if surf is not None and err.max() <= tol:  # the shortcut must stay on the surface
                err = np.maximum(err, surf(P[a + 1:b] - off))
        k = int(err.argmax())
        if err[k] > tol:
            c = a + 1 + k
            keep[chain[c]] = True
            st += [(a, c), (c, b)]
