"""Point-to-surface distance that can be trusted on sheet-metal CAD.

mathutils' BVHTree computes closest points in single precision. On the long
thin triangles sheet-metal parts are made of (a 400 x 1.7 mm strip is two of
them) the answer can be off by half a millimetre - more than the tolerance
being checked - so an exact copy of a part reads as out of tolerance, and
the tree's own range query does not even list the right triangle.

So the tree is only asked for a candidate triangle. The exact distance to
it, in double precision, is a true upper bound; every triangle whose
bounding box lies within that bound is then measured exactly as well.
"""

import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree


def _cross(a, b):
    """Row-wise cross product (np.cross is slow on the small arrays used here)."""
    return np.stack([a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
                     a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
                     a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]], 1)


def _seg(p, A, B):
    d = B - A
    L = np.einsum("ij,ij->i", d, d)
    t = np.clip(np.einsum("ij,ij->i", p - A, d) / np.where(L > 0, L, 1.0), 0.0, 1.0)
    return np.linalg.norm(p - (A + t[:, None] * d), axis=1)


def point_tri(p, A, B, C):
    """Exact distance from p (one point, or one per row) to triangles A, B, C."""
    n = _cross(B - A, C - A)
    nn = np.einsum("ij,ij->i", n, n)
    ok = nn > 0
    # the projection falls inside when p is on the inner side of all three edges
    inside = ok.copy()
    for P0, P1 in ((A, B), (B, C), (C, A)):
        inside &= np.einsum("ij,ij->i", _cross(P1 - P0, p - P0), n) >= 0
    plane = np.abs(np.einsum("ij,ij->i", p - A, n)) / np.sqrt(np.where(ok, nn, 1.0))
    edge = np.minimum(np.minimum(_seg(p, A, B), _seg(p, B, C)), _seg(p, C, A))
    return np.where(inside, plane, edge)


class Surface:
    """Distance of points to a triangle set, exact where it matters.

    A distance at or below `exact_above` is the exact distance to the
    tree's candidate: a true upper bound, and small. Anything larger is the
    exact minimum over every triangle that could be nearer.
    """

    def __init__(self, V, T, exact_above=0.0):
        self.V = np.asarray(V, float).reshape(-1, 3)
        self.T = np.asarray(T, np.int64).reshape(-1, 3)
        self.tree = BVHTree.FromPolygons(self.V.tolist(), self.T.tolist())
        X = self.V[self.T]
        self.lo, self.hi = X.min(1), X.max(1)
        self.exact_above = exact_above
        self.cells = None

    def __call__(self, P, exact_above=None):
        P = np.asarray(P, float).reshape(-1, 3)
        lim = self.exact_above if exact_above is None else exact_above
        if not len(self.T):
            return np.full(len(P), np.inf)
        idx = np.fromiter((self.tree.find_nearest(Vector(p))[2] for p in P), np.int64, len(P))
        X = self.V[self.T[idx]]
        d = point_tri(P, X[:, 0], X[:, 1], X[:, 2])
        for i in np.flatnonzero(d > lim):
            c = self._near_boxes(P[i], d[i])
            if len(c):
                Y = self.V[self.T[c]]
                d[i] = min(d[i], float(point_tri(P[i], Y[:, 0], Y[:, 1], Y[:, 2]).min()))
        return d

    def exact_max(self, P, d):
        """The exact largest distance of points P, given upper bounds d (as
        returned by a call): refine from the largest bound down, and stop
        once no remaining bound can beat what has been found."""
        best = 0.0
        for i in np.argsort(-d):
            if d[i] <= best:
                break
            c = self._near_boxes(P[i], d[i])
            Y = self.V[self.T[c]]
            e = float(point_tri(P[i], Y[:, 0], Y[:, 1], Y[:, 2]).min()) if len(c) else d[i]
            best = max(best, min(d[i], e))
        return best

    def _near_boxes(self, p, r):
        """Triangles whose bounding box comes within r of p."""
        if self.cells is None:
            self._index()
        a = np.floor((p - r) / self.h).astype(np.int64)
        b = np.floor((p + r) / self.h).astype(np.int64)
        if np.prod(b - a + 1) > 512:                     # a wide query: scan everything
            c = np.arange(len(self.T))
        else:
            got = [self.cells.get((i, j, k)) for i in range(a[0], b[0] + 1)
                   for j in range(a[1], b[1] + 1) for k in range(a[2], b[2] + 1)]
            got = [t for t in got if t is not None]
            # large boxes live only in `big`, so the two never overlap
            c = np.concatenate([np.unique(np.concatenate(got)) if got else self.big[:0], self.big])
        keep = np.all((self.lo[c] <= p + r) & (self.hi[c] >= p - r), axis=1)
        return c[keep]

    def _index(self):
        """Uniform grid over the bounding boxes; boxes spanning many cells
        (slivers, large plates) go to a list that is always checked."""
        ext = (self.hi - self.lo).max(1)
        self.h = max(float(np.median(ext)) * 2.0, 1e-9)
        a = np.floor(self.lo / self.h).astype(np.int64)
        b = np.floor(self.hi / self.h).astype(np.int64)
        big = np.prod(b - a + 1, axis=1) > 27
        self.big = np.flatnonzero(big)
        cells = {}
        for t in np.flatnonzero(~big):
            for i in range(a[t, 0], b[t, 0] + 1):
                for j in range(a[t, 1], b[t, 1] + 1):
                    for k in range(a[t, 2], b[t, 2] + 1):
                        cells.setdefault((i, j, k), []).append(t)
        self.cells = {key: np.array(v, np.int64) for key, v in cells.items()}
