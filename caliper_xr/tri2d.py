"""2D building blocks of the rebuild: outline loops of a piece, constrained
Delaunay triangulation, and point location in a flattened triangulation."""

import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree
from mathutils.geometry import delaunay_2d_cdt


class Fail(Exception):
    """Carries why a piece could not be rebuilt (for the report)."""

    def __init__(self, why="failed", restore=()):
        super().__init__(why)
        self.why = why
        self.restore = list(restore)


def frame(d):
    a = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = a - d * (a @ d)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(d, e1)


def outline_loops(cm, faces, fset, bnd):
    """Directed boundary loops (region on the left) and inner crease edges."""
    T, TE, ef = cm.T, cm.TE, cm.ef
    out, cons = {}, set()
    for f in faces:
        for k in range(3):
            e = TE[f, k]
            if not bnd[e]:
                continue
            a, b = int(T[f, k]), int(T[f, (k + 1) % 3])
            g0, g1 = ef[e]
            if g0 >= 0 and int(g0) in fset and int(g1) in fset:
                cons.add((min(a, b), max(a, b)))
            else:
                out.setdefault(a, []).append(b)
    loops = []
    while out:
        a = next(iter(out))
        loop, cur, closed = [a], a, False
        for _ in range(len(T) * 3):
            nxt = out.get(cur)
            if not nxt:
                break
            b = nxt.pop()
            if not nxt:
                del out[cur]
            if b == loop[0]:
                closed = True
                break
            loop.append(b)
            cur = b
        if not closed:
            raise Fail("open outline")
        loops.append(loop)
    return loops, cons


def reduce_constraints(cons, keep):
    adj = {}
    for a, b in cons:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    segs, seen = [], set()
    for s in [v for v in adj if keep[v]]:
        for nb in adj[s]:
            if (s, nb) in seen:
                continue
            prev, cur = s, nb
            seen.add((s, nb))
            while not keep[cur] and len(adj[cur]) == 2:
                nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
                prev, cur = cur, nxt
            seen.add((cur, prev))
            if keep[cur] and cur != s:
                segs.append((s, cur))
    return segs


def area2(P):
    return 0.5 * float(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))


def bary(p, a, b, c):
    v0, v1, v2 = b - a, c - a, p - a
    d00, d01, d11 = v0 @ v0, v0 @ v1, v1 @ v1
    d20, d21 = v2 @ v0, v2 @ v1
    den = d00 * d11 - d01 * d01
    if abs(den) < 1e-300:
        return np.array([1.0, 0, 0])
    v = (d11 * d20 - d01 * d21) / den
    w = (d00 * d21 - d01 * d20) / den
    return np.array([1 - v - w, v, w])


def tree2d(P2, tris):
    return BVHTree.FromPolygons([Vector((x, y, 0.0)) for x, y in P2], [list(t) for t in tris])


def probe(tree, x, y):
    return tree.ray_cast(Vector((x, y, 1.0)), Vector((0, 0, -1.0)), 2.0)[2]


def cdt(P2, loops, segs):
    coords = [Vector((float(x), float(y))) for x, y in P2]
    ext = float(np.ptp(np.asarray(P2), axis=0).max()) if len(P2) > 1 else 1.0
    try:
        vo, eo, fo, ov, oe, of = delaunay_2d_cdt(coords, segs, loops, 0, max(ext * 1e-11, 1e-18), True)
    except Exception:
        raise Fail("cdt error")
    if len(vo) != len(coords):
        # which input edges meet at the extra vertices: those segments cross
        extra = set(i for i in range(len(vo)) if not ov[i])
        crossing = set()
        for (a, b), orig in zip(eo, oe):
            if a in extra or b in extra:
                crossing.update(orig)
        raise Fail("outline crosses itself", crossing)
    remap = [ov[i][0] if ov[i] else -1 for i in range(len(vo))]
    out = []
    for f, orig in zip(fo, of):
        if len(orig) % 2 == 0 or len(f) != 3:
            continue
        a, b, c = (remap[i] for i in f)
        if min(a, b, c) < 0:
            raise Fail("new vertex")
        pa, pb, pc = P2[a], P2[b], P2[c]
        cr = (pb[0] - pa[0]) * (pc[1] - pa[1]) - (pb[1] - pa[1]) * (pc[0] - pa[0])
        if cr != 0:
            out.append((a, b, c) if cr > 0 else (a, c, b))
    return out
