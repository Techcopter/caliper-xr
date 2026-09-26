"""Outer shell: which faces border the air outside, and which are internal.

The geometry is sampled onto a voxel grid, openings narrower than the gap
are sealed, and the air reachable from the grid border is flood-filled. A
face is on the shell when the first air it meets along its normal, either
way, is that outside air - reached before any other surface. Internal faces
never meet it: inner bodies, the walls of closed cavities, faces between
touching bodies, whatever sits behind a seam or vent narrower than the gap.

Reachable is all-or-nothing for a cavity: one opening wider than the gap
makes the whole interior outside air. What is merely out of sight in there
is the camera test's job (visibility.py); the rebuild combines the two.

The grid cannot leak. Surfaces are sampled at half a voxel, so a chain of
face-adjacent air voxels from inside to outside would have to cross a
surface within a quarter voxel of a sample - inside one of its voxels.
"""

import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree

from .visibility import world_tris

MAX_VOXELS = 24_000_000
MAX_AXIS = 640
CHUNK = 50_000


def _lattice(n):
    """Barycentric weights of the vertices of an n x n subdivision."""
    i, j = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing="ij")
    keep = i + j <= n
    i, j = i[keep], j[keep]
    return (np.stack([n - i - j, i, j], 1) / n).astype(np.float32)


LATTICE = {n: _lattice(n) for n in range(1, 9)}


def _edges2(X):
    """Squared edge lengths (v0v1, v1v2, v2v0) of triangles X (m,3,3)."""
    e = np.stack([X[:, 1] - X[:, 0], X[:, 2] - X[:, 1], X[:, 0] - X[:, 2]], 1)
    return np.einsum("mij,mij->mi", e, e)


def _split(X, dmax, own):
    """Longest-edge bisection until no edge is longer than dmax; `own`
    follows every piece back to the triangle it came from."""
    out, oo = [], []
    while len(X):
        L = _edges2(X)
        k = L.argmax(1)
        big = L[np.arange(len(X)), k] > dmax * dmax
        out.append(X[~big])
        oo.append(own[~big])
        if not big.any():
            break
        X, k, own = X[big], k[big], own[big]
        X = X[np.arange(len(X))[:, None], (np.arange(3)[None] + k[:, None]) % 3]
        M = 0.5 * (X[:, 0] + X[:, 1])
        X = np.concatenate([np.stack([X[:, 0], M, X[:, 2]], 1), np.stack([M, X[:, 1], X[:, 2]], 1)])
        own = np.concatenate([own, own])
    return np.concatenate(out), np.concatenate(oo)


def _occupancy(tris, origin, h, shape):
    """Voxels holding surface. Coordinates go to voxel units first, so every
    triangle is sampled on a lattice no coarser than half a voxel."""
    occ = np.zeros(shape, bool)
    for X in tris:
        for a in range(0, len(X), CHUNK):
            Y = ((X[a:a + CHUNK] - origin) / h).astype(np.float32)
            Y, _ = _split(Y, 4.0, np.zeros(len(Y), np.int64))
            n = np.clip(np.ceil(np.sqrt(_edges2(Y).max(1)) / 0.5), 1, 8).astype(np.int64)
            for k in np.unique(n):
                Z = Y[n == k]
                for b in range(0, len(Z), CHUNK):
                    P = np.einsum("wj,mjc->mwc", LATTICE[int(k)], Z[b:b + CHUNK]).reshape(-1, 3)
                    ijk = np.clip(np.floor(P).astype(np.int64), 0, np.array(shape) - 1)
                    occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    return occ


def _grow(m, r):
    """Grow a voxel mask by r in every direction (a cube, one axis at a time)."""
    for ax in range(3):
        a = np.moveaxis(m, ax, 0)
        for _ in range(r):
            b = a.copy()
            b[1:] |= a[:-1]
            b[:-1] |= a[1:]
            a = b
        m = np.moveaxis(a, 0, ax)
    return np.ascontiguousarray(m)


def _outside(occ):
    """Air connected to the grid border through face-adjacent air voxels.
    The border layer is all air and all seeded, so a flat-index step that
    wraps off one row only ever lands on a border voxel already seen."""
    shape = occ.shape
    air = ~occ.ravel()
    seen = np.zeros(air.size, bool)
    b = np.zeros(shape, bool)
    b[[0, -1]] = True
    b[:, [0, -1]] = True
    b[:, :, [0, -1]] = True
    front = np.flatnonzero(b.ravel() & air)
    seen[front] = True
    steps = (1, shape[2], shape[1] * shape[2])
    while len(front):
        nb = np.concatenate([front + d for d in steps] + [front - d for d in steps])
        nb = nb[(nb >= 0) & (nb < air.size)]
        nb = np.unique(nb[air[nb] & ~seen[nb]])
        seen[nb] = True
        front = nb
    return seen.reshape(shape)


def _fan(N):
    """The normal, and six directions 60 degrees off it: a face in a concave
    corner has its normal running along the neighbouring wall, and only a
    tilted ray leaves that wall's voxels for the air beside it."""
    a = np.where(np.abs(N[:, :1]) < 0.9, [[1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0]])
    t1 = np.cross(N, a)
    t1 /= np.linalg.norm(t1, axis=1)[:, None]
    t2 = np.cross(N, t1)
    out = [N]
    for k in range(6):
        ang = k * np.pi / 3.0
        out.append(0.5 * N + 0.866 * (np.cos(ang) * t1 + np.sin(ang) * t2))
    return out


def _rays(P, N, own, n, g, fan=True):
    """For each of n triangles: does some sample of it (P, normals N, owner
    `own`) meet outside air first, along some direction of the fan (or the
    normal only), either side, before any surface?"""
    tree, occ, out, near, origin, h, reach = g
    steps = (np.arange(int(np.ceil(reach / (0.5 * h)))) + 0.5) * (0.5 * h)
    shape = np.array(occ.shape)
    # no outside air within reach: no direction can find any
    v = np.clip(np.floor((P - origin) / h).astype(np.int64), 0, shape - 1)
    live = near[v[:, 0], v[:, 1], v[:, 2]]
    won = np.zeros(n, bool)
    eps = 1e-3 * h
    for F in (_fan(N) if fan else [N]):
        for sign in (1.0, -1.0):
            for a in range(0, len(P), CHUNK):
                idx = a + np.flatnonzero(live[a:a + CHUNK] & ~won[own[a:a + CHUNK]])
                if not len(idx):
                    continue
                D = F[idx] * sign
                ijk = np.floor((P[idx, None] + D[:, None] * steps[None, :, None] - origin) / h).astype(np.int64)
                inb = np.all((ijk >= 0) & (ijk < shape), axis=2)
                ijk = np.clip(ijk, 0, shape - 1)
                i, j, k = ijk[..., 0], ijk[..., 1], ijk[..., 2]
                air = ~(occ[i, j, k] & inb)
                first = air.argmax(1)
                r = np.arange(len(idx))
                ok = air[r, first] & (out[i, j, k] | ~inb)[r, first]
                # the voxels say this way leads outside; a sheet thinner than
                # a voxel (the far side of a wall) may still stand in between.
                # One clear sample settles its triangle.
                for q in np.flatnonzero(ok):
                    o = own[idx[q]]
                    if won[o]:
                        continue
                    p, d, t = P[idx[q]], D[q], steps[first[q]]
                    hit = tree.ray_cast(Vector(p + d * eps), Vector(d), t)
                    won[o] = hit[0] is None or hit[3] + eps > t
    return won


def _classify(X, g):
    """Per triangle of X (m,3,3): does it border outside air?"""
    h = g[5]
    flags = np.zeros(len(X), bool)
    if not len(X):
        return flags
    n = np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0])
    ln = np.linalg.norm(n, axis=1)
    ok = ln > 0
    flags[~ok] = True                        # no normal to decide by: keep it
    n[ok] /= ln[ok, None]
    idx = np.flatnonzero(ok)
    flags |= _rays(X[idx].mean(1), n[idx], idx, len(X), g)
    # a large face its centre did not settle is sampled all over; there the
    # normal is enough - a corner only hides the samples along its edge
    big = idx[~flags[idx] & (_edges2(X[idx]).max(1) > 4.0 * h * h)]
    if len(big):
        Y, own = _split(X[big], 2.0 * h, big)
        flags |= _rays(Y.mean(1), n[own], own, len(X), g, fan=False)
    return flags


def _shell(tris, gap):
    """Shell flags for each triangle array in `tris`, analysed as one."""
    pts = [x.reshape(-1, 3) for x in tris if len(x)]
    if not pts:
        return [np.ones(len(x), bool) for x in tris], 0.0
    P = np.concatenate(pts)
    lo, hi = P.min(0), P.max(0)
    ext = np.maximum(hi - lo, 1e-9)
    h = max(gap / 2.0, float(ext.max()) / MAX_AXIS,
            (float(np.prod(ext)) / MAX_VOXELS) ** (1.0 / 3.0), float(ext.max()) * 1e-6)
    pad = 2
    origin = lo - pad * h
    shape = tuple(int(v) for v in np.ceil(ext / h).astype(np.int64) + 2 * pad + 1)
    occ = _occupancy(tris, origin, h, shape)
    out = _outside(occ)
    tree = BVHTree.FromPolygons(P.tolist(), np.arange(len(P)).reshape(-1, 3).tolist())
    reach = 3.6 * h
    g = (tree, occ, out, _grow(out, int(np.ceil(reach / h)) + 1), origin, h, reach)
    return [_classify(x, g) for x in tris], h


def analyse(objs, depsgraph, gap, each=False, progress=None):
    """Per object name, a bool per loop triangle: on the outer shell.
    `gap` is in world units; `each` gives every object its own shell.
    Returns (flags, voxel sizes used)."""
    tris = world_tris(objs, depsgraph)
    groups = [[i] for i in range(len(objs))] if each else [list(range(len(objs)))]
    flags, hs = {}, []
    for n, grp in enumerate(groups):
        if progress:
            progress(n, len(groups))
        got, h = _shell([tris[i] for i in grp], gap)
        for i, f in zip(grp, got):
            flags[objs[i].name] = f
        if h:
            hs.append(h)
    return flags, hs
