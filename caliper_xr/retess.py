"""Re-triangulate every recovered face lean, on the original surface.

Each piece has a one-to-one 2D domain: its plane (flat faces), its whole-face
UV unwrap (curved faces), or a projection (pieces of faces that would not
unwrap). It is rebuilt as a constrained Delaunay triangulation of its
simplified outline; points are then added only where the result is farther
than the tolerance from the source surface, measured in 3D along the surface
normal: a source vertex where one is there, otherwise a point placed exactly
on the source triangle at the worst spot. Flat faces get no interior points.
Every vertex lies on the source surface, and a piece is never allowed to come
out heavier than its source: then its source triangles are kept, its outline
is kept whole, and the neighbours are rebuilt against it, so no crack opens.
"""

import numpy as np

from .tri2d import (Fail, area2, bary, cdt, frame, outline_loops, probe,
                    reduce_constraints, tree2d)


class Result:
    def __init__(self):
        self.tris = []         # (a, b, c): < nv source vertex, >= nv added point
        self.tsub = []         # piece per triangle
        self.pts = []          # added points (on the source surface)
        self.pnrm = []         # their normals
        self.puv = []          # their UVs (or None)
        self.fallback = 0


class Piece:
    """How a piece is flattened: 'plane'/'proj' use direction d, 'uv' uses UV."""
    __slots__ = ("mode", "d")

    def __init__(self, mode, d=None):
        self.mode, self.d = mode, d


# --------------------------------------------------------------------------

def rebuild(cm, sub, pieces, UV, bnd, keep, tol, patch, vnorm, progress=None,
            state=None, only=None, near=None):
    """Rebuild pieces (all, or only those in `only`, reusing `state`).
    `near(points)` is the distance of 3D points to the source surface: no
    piece may leave it by more than tol (a fill over a hole, a shortcut
    across a notch). Returns (Result, state) so a later pass can redo a few."""
    nsub = len(pieces)
    order = np.argsort(sub, kind="stable")
    bounds = np.searchsorted(sub[order], np.arange(nsub + 1))
    if state is None:
        state = {"done": {}, "fell": set()}
    done, fell = state["done"], state["fell"]
    todo = set(range(nsub)) if only is None else set(only)
    for _pass in range(40):
        grown = set()
        for n, s in enumerate(sorted(todo)):
            if progress and n % 100 == 0:
                progress(n, len(todo))
            faces = order[bounds[s]:bounds[s + 1]]
            if len(faces) == 0 or s in fell:
                done[s] = ([tuple(int(v) for v in cm.T[f]) for f in faces], [])
                continue
            try:
                done[s] = _piece(cm, faces, pieces[s], UV, bnd, keep, tol, vnorm,
                                 int(patch[faces[0]]), near)
            except Fail as ex:
                tries = state.setdefault("tries", {})
                tries[s] = tries.get(s, 0) + 1
                fix = np.array([v for v in ex.restore if not keep[v]], np.int64)
                if len(fix) and tries[s] <= 12:
                    # restore the few outline vertices that caused it, and
                    # rebuild this piece and every piece sharing them
                    keep[fix] = True
                    touched = np.flatnonzero(np.isin(cm.T, fix).any(axis=1))
                    grown.update(int(x) for x in np.unique(sub[touched]))
                    grown.add(s)
                    continue
                fell.add(s)
                state.setdefault("why", {})[s] = ex.why
                done[s] = ([tuple(int(v) for v in cm.T[f]) for f in faces], [])
                grown |= _hold_outline(cm, sub, faces, bnd, keep)
        todo = grown - fell
        if not todo:
            break
    res = Result()
    nv = len(cm.V)
    for s in range(nsub):
        tris, extra = done.get(s, ([], []))
        base = nv + len(res.pts)
        for t in tris:
            res.tris.append(tuple(g if g >= 0 else base + (-g - 1) for g in t))
            res.tsub.append(s)
        for pos, nrm, uv in extra:
            res.pts.append(pos)
            res.pnrm.append(nrm)
            res.puv.append(uv)
    res.fallback = len(fell)
    return res, state


def force_source(cm, sub, bnd, keep, state, forced):
    """Put pieces back on their exact source triangles. Their outline is kept
    whole; returns the neighbouring pieces that must be rebuilt against it."""
    nsub = int(sub.max()) + 1
    order = np.argsort(sub, kind="stable")
    bounds = np.searchsorted(sub[order], np.arange(nsub + 1))
    fell = state["fell"]
    redo = set()
    for s in forced:
        fell.add(s)
        state.setdefault("why", {})[s] = "tolerance"
        redo |= _hold_outline(cm, sub, order[bounds[s]:bounds[s + 1]], bnd, keep)
    return (redo - fell) | set(forced)


def _hold_outline(cm, sub, faces, bnd, keep):
    """Keep a piece's whole outline; returns the pieces sharing what was added."""
    vs = np.unique(cm.T[faces])
    on_b = vs[_on_boundary(cm, faces, bnd, vs)]
    newly = on_b[~keep[on_b]]
    if not len(newly):
        return set()
    keep[newly] = True
    touched = np.flatnonzero(np.isin(cm.T, newly).any(axis=1))
    return set(int(x) for x in np.unique(sub[touched]))


def _on_boundary(cm, faces, bnd, vs):
    es = np.unique(cm.TE[faces].ravel())
    bv = np.unique(cm.E[es[bnd[es]]].ravel())
    return np.isin(vs, bv)


def _piece(cm, faces, piece, UV, bnd, keep, tol, vnorm, patch_id, near=None):
    fset = set(int(f) for f in faces)
    loops, cons = outline_loops(cm, faces, fset, bnd)
    if not loops:
        raise Fail("no outline")
    V = cm.V
    tri3 = V[cm.T[faces]]                                   # (k,3,3)
    if piece.mode == 'uv':
        uvmap = {}
        for f, uvs in zip(faces, UV[faces]):
            for k in range(3):
                uvmap[int(cm.T[f, k])] = uvs[k]
        src2 = UV[faces]

        def coord(ids):
            return np.array([uvmap[int(i)] for i in ids]).reshape(-1, 2)
    else:
        e1, e2 = frame(piece.d)
        o = V[cm.T[faces[0], 0]]
        B = np.stack([e1, e2], 1)
        src2 = np.stack([(tri3[:, j] - o) @ B for j in range(3)], 1)

        def coord(ids):
            return (V[np.asarray(ids, np.int64)] - o) @ B

    # a mirrored island turns counter-clockwise 2D triangles inside out
    sa = ((src2[:, 1, 0] - src2[:, 0, 0]) * (src2[:, 2, 1] - src2[:, 0, 1])
          - (src2[:, 1, 1] - src2[:, 0, 1]) * (src2[:, 2, 0] - src2[:, 0, 0]))
    mirrored = (sa < 0).sum() > (sa > 0).sum()

    kloops, gaps = [], []          # gaps[i][j]: dropped vertices after kloops[i][j]
    for lp in loops:
        k = [v for v in lp if keep[v]]
        if len(k) >= 3:
            a = area2(coord(k))
            if abs(a) > 0:
                full = lp if a > 0 else lp[::-1]
                k = k if a > 0 else k[::-1]
                start = full.index(k[0])
                full = full[start:] + full[:start]
                g, cur = [], []
                for v in full[1:] + [full[0]]:
                    if keep[v]:
                        g.append(cur)
                        cur = []
                    else:
                        cur.append(v)
                kloops.append(k)
                gaps.append(g)
    if not kloops:
        return [], []
    segs = reduce_constraints(cons, keep)

    # every outline segment is also an explicit, numbered CDT edge, so a
    # crossing reported by the CDT maps straight back to (loop, segment)
    seg_of = [None] * len(segs)
    for fi, k in enumerate(kloops):
        for j in range(len(k)):
            seg_of.append((fi, j))

    def farthest(fi, j):
        """The dropped vertex farthest from the shortcut kloops[fi][j] -> next."""
        drop = gaps[fi][j] if j < len(gaps[fi]) else []
        if not drop:
            return None
        a, b = kloops[fi][j], kloops[fi][(j + 1) % len(kloops[fi])]
        pa, pb = coord([a])[0], coord([b])[0]
        q, d = coord(drop), pb - pa
        off = np.abs((q[:, 0] - pa[0]) * d[1] - (q[:, 1] - pa[1]) * d[0])
        return drop[int(off.argmax())]

    def restore_for_edges(edge_ids):
        """Dropped vertices that, put back, uncross the given CDT input edges."""
        out = [farthest(*seg_of[e]) for e in edge_ids if e < len(seg_of) and seg_of[e] is not None]
        return [v for v in out if v is not None]

    def restore_near(points2):
        """The dropped outline vertex nearest each 2D point."""
        drop = [v for g in gaps for seg in g for v in seg]
        if not drop:
            return []
        Q = coord(drop)
        return list({drop[int(np.argmin(np.hypot(*(Q - p).T)))] for p in points2})

    if near is not None:
        cut = [(fi, j) for fi, g in enumerate(gaps) for j in range(len(g)) if g[j]]
        if cut:
            mid = np.array([V[kloops[fi][j]] + V[kloops[fi][(j + 1) % len(kloops[fi])]]
                            for fi, j in cut]) * 0.5
            off = near(mid) > tol
            if off.any():
                raise Fail("outline leaves the surface",
                           [farthest(fi, j) for (fi, j), o in zip(cut, off) if o])

    pid, ids = {}, []
    pinned = set()                 # source edges forced back where wiring matters
    newp = []                      # points added on the source surface: (pos, nrm, uv, xy)

    def add(g):
        g = int(g)
        if g not in pid:
            pid[g] = len(ids)
            ids.append(g)
        return pid[g]

    def add_new(f, w):
        """A point on source triangle f at barycentric w: exactly on the surface."""
        w = np.clip(w, 0.0, None)
        w = w / w.sum()
        fg = faces[f]
        nrm = w @ cm.CN[fg]
        nrm = nrm / max(np.linalg.norm(nrm), 1e-30)
        uv = (w @ cm.UV[fg]) if getattr(cm, "UV", None) is not None else None
        newp.append((w @ tri3[f], nrm, uv, w @ src2[f]))
        g = -len(newp)
        pid[g] = len(ids)
        ids.append(g)

    def coords(ids_):
        out = np.array([newp[-g - 1][3] if g < 0 else (0.0, 0.0) for g in ids_]).reshape(-1, 2)
        own = [i for i, g in enumerate(ids_) if g >= 0]
        if own:
            out[own] = coord([ids_[i] for i in own])
        return out

    def positions(ids_):
        return np.array([V[g] if g >= 0 else newp[-g - 1][0] for g in ids_])

    for lp in kloops:
        for g in lp:
            add(g)
    for a, b in segs:
        add(a)
        add(b)

    def triangulate():
        P2 = coords(ids)
        fl = [[pid[g] for g in lp] for lp in kloops]
        el = [(pid[a], pid[b]) for a, b in segs] + [(f[j], f[(j + 1) % len(f)]) for f in fl for j in range(len(f))]
        el += [(pid[a], pid[b]) for a, b in pinned if a in pid and b in pid]
        try:
            tris = cdt(P2, fl, el)
        except Fail as ex:
            if ex.restore:
                raise Fail(ex.why, restore_for_edges(ex.restore))
            raise
        return P2, tris

    def emit(tris):
        out = [(ids[a], ids[b], ids[c]) for a, b, c in tris]
        out = [(a, c, b) for a, b, c in out] if mirrored else out
        return out, [(q[0], q[1], q[2]) for q in newp]

    if piece.mode == 'plane':
        P2, tris = triangulate()
        if len(tris) > len(faces):
            raise Fail("plane heavier")                   # never heavier than the source
        if near is not None and tris:
            loc = np.array(tris)
            off = near(positions(ids)[loc].mean(1)) > tol   # a fill over a hole or a notch
            if off.any():
                raise Fail("fill leaves the surface", restore_near(P2[loc[off]].mean(1)))
        return emit(tris)

    allv = np.unique(cm.T[faces])
    inner = allv[~_on_boundary(cm, faces, bnd, allv)]
    ring = set(v for lp in loops for v in lp)
    if len(inner) == 0 and all(keep[v] for v in ring):
        return [tuple(int(v) for v in cm.T[f]) for f in faces], []   # already the leanest

    nin = len(inner)
    innerset = set(int(v) for v in inner)
    src_edges = {}
    for f in faces:
        for k in range(3):
            a, b = int(cm.T[f, k]), int(cm.T[f, (k + 1) % 3])
            if a in innerset or b in innerset:
                src_edges[(min(a, b), max(a, b))] = f
    s2 = np.vstack([coord(inner) if nin else np.zeros((0, 2)), src2.mean(1)])
    s3 = np.vstack([V[inner] if nin else np.zeros((0, 3)), tri3.mean(1)])
    sn = np.vstack([np.array([vnorm.get((int(v), patch_id), (0, 0, 1.0)) for v in inner]).reshape(-1, 3),
                    cm.N[faces]])
    stree = tree2d(src2.reshape(-1, 2), np.arange(len(faces) * 3).reshape(-1, 3))
    third = np.array([1 / 3, 1 / 3, 1 / 3])
    placed = set()                                         # (face, bary) already added

    for _round in range(120):
        P2, tris = triangulate()
        if not tris:
            raise Fail("empty")
        if len(tris) >= len(faces):
            raise Fail("heavier")                         # never heavier than the source
        loc = np.array(tris)
        X3 = positions(ids)
        ntree = tree2d(P2, loc)
        bad, best = {}, {}                                 # tri -> (err, xy, face, bary)
        for si in range(len(s2)):
            x, y = s2[si]
            t = probe(ntree, x, y)
            if t is None:
                continue
            a, b, c = loc[t]
            w = bary(s2[si], P2[a], P2[b], P2[c])
            q = w[0] * X3[a] + w[1] * X3[b] + w[2] * X3[c]
            err = abs((s3[si] - q) @ sn[si])
            if err > tol and err > bad.get(t, (0.0,))[0]:
                bad[t] = (err, s2[si], si - nin if si >= nin else None, third)
            if si < nin and int(inner[si]) not in pid and err > best.get(t, (-1.0, 0))[0]:
                best[t] = (err, si)
        cen2 = P2[loc].mean(1)
        cen3 = X3[loc].mean(1)
        for t in range(len(loc)):
            ft = probe(stree, cen2[t][0], cen2[t][1])
            if ft is not None:
                w = bary(cen2[t], src2[ft, 0], src2[ft, 1], src2[ft, 2])
                spot = (abs((w @ tri3[ft] - cen3[t]) @ cm.N[faces[ft]]), cen2[t], ft, w)
            elif near is not None:                         # beyond the source region
                spot = (float(near(cen3[t:t + 1])[0]), cen2[t], None, None)
            else:
                continue
            if spot[0] > tol and spot[0] > bad.get(t, (0.0,))[0]:
                bad[t] = spot
        if not bad:
            return emit(tris)
        # 1. source vertices inside the bad triangles
        fresh = sorted({int(inner[best[t][1]]) for t in bad if t in best})
        for g in fresh:
            add(g)
        if fresh:
            continue
        # 2. wiring: pin the source's own edges there, once
        before = len(pinned)
        badset = set(bad)
        for (a, b), f in src_edges.items():
            if (a, b) in pinned or a not in pid or b not in pid:
                continue
            pa, pb = P2[pid[a]], P2[pid[b]]
            mid = 0.5 * (pa + pb)
            if probe(ntree, mid[0], mid[1]) in badset:
                pinned.add((a, b))
        if len(pinned) > before:
            continue
        # 3. no source vertex left there: add the worst spot itself, a point
        #    exactly on the source triangle it lies in
        before = len(placed)
        for _err, _xy, f, w in bad.values():
            if f is not None and (int(f), tuple(np.round(w, 4))) not in placed:
                placed.add((int(f), tuple(np.round(w, 4))))
                add_new(int(f), np.asarray(w, float))
        if len(placed) > before:
            continue
        # 4. what is left sits against the outline: put back the outline
        #    vertices nearest the worst spots
        spots = [v[1] for v in sorted(bad.values(), key=lambda x: -x[0])[:16]]
        raise Fail("residual", restore_near(spots))
    raise Fail("rounds")
