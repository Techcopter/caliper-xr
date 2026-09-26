"""Pure numpy geometry: hole filling, Douglas-Peucker, window filters.

Nothing in here touches bpy, so it can be unit-tested outside Blender.
All lengths are millimetres unless a name says otherwise.
"""

import numpy as np


# --------------------------------------------------------------------------
# Morphological filters (van Herk / Gil-Werman, O(n) per line)
# --------------------------------------------------------------------------

def _running(x, w, fn, circular):
    """Max or min over a window of 2w+1 along axis 0."""
    n = x.shape[0]
    if w <= 0 or n == 0:
        return x.copy()
    if circular:
        if w >= n:
            return np.broadcast_to(fn.reduce(x, axis=0), x.shape).copy()
        pad = np.concatenate([x[-w:], x, x[:w]], axis=0)
    else:
        pad = np.concatenate([np.repeat(x[:1], w, 0), x, np.repeat(x[-1:], w, 0)], axis=0)
    k = 2 * w + 1
    m = pad.shape[0]
    big = -(-m // k) * k
    ident = -np.inf if fn is np.maximum else np.inf
    full = np.full((big,) + pad.shape[1:], ident)
    full[:m] = pad
    blocks = full.reshape((big // k, k) + pad.shape[1:])
    g = fn.accumulate(blocks, axis=1).reshape(full.shape)
    h = fn.accumulate(blocks[:, ::-1], axis=1)[:, ::-1].reshape(full.shape)
    s = np.arange(n)
    return fn(h[s], g[s + k - 1])


def closing(x, w, circular):
    """Grey closing along axis 0: fills dips narrower than the window."""
    return _running(_running(x, w, np.maximum, circular), w, np.minimum, circular)


# --------------------------------------------------------------------------
# Hole filling on the (station x angle) radius grid
# --------------------------------------------------------------------------

def fill_holes(R, dtheta, dz_mm, max_hole_mm, max_hole_deg, mode):
    """Cap holes flush with the surrounding outer surface.

    R is (stations, angles) of radii in mm, NaN where the ray missed or went
    through to the far side. A ray that enters a bolt hole or port hits the
    bore wall, so a hole is a local *dip* in r. A closing fills a dip that is
    bounded in the filtered direction and leaves steps untouched.

    mode PITS:       fill only where BOTH directions bound the dip, i.e. holes
                     and pockets. Circumferential grooves and long axial slots
                     survive, because they are silhouette.
    mode AGGRESSIVE: fill where EITHER direction bounds it (slots, grooves too).

    The angular window is capped by max_hole_deg so the flats of a hex (60 deg
    wide) are never mistaken for a hole and rounded off.
    """
    S, A = R.shape
    valid = np.isfinite(R)
    X = np.where(valid, R, 0.0)

    cA = np.empty_like(X)
    cap = max(1, int(round(np.radians(max_hole_deg) / dtheta / 2)))
    for i in range(S):
        row_ok = valid[i]
        r = float(np.median(X[i][row_ok])) if row_ok.any() else 0.0
        w = cap if r <= 0 else min(cap, int(np.ceil(max_hole_mm / (r * dtheta) / 2)))
        cA[i] = closing(X[i], w, circular=True)

    wz = max(1, int(np.ceil(max_hole_mm / max(dz_mm, 1e-9) / 2)))
    cZ = closing(X, wz, circular=False)

    out = np.minimum(cA, cZ) if mode == 'PITS' else np.maximum(cA, cZ)
    out = np.maximum(out, X)

    # Anything still at zero was a miss nothing could bound: interpolate
    # round the circumference, then along the axis.
    out[out <= 0] = np.nan
    out = _interp_nans(out, axis=1, circular=True)
    out = _interp_nans(out, axis=0, circular=False)

    # Only dips deeper than the hole threshold count as holes; shallower
    # differences are sampling noise and are reported as nothing.
    thr = np.maximum(0.3, 0.02 * np.nan_to_num(out))
    mask = np.nan_to_num(out) - X > thr
    return out, mask


def _interp_nans(M, axis, circular):
    M = np.moveaxis(M.copy(), axis, 1)
    n = M.shape[1]
    idx = np.arange(n)
    for i in range(M.shape[0]):
        row = M[i]
        ok = np.isfinite(row)
        if ok.all() or not ok.any():
            continue
        if circular:
            row[~ok] = np.interp(idx[~ok], idx[ok], row[ok], period=n)
        else:
            row[~ok] = np.interp(idx[~ok], idx[ok], row[ok])
    return np.moveaxis(M, 1, axis)


# --------------------------------------------------------------------------
# Douglas-Peucker
# --------------------------------------------------------------------------

def dp_stations(z, R, tol, forced=()):
    """DP along the axis on whole r-vectors. Returns kept row indices.

    A station is dropped when every angle's radius is within tol of the
    straight line between its neighbours - so a plain cylinder keeps two
    rings and a cone keeps two rings.
    """
    n = len(z)
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    for f in forced:
        keep[f] = True
    idx = np.flatnonzero(keep)
    stack = list(zip(idx[:-1], idx[1:]))
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        t = (z[a + 1:b] - z[a]) / max(z[b] - z[a], 1e-12)
        interp = R[a] + t[:, None] * (R[b] - R[a])
        err = np.abs(R[a + 1:b] - interp).max(axis=1)
        m = int(err.argmax())
        if err[m] > tol:
            c = a + 1 + m
            keep[c] = True
            stack += [(a, c), (c, b)]
    return np.flatnonzero(keep)


def dp_polyline(P, tol):
    """Classic DP on an open 2D polyline. Returns kept indices."""
    n = len(P)
    if n < 3:
        return np.arange(n)
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        d = P[b] - P[a]
        L = np.hypot(*d)
        seg = P[a + 1:b] - P[a]
        if L < 1e-12:
            err = np.hypot(seg[:, 0], seg[:, 1])
        else:
            err = np.abs(seg[:, 0] * d[1] - seg[:, 1] * d[0]) / L
        m = int(err.argmax())
        if err[m] > tol:
            c = a + 1 + m
            keep[c] = True
            stack += [(a, c), (c, b)]
    return np.flatnonzero(keep)


def section_breaks(row, thetas, tol):
    """DP one closed cross-section as a 2D polyline, in two halves.

    DP on a closed curve from index 0 back to index 0 has a zero-length chord
    and bails, so it runs 0-180 and 180-360 separately.
    """
    A = len(row)
    P = np.column_stack([row * np.cos(thetas), row * np.sin(thetas)])
    h = A // 2
    first = dp_polyline(P[:h + 1], tol)
    second = dp_polyline(np.vstack([P[h:], P[:1]]), tol) + h
    idx = set(first.tolist()) | set(i % A for i in second.tolist())
    return idx


def count_regions(mask):
    """4-connected regions in a (station x angle) mask; angles wrap."""
    S, A = mask.shape
    seen = np.zeros_like(mask, bool)
    n = 0
    for i0, j0 in zip(*np.nonzero(mask)):
        if seen[i0, j0]:
            continue
        n += 1
        stack = [(i0, j0)]
        seen[i0, j0] = True
        while stack:
            i, j = stack.pop()
            for a, b in ((i + 1, j), (i - 1, j), (i, (j + 1) % A), (i, (j - 1) % A)):
                if 0 <= a < S and mask[a, b] and not seen[a, b]:
                    seen[a, b] = True
                    stack.append((a, b))
    return n


def is_round(row, tol):
    return float(np.nanmax(row) - np.nanmin(row)) <= tol


def choose_angles(rows, thetas, tol, segments):
    """Pick the ring angles (as dense-sample indices) shared by every ring.

    Round sections get the uniform segment count; non-round sections get
    their DP break points (hex corners, flats). Where both exist, a uniform
    angle that crowds a DP break point is dropped - the corner wins.
    """
    A = len(thetas)
    dp_idx = set()
    any_round = False
    for row in rows:
        if is_round(row, tol):
            any_round = True
        else:
            dp_idx |= section_breaks(row, thetas, tol)

    out = set(dp_idx)
    if any_round or not dp_idx:
        step = A / float(segments)
        uni = {int(round(i * step)) % A for i in range(segments)}
        guard = max(1, int(step / 3))
        for u in uni:
            if all(min((u - d) % A, (d - u) % A) > guard for d in dp_idx):
                out.add(u)

    # Merge DP points closer than two samples (noise on a flat).
    merged = []
    for i in sorted(out):
        if merged and i - merged[-1] <= 1:
            continue
        merged.append(i)
    if len(merged) > 2 and (merged[0] + A - merged[-1]) <= 1:
        merged.pop()
    return np.array(merged, dtype=int)


def segments_for(radius_mm, seg_small, seg_medium, seg_large):
    """The size table, low end first."""
    if radius_mm < 10.0:
        return seg_small
    if radius_mm < 60.0:
        return seg_medium
    return seg_large


def budget_for(radius_mm):
    if radius_mm < 10.0:
        return 150
    if radius_mm < 40.0:
        return 400
    return 2000
