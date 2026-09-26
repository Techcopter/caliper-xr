"""Exact point-to-surface distance on sheet-metal slivers.

mathutils' BVH works in single precision and misplaces the nearest point on
long thin triangles; caliper_xr.nearest must agree with a double-precision
brute force wherever a decision depends on it.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
from mathutils import Vector  # noqa: E402
from harness import MM, check, done  # noqa: E402

from caliper_xr import nearest  # noqa: E402


def brute(P, V, T):
    A, B, C = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    return np.array([nearest.point_tri(p, A, B, C).min() for p in P])


rng = np.random.default_rng(7)
worst_bvh = worst_exact = 0.0
for trial in range(40):
    # a strip 600-900 mm long and 1.5-4 mm wide, somewhere in a 1 m cell,
    # tilted a little: two triangles, as a CAD exporter writes a flange
    L, W = rng.uniform(600, 900), rng.uniform(1.5, 4.0)
    o = rng.uniform(-500, 500, 3)
    u = rng.normal(size=3)
    u /= np.linalg.norm(u)
    v = np.cross(u, rng.normal(size=3))
    v /= np.linalg.norm(v)
    V = np.array([o, o + u * L, o + u * L + v * W, o + v * W]) * MM
    T = np.array([[0, 1, 2], [0, 2, 3]])
    X = V[T]
    P = np.vstack([X.mean(1),                                           # on the strip
                   X.mean(1) + np.cross(u, v) * rng.uniform(0.05, 0.5) * MM,  # just off it
                   V])
    surf = nearest.Surface(V, T, exact_above=0.0)
    ref = brute(P, V, T)
    got = surf(P)
    bvh = np.array([surf.tree.find_nearest(Vector(p))[3] for p in P])
    worst_exact = max(worst_exact, float(np.abs(got - ref).max()))
    worst_bvh = max(worst_bvh, float(np.abs(bvh - ref).max()))

print(f"  info  mathutils BVH off by up to {worst_bvh / MM:.4f} mm on these strips")
check(worst_exact / MM < 1e-6, f"exact distance matches brute force (max error {worst_exact / MM:.2e} mm)")

# at or below exact_above a value may be an upper bound, never an underestimate
V = np.array([[0, 0, 0], [800, 0, 0], [800, 3, 0.8], [0, 3, 0.8]]) * MM
T = np.array([[0, 1, 2], [0, 2, 3]])
P = np.vstack([V[T].mean(1), rng.uniform([0, -5, -5], [800, 8, 5], (400, 3)) * MM])
lim = 0.1 * MM
got = nearest.Surface(V, T, exact_above=lim)(P)
ref = brute(P, V, T)
check(bool(np.all(got >= ref - 1e-12)), "never below the true distance")
over = ref > lim
check(bool(np.allclose(got[over], ref[over], atol=1e-12)), "exact above the threshold")
check(bool(np.all(got[~over] <= lim + 1e-12)), "within the threshold below it")
check(float(got[:2].max()) / MM < 1e-6, "an exact copy measures zero, not a third of a millimetre")
done("nearest")
