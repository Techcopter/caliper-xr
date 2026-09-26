"""Hidden geometry: CAD faces no one can see from outside the selection.

Primary method: the whole selection is drawn from many directions around a
sphere on the GPU, every triangle in its own ID colour, both sides, no
culling. A triangle that shows up in any view is visible; a CAD face is
hidden only if none of its triangles ever shows up. This sees through vents,
gaps and grazing angles down to about a pixel (about 1 mm on a 1 m part).

Fallback (no GPU available): dense ray sampling. It is slower and less
certain, so the report says which method ran.
"""

import math

import numpy as np
from mathutils import Matrix, Vector
from mathutils.bvhtree import BVHTree


_GPU_READY = False


def fibonacci(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)])


def world_tris(objs, depsgraph):
    """World-space triangles of every object, in loop-triangle order."""
    out = []
    for ob in objs:
        ev = ob.evaluated_get(depsgraph)
        me = ev.to_mesh()
        try:
            co = np.empty(len(me.vertices) * 3)
            me.vertices.foreach_get("co", co)
            co = co.reshape(-1, 3)
            M = np.array(ob.matrix_world)
            W = co @ M[:3, :3].T + M[:3, 3]
            lt = me.loop_triangles
            T = np.empty(len(lt) * 3, np.int64)
            lt.foreach_get("vertices", T)
            out.append(W[T.reshape(-1, 3)])
        finally:
            ev.to_mesh_clear()
    return out


def _look(eye, target, up):
    z = eye - target
    z /= np.linalg.norm(z)
    x = np.cross(up, z)
    if np.linalg.norm(x) < 1e-6:
        x = np.cross(np.array([0.0, 1.0, 0.0]), z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    R = np.array([x, y, z])
    V = np.eye(4)
    V[:3, :3] = R
    V[:3, 3] = -R @ eye
    return Matrix(V.tolist())


def _ortho(r, near, far):
    P = np.eye(4)
    P[0, 0] = P[1, 1] = 1.0 / r
    P[2, 2] = -2.0 / (far - near)
    P[2, 3] = -(far + near) / (far - near)
    return Matrix(P.tolist())


def visible_gpu(objs, depsgraph, n_views=128, res=1536):
    """Per object, a bool per loop triangle: seen in at least one view.
    Returns None when no GPU is available."""
    try:
        import bpy
        import gpu
        from gpu_extras.batch import batch_for_shader
        global _GPU_READY
        if bpy.app.background and hasattr(gpu, "init") and not _GPU_READY:
            gpu.init()                     # Blender 5: GPU drawing from the command line
            _GPU_READY = True
        off = gpu.types.GPUOffScreen(res, res)
    except Exception:
        return None
    tris = world_tris(objs, depsgraph)
    counts = [len(t) for t in tris]
    allt = np.concatenate(tris) if tris else np.zeros((0, 3, 3))
    M = len(allt)
    ids = np.arange(1, M + 1, dtype=np.int64)
    rgb = np.stack([ids & 255, (ids >> 8) & 255, (ids >> 16) & 255], 1) / 255.0
    col = np.repeat(np.column_stack([rgb, np.ones(M)]), 3, axis=0).astype(np.float32)
    pos = allt.reshape(-1, 3).astype(np.float32)
    lo, hi = pos.min(0), pos.max(0)
    c = (lo + hi) / 2
    r = float(np.linalg.norm(hi - lo)) / 2 * 1.02 or 1.0
    shader = gpu.shader.from_builtin('FLAT_COLOR')
    batch = batch_for_shader(shader, 'TRIS', {"pos": pos, "color": col})
    seen = np.zeros(M + 1, bool)
    up = np.array([0.0, 0.0, 1.0])
    proj = _ortho(r, 0.01 * r, 4.0 * r)
    try:
        with off.bind():
            fb = gpu.state.active_framebuffer_get()
            gpu.state.depth_test_set('LESS_EQUAL')
            gpu.state.depth_mask_set(True)
            gpu.state.face_culling_set('NONE')
            for d in fibonacci(n_views):
                fb.clear(color=(0.0, 0.0, 0.0, 0.0), depth=1.0)
                with gpu.matrix.push_pop():
                    gpu.matrix.load_matrix(_look(c + d * 2.0 * r, c, up))
                    gpu.matrix.load_projection_matrix(proj)
                    batch.draw(shader)
                buf = fb.read_color(0, 0, res, res, 4, 0, 'UBYTE')
                px = np.frombuffer(buf, dtype=np.uint8).reshape(-1, 4).astype(np.int64)
                idx = px[:, 0] | (px[:, 1] << 8) | (px[:, 2] << 16)
                seen[idx] = True
            gpu.state.depth_test_set('NONE')
    finally:
        off.free()
    seen = seen[1:]
    out, s = {}, 0
    for ob, n in zip(objs, counts):
        out[ob.name] = seen[s:s + n]
        s += n
    return out


def occluder(objs, depsgraph):
    """One BVH over all objs in world space (fallback method)."""
    tris = world_tris(objs, depsgraph)
    allt = np.concatenate(tris)
    V = allt.reshape(-1, 3)
    T = np.arange(len(V)).reshape(-1, 3)
    diag = float(np.linalg.norm(V.max(0) - V.min(0))) or 1.0
    return BVHTree.FromPolygons(V.tolist(), T.tolist()), diag


def hidden_from_faces(patch, visible):
    """A patch is hidden only if none of its triangles was ever seen."""
    npatch = int(patch.max()) + 1
    seen = np.zeros(npatch, bool)
    np.logical_or.at(seen, patch, visible)
    return ~seen


def hidden_patches(cm, patch, matrix_world, tree, diag, n_faces=48, n_dirs=256, seed=0):
    """Fallback: dense ray sampling. True per patch when no sample is seen."""
    Mw = np.array(matrix_world)
    R, t = Mw[:3, :3], Mw[:3, 3]
    Rn = np.linalg.inv(R).T
    D = fibonacci(n_dirs)
    far = diag * 10.0
    eps = diag * 2e-6
    rng = np.random.default_rng(seed)
    npatch = int(patch.max()) + 1
    order = np.argsort(patch, kind="stable")
    bounds = np.searchsorted(patch[order], np.arange(npatch + 1))
    hidden = np.zeros(npatch, bool)
    bary = np.array([[1 / 3, 1 / 3, 1 / 3], [0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]])
    for p in range(npatch):
        faces = order[bounds[p]:bounds[p + 1]]
        if len(faces) > n_faces:
            big = faces[np.argsort(-cm.A[faces])[:8]]
            rest = rng.choice(faces, n_faces - 8, replace=False)
            faces = np.unique(np.concatenate([big, rest]))
        seen = False
        for f in faces:
            tri = cm.V[cm.T[f]]
            n = Rn @ cm.N[f]
            n /= max(np.linalg.norm(n), 1e-30)
            dirs = D[np.argsort(-np.abs(D @ n))]
            for w in bary:
                p0 = (w @ tri) @ R.T + t
                for d in dirs:
                    side = 1.0 if d @ n >= 0 else -1.0
                    o = p0 + n * (eps * side) + d * eps
                    if tree.ray_cast(Vector(o), Vector(d), far)[0] is None:
                        seen = True
                        break
                if seen:
                    break
            if seen:
                break
        hidden[p] = not seen
    return hidden
