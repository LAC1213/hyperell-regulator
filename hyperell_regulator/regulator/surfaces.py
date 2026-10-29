"""Triangulated, consistently labelled sheets of the slit phi-plane.

Display geometry is independent of the saved quadrature. Triangulation edges
recover the entire cut, including its detours and the ray to the window edge.
Labels are continued across triangle edges only in the slit complement and
checked against both banks of the integrator's cut.
"""

from collections import deque

import numpy as np

from .cover import Cover
from .explorer import parse_request
from .stokes import cut_data, cut_scale, follow, racetrack


def _continue(cov, start, end, x, y):
    try:
        return follow(cov, start, end, x, y)
    except RuntimeError:
        from ._legs import leg_in_w
        xx, yy, _, _ = leg_in_w(cov, start, end, x, y,
                                np.zeros(len(x), complex), integrate=False)
        return xx, yy


def surface_mesh(data, density=20):
    """Return one finite-window mesh per globally labelled sheet.

    Faces approach their vertices by a small relative inset, avoiding exact
    ramification points and finite-value poles of the display projection.
    All vertices are actual continued points of the curve, not interpolated
    heights or independently sorted fibres.
    """
    import triangle

    F, phi, _ = parse_request(data["input"])
    cov = Cover(F, *phi)
    nodes, kinds, direction = cov.cut_path()
    upper, lower, _, _ = cut_data(cov, nodes, direction)
    integration_direction = direction
    # Cardinal rays inherit tiny trig roundoff. Qhull must see exactly
    # collinear real nodes, rather than an artificial 1e-15-wide sliver.
    direction = complex(0 if abs(direction.real) < 1e-12 else direction.real,
                        0 if abs(direction.imag) < 1e-12 else direction.imag)
    critical = nodes[np.array(kinds) != "none"]
    base = np.array([complex(*p) for edge in data["edges"] for p in edge["base"]])
    seed = racetrack(nodes, direction)[0][0]
    base = np.r_[base, seed]
    span = max(np.ptp(base.real), np.ptp(base.imag))
    pad = span*.13
    bounds = [base.real.min()-pad, base.real.max()+pad,
              base.imag.min()-pad, base.imag.max()+pad]
    center = complex((bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2)
    distances = []
    for value, component, low, high in ((nodes[-1].real, direction.real, bounds[0], bounds[1]),
                                       (nodes[-1].imag, direction.imag, bounds[2], bounds[3])):
        if abs(component) > 1e-12:
            distances.append(((high if component > 0 else low)-value)/component)
    end = nodes[-1]+min(v for v in distances if v > 0)*direction
    cut = np.r_[nodes, end]
    gridx, gridy = np.meshgrid(np.linspace(bounds[0], bounds[1], density),
                              np.linspace(bounds[2], bounds[3], density))
    sites = list((gridx+1j*gridy).ravel())
    for z in critical:
        gap = min(abs(z-v) for v in critical if abs(z-v) > span*1e-13)
        for fraction in (.025, .12, .35):
            sites.extend(z+fraction*gap*np.exp(1j*np.linspace(0, 2*np.pi, 17)[:-1]))
    points, lookup = [], {}

    def add(z):
        q = (z-center)/span
        key = (round(q.real, 14), round(q.imag, 14))
        if key not in lookup:
            lookup[key] = len(points)
            points.append([q.real, q.imag])
        return lookup[key]

    for z in sites:
        on_cut = False
        for a, b in zip(cut[:-1], cut[1:]):
            t = np.clip(((z-a)/(b-a)).real, 0, 1)
            if abs(z-(a+t*(b-a))) < span*1e-13:
                on_cut = True
                break
        if not on_cut:
            add(z)
    constraints = []
    step = span/density/3
    for a, b in zip(cut[:-1], cut[1:]):
        samples = np.linspace(a, b, max(2, int(np.ceil(abs(b-a)/step))+1))
        constraints.extend((add(a), add(b)) for a, b in zip(samples[:-1], samples[1:]))
    # Constrained triangulation preserves tiny critical-value clusters that
    # Qhull can discard as coplanar. Its robust predicates enforce the slit
    # directly; no repeated edge bisection or coordinate perturbation.
    triangulated = triangle.triangulate(dict(vertices=np.array(points), segments=np.array(constraints)), "pcnQ")
    normalized = triangulated["vertices"]
    faces, neighbors = triangulated["triangles"], triangulated["neighbors"]
    blocked = {tuple(sorted(pair)) for pair in triangulated["segments"]}
    coordinates = center+span*(normalized[:, 0]+1j*normalized[:, 1])
    centers = coordinates[faces].mean(axis=1)
    origins = normalized[faces[:, 0]]
    ab = normalized[faces[:, 1]]-origins
    ac = normalized[faces[:, 2]]-origins
    determinants = ab[:, 0]*ac[:, 1]-ab[:, 1]*ac[:, 0]
    if np.any(determinants == 0):
        raise RuntimeError("A surface triangle has zero area")

    def face_at(z):
        q = (z-center)/span
        delta = np.array([q.real, q.imag])-origins
        s = (delta[:, 0]*ac[:, 1]-delta[:, 1]*ac[:, 0])/determinants
        t = (ab[:, 0]*delta[:, 1]-ab[:, 1]*delta[:, 0])/determinants
        matches = np.flatnonzero((s >= -1e-10) & (t >= -1e-10) & (s+t <= 1+1e-10))
        return int(matches[0]) if len(matches) else -1

    first = face_at(seed)
    if first < 0:
        raise RuntimeError("Sheet-label seed is outside the surface window")
    x, y = cov.fibre(seed)
    labels = {first: _continue(cov, seed, centers[first], x, y)}
    queue = deque([first])
    cycle_error = 0.0
    cycle_count = 0
    while queue:
        index = queue.popleft()
        x, y = labels[index]
        face = faces[index]
        for opposite, neighbor in enumerate(neighbors[index]):
            if neighbor < 0:
                continue
            edge = tuple(sorted(int(face[k]) for k in range(3) if k != opposite))
            if edge in blocked or (neighbor in labels and cycle_count >= 80):
                continue
            midpoint = (coordinates[edge[0]]+coordinates[edge[1]])/2
            xx, yy = _continue(cov, centers[index], midpoint, x, y)
            xx, yy = _continue(cov, midpoint, centers[neighbor], xx, yy)
            if neighbor in labels:
                ex, ey = labels[neighbor]
                defect = np.max((abs(xx-ex)+abs(yy-ey))/(1+abs(ex)+abs(ey)))
                cycle_error = max(cycle_error, float(defect))
                cycle_count += 1
            else:
                labels[int(neighbor)] = (xx, yy)
                queue.append(int(neighbor))
    if len(labels) != len(faces) or cycle_error > 1e-6:
        raise RuntimeError("Inconsistent sheet continuation over the slit mesh: %g" % cycle_error)
    boundary_error = 0.0
    for j, a in enumerate(nodes):
        b = nodes[j+1] if j+1 < len(nodes) else data["edges"][-1]["base"][-1]
        b = complex(*b) if isinstance(b, list) else complex(b)
        w = ((a+b)/2 if j+1 < len(nodes)
             else a+cut_scale(nodes)[0]*integration_direction)
        normal = 1j*(b-a)/abs(b-a)
        radius = min(abs(w-z) for z in critical)*1e-5
        x, y = cov.fibre(w)
        for sign, permutation in ((1, upper[j]), (-1, lower[j])):
            target = w+sign*radius*normal
            index = face_at(target)
            if index < 0:
                raise RuntimeError("Cut bank is outside the mesh")
            rx, ry = labels[index]
            mx, my = _continue(cov, centers[index], target, rx, ry)
            ex, ey = _continue(cov, w, target, x, y)
            ex, ey = ex[permutation], ey[permutation]
            defect = np.max((abs(mx-ex)+abs(my-ey))/(1+abs(ex)+abs(ey)))
            boundary_error = max(boundary_error, float(defect))
    if boundary_error > 1e-6:
        raise RuntimeError("Mesh sheets disagree with the cut's bank labels: %g" % boundary_error)
    vertices = [[] for _ in range(cov.N)]
    triangles = []
    curve_error = function_error = 0.0
    for index, face in enumerate(faces):
        x, y = labels[index]
        projected = []
        for target in coordinates[face]:
            # Separate cut banks and avoid exactly singular fibres. These
            # subpixel insets also keep finite-value infinities out of GL.
            target += (centers[index]-target)*1e-6
            xx, yy = _continue(cov, centers[index], target, x, y)
            fx = np.polyval(cov.F, xx)
            ce = float(np.max(abs(yy*yy-fx)/(1+abs(fx))))
            fe = float(np.max(abs(cov.phi_stable(xx, yy)-target)/(1+abs(target))))
            curve_error, function_error = max(curve_error, ce), max(function_error, fe)
            if not np.isfinite(ce+fe) or max(ce, fe) > 1e-6:
                raise RuntimeError("Surface vertex residual is too large (%g, %g)" % (ce, fe))
            height = np.arcsinh((np.exp(-.55j)*(yy if cov.pure_x else xx)).real)
            projected.append([[float(target.real), float(target.imag), float(h)] for h in height])
        start = len(vertices[0])
        triangles.append([start, start+1, start+2])
        for sheet in range(cov.N):
            vertices[sheet].extend(p[sheet] for p in projected)
    return dict(schema_version=1, degree=cov.N, bounds=[float(v) for v in bounds],
                triangles=triangles, sheets=[dict(sheet=i+1, vertices=v) for i, v in enumerate(vertices)],
                validation=dict(triangles_per_sheet=len(triangles), checked_continuation_cycles=cycle_count,
                                max_relative_sheet_cycle_defect=cycle_error,
                                max_relative_sheet_boundary_defect=boundary_error,
                                max_relative_curve_residual=curve_error,
                                max_relative_function_residual=function_error))
