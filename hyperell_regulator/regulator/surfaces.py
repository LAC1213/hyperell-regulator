"""Two slit-sheet display meshes from explicit square-root factors.

The constrained triangulation stays independent of numerical integration.
Sheet labels are checked around mesh cycles and against both cut banks.
"""
from collections import deque
import numpy as np
from .explorer import parse_request
from .stokes import Cut, cut_scale

def _bank_radius(w, critical, cut, index):
    """Stay on the chosen cut bank, including between nearby parallel slits."""
    radius = min(abs(w-z) for z in critical)*1e-5
    for other, (left, right) in enumerate(zip(cut[:-1], cut[1:])):
        if other == index:
            continue
        fraction = np.clip(((w-left)/(right-left)).real, 0, 1)
        separation = abs(w-(left+fraction*(right-left)))
        radius = min(radius, separation*.1)
    return radius


def _ray_endpoint(start, direction, bounds):
    """End the slit exactly on the window boundary, without a rounding gap."""
    crossings = []
    for axis, (value, component, low, high) in enumerate(
            ((start.real, direction.real, bounds[0], bounds[1]),
             (start.imag, direction.imag, bounds[2], bounds[3]))):
        if abs(component) > 1e-12:
            wall = high if component > 0 else low
            reach = (wall-value)/component
            if reach > 0:
                crossings.append((reach, axis, wall))
    reach, axis, wall = min(crossings)
    end = start+reach*direction
    return complex(wall, end.imag) if axis == 0 else complex(end.real, wall)


def _label_seed(nodes):
    """The integrator's first racetrack point, without constructing its path."""
    _, gap = cut_scale(nodes)
    tangent = (nodes[1]-nodes[0])/abs(nodes[1]-nodes[0])
    return nodes[0]+.02j*gap*tangent


def surface_mesh(data, density=20):
    import triangle

    F, centers, prec = parse_request(data['input'])
    cut_data = Cut(F, centers, prec)
    nodes = np.array([complex(z) for z in cut_data.nodes])
    direction = complex(cut_data.direction)
    upper, lower = cut_data.upper, cut_data.lower
    integration_direction = direction
    direction = complex(0 if abs(direction.real) < 1e-12 else direction.real,
                        0 if abs(direction.imag) < 1e-12 else direction.imag)
    critical = np.array([complex(r) for r in cut_data.roots])
    leading = complex(F.leading_coefficient())
    def continue_y(start, end, y):
        return y*np.prod(np.sqrt((end-critical)/(start-critical)))
    base = np.array([complex(*p) for edge in data["edges"] for p in edge["base"]])
    seed = _label_seed(nodes)
    base = np.r_[base, seed]
    span = max(np.ptp(base.real), np.ptp(base.imag))
    pad = span*.13
    bounds = [base.real.min()-pad, base.real.max()+pad,
              base.imag.min()-pad, base.imag.max()+pad]
    center = complex((bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2)
    end = _ray_endpoint(nodes[-1], direction, bounds)
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
    y = np.sqrt(leading)*np.prod(np.sqrt(seed-critical))
    labels = {first: continue_y(seed, centers[first], y)}
    queue = deque([first])
    cycle_error = 0.0
    cycle_count = 0
    visited_edges = set()
    while queue:
        index = queue.popleft()
        y = labels[index]
        face = faces[index]
        for opposite, neighbor in enumerate(neighbors[index]):
            if neighbor < 0:
                continue
            edge = tuple(sorted(int(face[k]) for k in range(3) if k != opposite))
            pair = tuple(sorted((index, int(neighbor))))
            if edge in blocked or pair in visited_edges:
                continue
            visited_edges.add(pair)
            midpoint = (coordinates[edge[0]]+coordinates[edge[1]])/2
            yy = continue_y(centers[index], midpoint, y)
            yy = continue_y(midpoint, centers[neighbor], yy)
            if neighbor in labels:
                ey = labels[neighbor]
                defect = abs(yy-ey)/(1+abs(ey))
                cycle_error = max(cycle_error, float(defect))
                cycle_count += 1
                if defect > 1e-6:
                    raise RuntimeError("Inconsistent sheet continuation over the slit mesh: %g" % defect)
            else:
                labels[int(neighbor)] = yy
                queue.append(int(neighbor))
    if len(labels) != len(faces):
        raise RuntimeError("Cut mesh has disconnected regions at floating point precision")
    if cycle_error > 1e-6:
        raise RuntimeError("Inconsistent sheet continuation over the slit mesh: %g" % cycle_error)
    boundary_error = 0.0
    for j, a in enumerate(nodes):
        b = nodes[j+1] if j+1 < len(nodes) else data["edges"][-1]["base"][-1]
        b = complex(*b) if isinstance(b, list) else complex(b)
        w = ((a+b)/2 if j+1 < len(nodes)
             else a+cut_scale(nodes)[0]*integration_direction)
        normal = 1j*(b-a)/abs(b-a)
        radius = _bank_radius(w, critical, cut, j)
        y = complex(cut_data.reference_y[j])
        for sign, permutation in ((1, upper[j]), (-1, lower[j])):
            target = w+sign*radius*normal
            index = face_at(target)
            if index < 0:
                raise RuntimeError("Cut bank is outside the mesh")
            my = continue_y(centers[index], target, labels[index])
            ey = continue_y(w, target, y)*(1 if permutation[0] == 0 else -1)
            defect = abs(my-ey)/(1+abs(ey))
            boundary_error = max(boundary_error, float(defect))
    if boundary_error > 1e-6:
        raise RuntimeError("Mesh sheets disagree with the cut's bank labels: %g" % boundary_error)
    vertices = [[] for _ in range(2)]
    triangles = []
    curve_error = 0.0
    for index, face in enumerate(faces):
        y = labels[index]
        projected = []
        for target in coordinates[face]:
            # Separate cut banks and avoid exactly singular fibres. These
            # subpixel insets also keep finite-value infinities out of GL.
            target += (centers[index]-target)*1e-6
            yy = continue_y(centers[index], target, y)
            fx = leading*np.prod(target-critical)
            ce = float(abs(yy*yy-fx)/(1+abs(fx)))
            curve_error = max(curve_error, ce)
            if not np.isfinite(ce) or ce > 1e-6:
                raise RuntimeError('Surface curve residual is too large')
            height = float(np.arcsinh((np.exp(-.55j)*yy).real))
            projected.append([[float(target.real), float(target.imag), height],
                              [float(target.real), float(target.imag), -height]])
        start = len(vertices[0])
        triangles.append([start, start+1, start+2])
        for sheet in range(2):
            vertices[sheet].extend(p[sheet] for p in projected)
    return dict(schema_version=1, degree=2, bounds=[float(v) for v in bounds],
                triangles=triangles, sheets=[dict(sheet=i+1, vertices=v) for i, v in enumerate(vertices)],
                validation=dict(triangles_per_sheet=len(triangles), checked_continuation_cycles=cycle_count,
                                max_relative_sheet_cycle_defect=cycle_error,
                                max_relative_sheet_boundary_defect=boundary_error,
                                max_relative_curve_residual=curve_error))
