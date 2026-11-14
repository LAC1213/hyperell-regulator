"""Euclidean minimum spanning trees and their oriented slit boundaries."""
import numpy as np


def minimum_spanning_tree(points):
    """Deterministic Kruskal on the complete Euclidean graph of marked points."""
    z = np.asarray(points, complex)
    if len(set(z)) != len(z):
        raise ValueError('cut vertices cannot be distinguished in float64 geometry')
    parent = list(range(len(z)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    edges = []
    candidates = sorted((abs(z[i]-z[j]), i, j)
                        for i in range(len(z)) for j in range(i))
    for _, i, j in candidates:
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b
            edges.append((i, j))
    return edges


def orient_to_root(size, edges, root):
    """Orient every finite edge towards the vertex carrying the infinity ray."""
    neighbors = [[] for _ in range(size)]
    for a, b in edges:
        neighbors[a].append(b)
        neighbors[b].append(a)
    parents = {root: None}
    queue = [root]
    for a in queue:
        for b in sorted(neighbors[a]):
            if b not in parents:
                parents[b] = a
                queue.append(b)
    if len(parents) != size or len(edges) != size-1:
        raise ValueError('cut edges must form a spanning tree')
    return [(a, parents[a]) for a in range(size) if a != root]+[(root, None)]


def boundary_word(points, edges, direction):
    """Walk around an embedded tree using clockwise incident-edge order.

    Each directed bank occurs exactly once. The final virtual vertex is
    infinity; its two banks are joined by the exterior circle, not a cap.
    """
    z = list(map(complex, points))
    root = edges[-1][0]
    infinity = len(z)
    scale = max(abs(a-b) for a in z for b in z) or 1
    z.append(z[root]+60*scale*direction)
    neighbors = [[] for _ in z]
    lookup = {}
    for i, (a, b) in enumerate(edges):
        b = infinity if b is None else b
        neighbors[a].append(b)
        neighbors[b].append(a)
        lookup[a, b] = i, 1
        lookup[b, a] = i, -1
    for a, row in enumerate(neighbors):
        row.sort(key=lambda b: np.angle(z[b]-z[a]))
    leaf = next(a for a in range(infinity) if len(neighbors[a]) == 1)
    first = (leaf, neighbors[leaf][0])
    word, banks = [], []
    a, b = first
    while True:
        word.append(lookup[a, b])
        banks.append((a, b))
        row = neighbors[b]
        c = row[(row.index(a)-1) % len(row)]
        a, b = b, c
        if (a, b) == first:
            break
        if len(word) > 2*len(edges):
            raise ValueError('tree boundary does not close')
    if len(word) != 2*len(edges):
        raise ValueError('tree boundary does not visit every bank')
    return np.asarray(z), neighbors, word, banks


def tree_contour(points, edges, direction, references, eps, span, dens=30):
    """A narrow contour with a sample on each bank of every tree edge."""
    z, neighbors, word, banks = boundary_word(points, edges, direction)
    infinity = len(points)
    root = edges[-1][0]
    def line(path, end, count):
        path.extend(path[-1]+(end-path[-1])*np.linspace(0, 1, count)[1:])
    def arc(path, center, radius, sweep, count):
        angle = np.angle(path[-1]-center)
        path.extend(center+radius*np.exp(1j*(angle+sweep*np.linspace(0, 1, count)[1:])))
    a, b = banks[0]
    tangent = (z[b]-z[a])/abs(z[b]-z[a])
    path, marks = [z[a]+1j*eps*tangent], {}
    for k, ((a, b), (edge, side)) in enumerate(zip(banks, word)):
        tangent = (z[b]-z[a])/abs(z[b]-z[a])
        mid = complex(references[edge])+1j*eps*tangent
        line(path, mid, max(8, int(dens*abs(mid-path[-1])/span)))
        marks[edge, side] = len(path)-1
        end = z[b]+1j*eps*tangent
        line(path, end, max(8, int(dens*abs(end-path[-1])/span)))
        _, c = banks[(k+1) % len(banks)]
        following = (z[c]-z[b])/abs(z[c]-z[b])
        if b == infinity:
            radius = abs(path[-1]-z[root])
            target = z[b]+1j*eps*following
            sweep = np.mod(np.angle(target-z[root])-np.angle(path[-1]-z[root]), 2*np.pi)
            arc(path, z[root], radius, sweep, 900)
        else:
            gap = np.mod(np.angle(z[a]-z[b])-np.angle(z[c]-z[b]), 2*np.pi)
            if len(neighbors[b]) == 1:
                gap = 2*np.pi
            turn = np.pi-gap
            if abs(turn) < 1e-12:
                continue
            if turn < 0:
                arc(path, z[b], eps, turn, 60)
            else:
                # Intersection of the two offset banks in the exterior wedge.
                rhs = 1j*eps*(following-tangent)
                A = np.array([[tangent.real, -following.real],
                              [tangent.imag, -following.imag]])
                t = np.linalg.solve(A, [rhs.real, rhs.imag])[0]
                line(path, z[b]+1j*eps*tangent+t*tangent, 6)
    if abs(path[-1]-path[0]) > 1e-10*span:
        raise ValueError('tree contour does not close')
    return np.asarray(path), marks, word
