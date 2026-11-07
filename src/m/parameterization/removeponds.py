import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from reinitializelevelset import reinitializelevelset

def removeponds(md, min_vertices=3, min_area=1e2, raise_bed=1,
                       keep_boundary_connected=True, verbose=True):
    """Remove small isolated ocean regions (ponds) from md.mask.ocean_levelset.

    A connected ocean region is a pond if it has fewer than min_vertices
    vertices (or an area below min_area [m^2], if given) and, when
    keep_boundary_connected is True, does not touch the domain boundary.

    At removed vertices the levelset is flipped to grounded and reinitialized,
    base is raised by raise_bed [m], bed is set equal to base, and thickness
    is recomputed as surface - base (so it drops by raise_bed).

    """
    ls = np.asarray(md.mask.ocean_levelset, dtype=float).flatten().copy()
    nv = ls.size

    elems = np.asarray(md.mesh.elements).astype(int) - 1      # 1-based -> 0-based
    x = np.asarray(md.mesh.x, dtype=float).flatten()
    y = np.asarray(md.mesh.y, dtype=float).flatten()
    onbound = np.asarray(md.mesh.vertexonboundary).flatten().astype(bool)

    ocean = ls < 0

    # edges of all triangles; keep only those with both ends in the ocean
    edges = np.vstack((elems[:, [0, 1]], elems[:, [1, 2]], elems[:, [2, 0]]))
    edges = edges[ocean[edges[:, 0]] & ocean[edges[:, 1]]]

    graph = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                       shape=(nv, nv))
    ncomp, labels = connected_components(graph, directed=False)

    # per-component statistics (non-ocean vertices are their own singleton
    # components, but are excluded by the ocean mask below)
    count = np.bincount(labels[ocean], minlength=ncomp)

    # lumped vertex areas (each triangle contributes 1/3 of its area to each vertex)
    x1, x2, x3 = x[elems[:, 0]], x[elems[:, 1]], x[elems[:, 2]]
    y1, y2, y3 = y[elems[:, 0]], y[elems[:, 1]], y[elems[:, 2]]
    elem_area = 0.5 * np.abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    vert_area = np.bincount(elems.ravel(), weights=np.repeat(elem_area / 3.0, 3),
                            minlength=nv)
    area = np.bincount(labels[ocean], weights=vert_area[ocean], minlength=ncomp)

    touches_boundary = np.zeros(ncomp, dtype=bool)
    touches_boundary[labels[ocean & onbound]] = True

    # decide which components are ponds
    small = count < min_vertices
    if min_area is not None:
        small |= area < min_area
    pond = small.copy()
    if keep_boundary_connected:
        pond &= ~touches_boundary

    remove = ocean & pond[labels]
    ls[remove] = np.abs(ls[remove])

    if verbose:
        n_ponds = np.unique(labels[remove]).size
        print(f"      → Removed {n_ponds} ocean ponds ({remove.sum()} vertices)")

    if remove.any():
        md.mask.ocean_levelset = ls.reshape(np.shape(md.mask.ocean_levelset))
        md.mask.ocean_levelset = reinitializelevelset(md,md.mask.ocean_levelset)
        md.geometry.base[remove] += raise_bed
        md.geometry.surface[remove] += raise_bed # preserve thickness
        md.geometry.bed[remove] = md.geometry.base[remove] 
    return md