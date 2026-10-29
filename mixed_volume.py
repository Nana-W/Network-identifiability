"""Mixed volume of the eliminated form of the paper's example.

Inclusion-exclusion over Minkowski sums; the normalization is the one
of the paper (n generic affine equations have mixed volume 1), and the
sanity check below confirms it. Result: MV = 2 = d.
"""
import itertools
import numpy as np
from scipy.spatial import ConvexHull

def vol(points, n):
    P = np.array(points, float)
    if len(P) < n + 1:
        return 0.0
    try:
        return ConvexHull(P).volume
    except Exception:          # lower-dimensional
        return 0.0

def minkowski(supports, ids, n):
    pts = [tuple([0] * n)]
    for i in ids:
        pts = list({tuple(a + b for a, b in zip(p, q))
                    for p in pts for q in supports[i]})
    return pts

def mixed_volume(supports):
    n = len(supports)
    total = 0.0
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            total += (-1) ** (n - k) * vol(minkowski(supports, S, n), n)
    return total

# supports of h(theta) - y*, variables (x, y, w, z)
example = [
    [(0,0,0,0), (1,0,0,0), (0,1,0,0), (1,1,0,0)],            # x+y+xy
    [(0,0,0,0), (1,0,0,0), (1,1,0,0)],                        # 2x+xy
    [(0,0,0,0), (1,0,0,0), (1,1,0,0), (0,0,1,0)],             # x+xy+w
    [(0,0,0,0), (1,0,0,0), (0,1,0,0), (1,1,0,0), (0,0,0,1)],  # 2x+y+xy+z
]
print("mixed volume of the eliminated form =", round(mixed_volume(example), 6))

simplex = [[(0,0,0,0), (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]] * 4
print("normalization check, 4 standard simplices =",
      round(mixed_volume(simplex), 6), "(must be 1)")
