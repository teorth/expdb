
# Some code for benchmarking frequently used Polytope functions
import parent 

import itertools
from polytope import Polytope
from region import Region
import random as rd
import time

rd.seed(1007)

def benchmark_empty(N, dim):

    # Emptyness check for non-empty polytopes
    polys = []
    for i in range(N):
        verts = [
            [rd.randint(0, 100) for j in range(dim + 1)] 
            for k in range(dim * 2)
        ]
        polys.append(Polytope.from_V_rep(verts))
    
    start = time.perf_counter()
    for i in range(N):
        polys[i].is_empty(include_boundary=False)
    end = time.perf_counter()

    print(f"{N} nonempty polytopes (dimension {dim}) emptiness check took {end - start} sec.")


def benchmark_union(N):

    # 5 dimensional boxes 
    polys = []
    for combo in itertools.product(range(N), repeat=5):
        b = [(c, c + 1) for c in combo]
        polys.append(Polytope.rect(b[0], b[1], b[2], b[3], b[4]))
    
    # Shuffle the list 
    rd.shuffle(polys)

    region = Region.disjoint_union(
        [Region.from_polytope(p) for p in polys]
    )

    starting_count = len(polys)
    start = time.perf_counter()
    result = Region._simplify_union_of_polys(polys)
    end = time.perf_counter()
    print(f"Computing union of {N ** 5} 5-dimensional polytopes took {end - start} sec.")
    print(f"Starting boxes: {starting_count}, ending boxes: {len(result)}")

    starting_count = len(region.child)
    start = time.perf_counter()
    region.simplify()
    end = time.perf_counter()
    print(f"Computing union of {N ** 5} 5-dimensional polytopes took {end - start} sec.")
    print(f"Starting boxes: {starting_count}, ending boxes: {len(region.child)}")



benchmark_empty(10, 5)
benchmark_union(3)