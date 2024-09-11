
# Some code for benchmarking frequently used Polytope functions
import parent 

from polytope import Polytope
import random as rd
import time

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


benchmark_empty(100, 5)