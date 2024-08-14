# Test cases for Affine2 implementation
from fractions import Fraction as frac
import large_values as lv
from polytope import *

def test_edge_case():
    fns = [
        [12, -16, 1],
        [frac(5,2), -4, 1],
        [3, -4, 0],
        [0, -4, 2],
        [3, -6, 2]
    ]
    region = Polytope.rect((frac(1,2), frac(1)), (frac(1), frac(3)))
    m = lv.max_of(fns, region)
    print(m.check((1/2, 1), (1,3)))
    
test_edge_case()
