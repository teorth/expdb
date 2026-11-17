# Test cases for Affine2 implementation
from fractions import Fraction as frac
from functions import Affine2, Piecewise
from polytope import *


# The pointwise minimum of a list of affine functions over a common domain,
# as a Piecewise function.
def min_of(fns, region):
    f = Piecewise([Affine2(fns[0], region, 0)])
    for i in range(1, len(fns)):
        f = f.min_with(Piecewise([Affine2(fns[i], region, i)]))
    return f


def test_edge_case():
    fns = [
        [12, -16, 1],
        [frac(5,2), -4, 1],
        [3, -4, 0],
        [0, -4, 2],
        [3, -6, 2]
    ]
    region = Polytope.rect((frac(1,2), frac(1)), (frac(1), frac(3)))
    m = min_of(fns, region)
    # The pieces tile the region without overlapping
    assert m.check((frac(1,2), frac(1)), (frac(1), frac(3)))
    # and each piece really carries the minimum of the five functions
    for i in range(11):
        for j in range(11):
            x = frac(1,2) + frac(i, 22)
            y = frac(1) + frac(2 * j, 10)
            assert m.at([x, y]) == min(f[0] + f[1] * x + f[2] * y for f in fns)

test_edge_case()


def test_equal_affine_functions_with_optional_labels():
    domain = Polytope.rect((0, 1), (0, 1))
    for left, right, expected in [(None, None, None), (None, 2, 2), (2, None, 2), (2, 1, 1)]:
        a = Affine2([1, 2, 3], domain, left)
        b = Affine2([1, 2, 3], domain, right)
        pieces = a.min_with([b])
        assert len(pieces) == 1
        assert pieces[0].label == expected
        assert pieces[0].at([frac(1, 2), frac(1, 2)]) == frac(7, 2)


test_equal_affine_functions_with_optional_labels()
