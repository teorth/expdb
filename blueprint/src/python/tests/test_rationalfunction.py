# Test cases for RationalFunction implementation
from fractions import Fraction as frac
from functions import RationalFunction2 as RF, Interval

def test_repr():
    domain = Interval(0, 1)
    I = frac(1)
    assert str(RF([0], domain)) == "0 on [0,1)"
    assert str(RF(([0], [1, 1, 1]), domain)) == "0 on [0,1)"
    assert str(RF(([I], [I]), domain)) == "1 on [0,1)"
    assert str(RF(([I], [frac(2)]), domain)) == "1/2 on [0,1)"
    assert str(RF(([I, I], [I]), domain)) == "x + 1 on [0,1)"
    assert str(RF(([I, I, I], [I, I]), domain)) == "(x^2 + x + 1)/(x + 1) on [0,1)"

test_repr()
