from functions import Interval, RationalFunction2


def test_multi_coefficient_polynomial_input():
    function = RationalFunction2([1, 2, 3], Interval(0, 1))
    assert function._num_coeffs == (1, 2, 3)
    assert function._den_coeffs == (1,)
    assert str(function) == "x^2 + 2x + 3 on [0,1)"


def test_rational_pair_input_unchanged():
    function = RationalFunction2(([1, 2], [3, 4]), Interval(0, 1))
    assert function._num_coeffs == (1, 2)
    assert function._den_coeffs == (3, 4)


def test_single_coefficient_input_unchanged():
    function = RationalFunction2([2], Interval(0, 1))
    assert function._num_coeffs == (2,)
    assert function._den_coeffs == (1,)
