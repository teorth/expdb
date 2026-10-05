from functions import Interval, RationalFunction2


def function(numerator, denominator=(1,)):
    return RationalFunction2((list(numerator), list(denominator)), Interval(0, 1))


def test_monomials_have_their_actual_degree():
    linear = function((1, 0))
    quadratic = function((1, 0, 0))
    assert not linear.is_constant()
    assert linear.is_affine()
    assert not quadratic.is_constant()
    assert not quadratic.is_affine()


def test_leading_zeros_do_not_raise_degree():
    constant = function((0, 0, 2), (0, 1))
    assert constant.is_constant()
    assert constant.is_affine()
    assert str(constant) == "2.0 on [0,1)"
    assert function((0, 1, 2)).is_affine()


def test_nonconstant_denominator_is_not_affine():
    reciprocal = function((1,), (1, 0))
    assert not reciprocal.is_constant()
    assert not reciprocal.is_affine()
