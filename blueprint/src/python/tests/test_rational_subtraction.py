from fractions import Fraction

from functions import RationalFunction


def test_subtraction_operator():
    f = RationalFunction([1, 1], [1, 2])
    g = RationalFunction([2, 1], [1, 3])
    for other in (g, 2, Fraction(1, 3)):
        result = f - other
        for x in (Fraction(0), Fraction(1), Fraction(3, 2)):
            expected = f.at(x) - (other.at(x) if isinstance(other, RationalFunction) else other)
            assert result.at(x) == expected
    assert (f - f).at(Fraction(1)) == 0


test_subtraction_operator()
