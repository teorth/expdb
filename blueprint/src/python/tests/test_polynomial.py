
from polynomial import Polynomial as P

def test_degree():
    assert P(()).degree() == float("-inf")
    assert P((1,)).degree() == 0
    assert P((-1, 1)).degree() == 1
    assert P((1, 0, 2, 3)).degree() == 3

def test_repr():
    assert str(P((1,))) == "1"
    assert str(P((0, 1))) == "x"
    assert str(P((1, -1, -1))) == "1 - x - x^2"
    assert str(P((0, 0, 1, 2))) == "x^2 + 2x^3"

def test_addition():
    assert P((1, 1)) + P((0, 1)) == P((1, 2))
    assert P((1, 1, 1)) + P((1, 1, 1, 1)) == P((2, 2, 2, 1))
    assert P((1, -1, -1)) + P((-1, 1, 1)) == P(())
    assert P((0, 0, -3, 5)) + P((-1, 2, 3, -5)) == P((-1, 2))
    assert P((1, 2, 3, 4)) + P(()) == P((1, 2, 3, 4))

def test_subtract():
    assert P((1, 2, 3)) - P((1, 2)) == P((0, 0, 3))
    assert P((1,)) - P((1, 2, 3, 4)) == P((0, -2, -3, -4))

def test_multiply():
    assert P((1, 1)) * P((1, 1)) == P((1, 2, 1))
    assert P(()) * P((0, -1, 2, 3)) == P(())
    assert P((1,)) * P((1, 2, 3, 4)) == P((1, 2, 3, 4))
    assert P((-1, 4, 5)) * P((3, 0, -1)) == P((-3, 12, 16, -4, -5))

test_degree()
test_repr()
test_addition()
test_subtract()
test_multiply()

print("All test cases passed")
