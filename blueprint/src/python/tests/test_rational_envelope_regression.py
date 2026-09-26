from functions import RationalFunction as R, Interval
from fractions import Fraction as F
import sympy

def run_regressions():
    domain = Interval(-2, 2, True, True)
    for method, value in ((R.min, 2000000), (R.max, -2000000)):
        result = method([(R([value]), domain)], domain)
        assert len(result) == 1 and result[0][2] == 0
        assert result[0][1] == domain and result[0][0].at(0) == value
    for method, choose, default in ((R.min, min, sympy.oo), (R.max, max, -sympy.oo)):
        inputs = [(R([1], [1, 0]), Interval(-2, 2, False, True)),
                  (R([0]), Interval(-1, 1, True, False)),
                  (R([-3]), Interval(2, 2, True, True))]
        pieces = method(inputs, domain)
        for x in [F(k, 4) for k in range(-8, 9)]:
            active = [(f, i) for f, cell, i in pieces if cell.contains(x)]
            assert len(active) == 1
            values = [f.at(x) for f, cell in inputs
                      if cell.contains(x) and sympy.sympify(f.den).subs(R.x, x) != 0]
            expected = choose(values) if values else default
            assert active[0][0].at(x) == expected
            assert (active[0][1] == -1) == (not values)
        assert all(len(item) == 2 for item in method(inputs, domain, False))
    singleton = Interval(1, 1, True, True)
    assert R.min([(R([7]), singleton)], singleton)[0][1] == singleton
    assert R.min([], Interval(1, 1)) == []

run_regressions()
