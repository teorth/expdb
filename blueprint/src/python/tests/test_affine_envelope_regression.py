from functions import Affine, Interval
from fractions import Fraction as F

def run_regressions():
    f = Affine(1, 0, Interval(0, 4, True, True))
    g = Affine(-1, 4, Interval(1, 3, False, True))
    originals = [v.domain.deep_copy() for v in (f, g)]
    for method, choose in ((f.min_with, min), (f.max_with, max)):
        domain = Interval(F(1, 2), F(7, 2), False, True)
        pieces = method([g], domain)
        for x in [F(k, 4) for k in range(17)]:
            active = [p for p in pieces if p.domain.contains(x)]
            if not domain.contains(x):
                assert not active
            else:
                assert len(active) == 1
                values = [v.at(x) for v in (f, g) if v.domain.contains(x)]
                assert active[0].at(x) == choose(values)
        assert f.domain == originals[0] and g.domain == originals[1]
        assert all(p is not f and p is not g for p in pieces)
    only = f.min_with([], Interval(1, 2, False, False))
    assert len(only) == 1 and only[0].domain == Interval(1, 2, False, False)
    assert f.domain == originals[0]

run_regressions()
