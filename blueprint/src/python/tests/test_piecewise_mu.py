from fractions import Fraction as F
from bound_mu import best_mu_bound, best_mu_bound_piecewise, classical_bound_mu
from functions import Affine, Interval
from hypotheses import Hypothesis_Set


def test_piecewise_mu_domains_and_pointwise_values():
    for stronger in ([], [classical_bound_mu(F(1, 2), F(1, 6))]):
        for left, right in [(-2, -1), (2, 3), (-1, 2), (0, 1), (F(1, 4), F(3, 4)), (F(1, 2), F(1, 2))]:
            for lower, upper in [(True, True), (False, False), (True, False), (False, True)]:
                domain = Interval(left, right, lower, upper)
                h = Hypothesis_Set(stronger)
                pieces = best_mu_bound_piecewise(domain, h)
                assert all(isinstance(p, Affine) for p in pieces)
                for x in [left - 1, left, (F(left) + F(right)) / 2, right, right + 1]:
                    containing = [p for p in pieces if p.domain.contains(x)]
                    assert bool(containing) == domain.contains(x)
                    if containing:
                        expected = best_mu_bound(x, Hypothesis_Set(stronger)).data.mu
                        assert all(p.at(x) == expected for p in containing)


test_piecewise_mu_domains_and_pointwise_values()
