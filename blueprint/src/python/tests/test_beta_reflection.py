from fractions import Fraction as F
from bound_beta import apply_reflection_beta, classical_bound_beta
from functions import Affine, Interval


def test_reflection_preserves_domain_membership():
    for lower in (False, True):
        for upper in (False, True):
            source = classical_bound_beta(Affine(F(2, 3), F(1, 7),
                Interval(F(1, 4), F(1, 2), lower, upper)))
            result = apply_reflection_beta([source])[0]
            reflected = result.data.bound
            assert reflected.domain == Interval(F(1, 2), F(3, 4), upper, lower)
            assert result.dependencies == {source}
            for x in (F(1, 2), F(5, 8), F(3, 4)):
                assert reflected.domain.contains(x) == source.data.bound.domain.contains(1 - x)
                if reflected.domain.contains(x):
                    assert reflected.at(x) == x - F(1, 2) + source.data.bound.at(1 - x)
            twice = apply_reflection_beta([result])[0].data.bound
            assert twice == source.data.bound


test_reflection_preserves_domain_membership()
