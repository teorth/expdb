import copy
from fractions import Fraction as F

import bound_mu as mu
import exponent_pair as ep
from hypotheses import Hypothesis_Set
from reference import Reference


def hypotheses():
    return Hypothesis_Set([
        ep.trivial_exp_pair,
        ep.literature_exp_pair(F(1, 2), F(1, 2), Reference.make("Test", 2024)),
        ep.literature_exp_pair(F(1, 6), F(2, 3), Reference.make("Test", 2024)),
    ])


def test_both_hulls_can_be_computed_in_either_order():
    for mu_first in (False, True):
        h = hypotheses()
        if mu_first:
            assert mu.best_mu_bound(F(1, 2), h).data.mu == F(1, 6)
        vertices = ep.compute_convex_hull(h)
        assert all(isinstance(v.data, ep.Exp_pair) for v in vertices)
        assert mu.best_mu_bound(F(1, 2), h).data.mu == F(1, 6)
        assert ep.compute_convex_hull(h) == vertices


def test_mutation_invalidates_both_hulls():
    for add in (lambda h, x: h.add_hypothesis(x),
                lambda h, x: h.add_hypotheses([x])):
        h = hypotheses()
        ep.compute_convex_hull(h)
        mu.best_mu_bound(F(1, 2), h)
        stronger = ep.literature_exp_pair(F(0), F(1, 2), Reference.make("Test", 2024))
        add(h, stronger)
        assert stronger in ep.compute_convex_hull(h)
        assert mu.best_mu_bound(F(1, 2), h).data.mu == 0


def test_recomputing_copy_does_not_corrupt_original_cache():
    original = hypotheses()
    assert mu.best_mu_bound(F(1, 2), original).data.mu == F(1, 6)
    derived = copy.copy(original)
    derived.add_hypothesis(mu.classical_bound_mu(F(1, 2), 0))
    assert mu.best_mu_bound(F(1, 2), derived).data.mu == 0
    assert mu.best_mu_bound(F(1, 2), original).data.mu == F(1, 6)


test_both_hulls_can_be_computed_in_either_order()
test_mutation_invalidates_both_hulls()
test_recomputing_copy_does_not_corrupt_original_cache()
