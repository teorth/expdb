from fractions import Fraction as frac

from constants import Constants
from functions import Affine, Interval
from polytope import Polytope
from region import Region
from hypotheses import *
from literature import *
import exponent_pair as ep
import zero_density_estimate as zd


def run_ep_to_zd_regression_tests():
    """Regression tests for three defects in zero_density_estimate.py.

    BUG 1  ep_to_zd passed a list to bourgain_ep_to_zd's `hypotheses=` slot
           (-> AttributeError: 'list' object has no attribute 'add_hypotheses'),
           and did not seed literature beta bounds (-> empty hull).
    BUG 2  ivic_ep_to_zd built "6/(4x + 0)" for m=2 (no '*', dangling '+ 0')
           (-> sympy SyntaxError).
    BUG 3  ivic_ep_to_zd dereferenced dep when None (latent until BUG 1 fixed).
    """
    lit = literature

    # --- BUG 2: ivic_ep_to_zd parses and matches 3m/((3m-2)s + (2-m)) ---
    hs = Hypothesis_Set()
    for ht in (
        "Exponent pair",
        "Exponent pair bound",
        "Beta bound",
        "Upper bound on beta",
    ):
        hs.add_hypotheses(lit.list_hypotheses(hypothesis_type=ht))
    hs.add_hypotheses(ep.compute_exp_pairs(hs, search_depth=3, prune=True))
    ephs = hs.list_hypotheses(hypothesis_type="Exponent pair")
    for m in (2, 3, 4):
        z = zd.ivic_ep_to_zd(ephs, m=m)
        z.data._ensure_bound_is_computed()  # raised SyntaxError pre-patch (m=2)
        s = frac(9, 10)
        expected = frac(3 * m) / ((3 * m - 2) * s + (2 - m))
        assert abs(float(z.data.at(s)) - float(expected)) < 1e-12

    # --- BUG 1 + 3: ep_to_zd runs and returns parseable estimates ---
    hs2 = Hypothesis_Set()
    hs2.add_hypotheses(lit.list_hypotheses(hypothesis_type="Zero density estimate"))
    out = zd.ep_to_zd(hs2)  # raised AttributeError pre-patch
    assert len(out) > 0
    for h in out:
        h.data._ensure_bound_is_computed()


def run_density_normalization_regression_tests():
    import additive_energy as ad
    import derived
    import examples
    import large_values as lv
    from unittest.mock import patch

    def check_values(bounds, expected, samples):
        assert bounds
        for sigma in samples:
            values = [h.data.at(sigma) for h in bounds if h.data.interval.contains(sigma)]
            assert values == [expected(sigma)], (sigma, values, expected(sigma))
        assert not any(h.data.interval.contains(1) for h in bounds)

    # The named derivation must use the same A(sigma) convention as the
    # literature entry, not return the exponent (1-sigma) A(sigma).
    samples = [frac(1, 2), frac(3, 5), frac(2, 3), frac(7, 10)]
    ingham = literature.find_hypothesis(name="Ingham (1940) zero density estimate")
    bounds = derived.prove_zero_density_ingham_1940(verbose=False)
    check_values(bounds, ingham.data.at, samples)
    assert not any(h.data.interval.contains(frac(3, 4)) for h in bounds)
    best = zd.best_zero_density_estimate(Hypothesis_Set(bounds))
    check_values(best, ingham.data.at, samples)

    # The plotting caller previously compensated for the incorrect convention.
    # Capture its data without creating a figure; normalization must occur once.
    with patch.object(derived.np, "linspace", return_value=samples), \
            patch.object(derived, "plt") as plotting:
        derived.prove_zero_density(
            [], False, Interval(frac(1, 2), frac(3, 4)), "Ingham", plot=True
        )
    plotting.plot.assert_called_once_with(samples, [ingham.data.at(s) for s in samples])

    with patch.object(examples, "plt") as plotting:
        examples.plot([ingham], [ingham])
    assert plotting.plot.call_count == 2
    assert plotting.plot.call_args_list[0].args[:2] == plotting.plot.call_args_list[1].args[:2]

    # L2 alone gives B(sigma) = (7-2*sigma)/6 on the tau range [3,6].
    # Request a closed upper endpoint to ensure sigma=1 is not divided by zero.
    interval = Interval(frac(1, 2), 1, True, True)
    samples = [frac(1, 2), frac(2, 3), frac(3, 4), frac(9, 10)]
    expected = lambda s: (7 - 2 * s) / (6 * (1 - s))
    hypotheses = Hypothesis_Set([lv.large_value_estimate_L2])
    bounds = zd.lv_zlv_to_zd(hypotheses, interval, tau0=frac(3))
    check_values(bounds, expected, samples)

    # Exercise the affine-tau0 converter with the trivial estimate rho <= tau,
    # whose supremum rho/tau is identically 1.
    hypotheses = Hypothesis_Set([lv.classical_LV_estimate([[0, 0, 1]])])
    bounds = zd.lv_zlv_to_zd2(hypotheses, interval, Affine(0, frac(3), interval))
    check_values(bounds, lambda s: 1 / (1 - s), samples)

    # Lift the same L2 estimate to five dimensions, restrict its two tau ranges,
    # and exercise the energy-region converter through its real projection path.
    base = ad.lv_to_lver(Hypothesis_Set([lv.large_value_estimate_L2]))[0]

    def restrict_lver(tau_range, zeta=False):
        domain = Region.from_polytope(
            Polytope.rect((frac(1, 2), frac(1)), tau_range)
        ).lift([0, 1] + [(0, Constants.LV_DEFAULT_UPPER_BOUND)] * 3)
        data = ad.Large_Value_Energy_Region(Region.intersect([base.data.region, domain]))
        constructor = (ad.derived_zeta_large_value_energy_region if zeta else
                       ad.derived_large_value_energy_region)
        return constructor(data, "Restrict L2 to the tau range", {base})

    bounds = zd.lver_to_zd(
        restrict_lver((frac(3), frac(6))),
        restrict_lver((frac(2), frac(3)), zeta=True),
        interval,
    )
    check_values(bounds, expected, samples)

    for empty in (Interval(1, 1, True, True), Interval(frac(2, 3), frac(2, 3))):
        hypotheses = Hypothesis_Set([lv.large_value_estimate_L2])
        assert zd.lv_zlv_to_zd(hypotheses, empty, tau0=frac(3)) == []
        assert zd.lv_zlv_to_zd2(hypotheses, empty, Affine(0, frac(3), empty)) == []
        assert zd.lver_to_zd(None, None, empty) == []


run_density_normalization_regression_tests()
print("All zero-density normalization regression test cases passed.")
run_ep_to_zd_regression_tests()
