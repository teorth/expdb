from fractions import Fraction as frac

from hypotheses import Hypothesis_Set
import zero_density_estimate as zd


def run_density_tau_domain_regression_tests():
    import derived
    import large_values as lv
    from constants import Constants
    from functions import Affine, Interval
    from unittest.mock import patch

    # Compare the real named derivations using their shared output convention;
    # this tests the tau domain independently of zero-density normalization.
    method1 = derived.prove_zero_density_ingham_1940(verbose=False)
    method2 = derived.prove_zero_density_ingham_1940_v2(verbose=False)
    for sigma in (frac(3, 5), frac(2, 3), frac(7, 10)):
        values1 = [h.data.at(sigma) for h in method1 if h.data.interval.contains(sigma)]
        values2 = [h.data.at(sigma) for h in method2 if h.data.interval.contains(sigma)]
        assert len(values1) == 1
        assert values2 == values1, (sigma, values2, values1)

    # Inspect the actual domain passed to the real region computation, not a
    # reconstructed inequality. A nonconstant tau0 checks both affine terms.
    interval = Interval(frac(1, 2), frac(3, 4))
    tau0 = Affine(2, 1, interval)
    hypotheses = Hypothesis_Set([lv.large_value_estimate_L2])
    with patch.object(zd, "compute_large_value_region", wraps=zd.compute_large_value_region) as compute:
        bounds = zd.lv_zlv_to_zd2(hypotheses, interval, tau0)
    assert bounds
    assert compute.call_count == 2
    domain = compute.call_args_list[1].args[1]
    assert compute.call_args_list[1].kwargs["zeta"] is True
    for sigma in (interval.x0, frac(2, 3), interval.x1):
        upper = frac(4, 3) * (tau0.m * sigma + tau0.c)
        assert domain.contains([sigma, 2, 0])
        assert domain.contains([sigma, upper, Constants.LV_DEFAULT_UPPER_BOUND])
        assert domain.contains([sigma, (2 + upper) / 2, 1])
        assert not domain.contains([sigma, upper + frac(1, 100), 0])
        assert not domain.contains([sigma, 2 - frac(1, 100), 0])
    for point in (
        [interval.x0 - frac(1, 100), 2, 0],
        [interval.x1 + frac(1, 100), 2, 0],
        [frac(2, 3), 2, -1],
        [frac(2, 3), 2, Constants.LV_DEFAULT_UPPER_BOUND + 1],
    ):
        assert not domain.contains(point), point


run_density_tau_domain_regression_tests()
print("All zero-density tau-domain regression test cases passed.")
