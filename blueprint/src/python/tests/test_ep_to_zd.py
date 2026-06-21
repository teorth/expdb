from fractions import Fraction as frac

from ..hypotheses import *
from ..literature import *
from .. import exponent_pair as ep
from .. import zero_density_estimate as zd


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


run_ep_to_zd_regression_tests()
