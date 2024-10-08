# This file contains calculations related to prime gaps (Ch. 13)

from hypotheses import *
import zero_density_estimate as zd
import zero_density_energy_estimate as ze

# Compute the best estimate of \theta_{gap, 2} using Corollary 13.9
def compute_gap2(hypotheses, debug=False):
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")
    
    # Compute the best zero-density estimate
    zds = hypotheses.list_hypotheses(hypothesis_type="Zero density estimate")
    if debug:
        n_derived = sum(1 for z in zds if z.reference.is_derived())
        print(f"Collected {len(zds)} zero-density estimates (including {n_derived} derived estimates)")

    zd.best_zero_density_estimate(hypotheses, verbose=True)

    # Compute the best zero-density energy estimate
    zes = hypotheses.list_hypotheses(hypothesis_type="Zero density energy estimate")
    if debug:
        n_derived = sum(1 for z in zes if z.reference.is_derived())
        print(f"Collected {len(zes)} additive energy estimates (including {n_derived} derived estimates)")