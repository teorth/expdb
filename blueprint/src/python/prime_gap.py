# This file contains calculations related to prime gaps (Ch. 13)

from hypotheses import *
import zero_density_estimate as zd
import zero_density_energy_estimate as ze

# Compute the best estimate of \theta_{gap, 2} using Corollary 13.9
def compute_gap2(hypotheses, debug=False):
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")
    
    best_zd = zd.best_zero_density_estimate(hypotheses, verbose=False)
    for zdh in best_zd:
        print(zdh.data)
    best_ze = ze.compute_best_energy_bound(hypotheses)
    for zeh in best_ze:
        print(zeh.data)

    





