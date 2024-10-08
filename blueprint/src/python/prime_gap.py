# This file contains calculations related to prime gaps (Ch. 13)

from hypotheses import *
import zero_density_estimate as zd
import zero_density_energy_estimate as ze

# Compute the best estimate of \theta_{gap, 2} using Corollary 13.9
def compute_gap2(hypotheses):
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")
    
    # Compute the best zero-density estimate
    zds = hypotheses.list_hypotheses(hypothesis_type="Zero density estimate")
    for zdh in zds:
        print(zdh.data)

    # Compute the best zero-density energy estimate
    zes = hypotheses.list_hypotheses(hypothesis_type="Zero density energy estimate")
    for zeh in zes:
        print(zeh.data)