# This file contains calculations related to prime gaps (Ch. 13)

from fractions import Fraction as frac
from functions import Interval, RationalFunction as RF
from hypotheses import *
import numpy as np
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

    # Set of endpoints 
    crits = {frac(1,2), frac(1)}
    crits.update(zh.data.interval.x0 for zh in best_zd)
    crits.update(zh.data.interval.x1 for zh in best_zd)
    crits.update(zh.data.interval.x0 for zh in best_ze)
    crits.update(zh.data.interval.x1 for zh in best_ze)
    crits = list(crits)
    crits.sort()

    x = RF.x

    for i in range(1, len(crits)):
        interval = Interval(crits[i - 1], crits[i])
        mid = interval.midpoint()

        zdh = next(h for h in best_zd if h.data.interval.contains(mid))
        zeh = next(h for h in best_ze if h.data.interval.contains(mid))
        
        # For simplicity - just work with sympy objects for now (in the future 
        # these will be replaced with RationalFunction instances)
        A = zdh.data.bound.num / zdh.data.bound.den
        B = zeh.data.bound.num / zeh.data.bound.den

        alpha = 4 * x - 2 + 2 * (B * (1 - x) - 1) / (B - A)
        beta = 4 * x - 2 + (B * (1 - x) - 1) / A

        sup = 0
        for sigma in np.linspace(float(interval.x0), float(interval.x1), 100):
            v = max(alpha.subs(x, sigma), beta.subs(x, sigma))
            if sup < v: sup = v
        print(interval, v, alpha.simplify(), beta.simplify())





