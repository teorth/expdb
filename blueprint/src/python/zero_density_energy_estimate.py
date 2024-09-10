# code for bounding the function A*(\sigma) in the bound
#
# N*(\sigma - \delta, T) \ll T^{A*(\sigma)(1 - \sigma) + o(1)}
#
# for all \delta > 0, where N(\sigma, T) is the additive energy of 
# the imaginary parts of the zeroes \rho of the Riemann zeta-function 
# satisfying with real part \in [\sigma, 1] and imaginary part \in 
# [-T, T]. 

import additive_energy as ad
from fractions import Fraction as frac
from functions import Interval, RationalFunction as RF
from hypotheses import Hypothesis, Hypothesis_Set
import numpy as np
from reference import Reference

# class representing a zero-density energy estimate 
# For now, this is identical to the Zero_Density_Estimate class 
# TODO: set up an inheritance structure for these classes
class Zero_Density_Energy_Estimate:

    # parameters:
    #   - expr: an expression that represents a function of x
    #   - interval: the Interval object representing the range of validity of the
    #               bound.
    def __init__(self, expr, interval):
        if not isinstance(expr, str):
            raise ValueError("Parameter expr must be of type string")
        if not isinstance(interval, Interval):
            raise ValueError("Parameter interval must be of type Interval")
        self.expr = expr
        self.interval = interval

        # Do not compute this yet
        self.bound = None

    def __repr__(self):
        return f"A(x) \leq {self.expr} on {self.interval}"

    def _ensure_bound_is_computed(self):
        if self.bound is None:
            self.bound = RF.parse(self.expr)

    # Computes the zero-density estimate at a particular value of sigma
    def at(self, sigma):
        if not self.interval.contains(sigma):
            raise ValueError(f'Sigma value {sigma} is out of bounds {self.interval}')

        self._ensure_bound_is_computed()
        return self.bound.at(sigma)

    # -------------------------------------------------------------------------
    # Static methods

    def from_rational_func(rf, interval):
        if not isinstance(rf, RF):
            raise ValueError("Parameter rf must be of type RationalFunction")

        zde = Zero_Density_Energy_Estimate(str(rf), interval)
        zde.bound = rf
        return zde

# Returns a Hypothesis object representing a zero density energy theorem from 
# the literature
def literature_zero_density_energy_estimate(estimate, interval, ref, params=""):
    return Hypothesis(
            f"{ref.author()} ({ref.year()}) zero density estimate" + params,
            "Zero density energy estimate",
            Zero_Density_Energy_Estimate(estimate, interval),
            f"See [{ref.author()}, {ref.year()}]",
            ref,
        )

def derived_zero_density_energy_estimate(data, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived zero density energy estimate",
        "Zero density energy estimate",
        data,
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound

def trivial_zero_density_energy_estimate(data):
    bound = Hypothesis(
        "Trivial zero density energy estimate",
        "Zero density energy estimate",
        data,
        "Trivial",
        Reference.trivial(),
    )
    return bound

# Add trivial energy estimates to a set of hypotheses (taken from Lemma 12.2) 
def add_trivial_zero_density_energy_estimates(hypotheses):
    # add zero density energy estimates arising from A(sigma)
    zdes = hypotheses.list_hypotheses(hypothesis_type="Zero density estimate")
    for zde in zdes:
        hypotheses.add_hypothesis(
            trivial_zero_density_energy_estimate(
                Zero_Density_Energy_Estimate(f"3 * ({zde.data.expr})", zde.data.interval)
            )
        )
    
    # add zero density energy estimates arising from A*(1/2) = 6 and non-increasingness
    # of (1 - \sigma) A*(\sigma), which implies A*(\sigma) \leq 3/(1-\sigma) for \sigma
    # \geq 1/2
    hypotheses.add_hypothesis(
        trivial_zero_density_energy_estimate(
            Zero_Density_Energy_Estimate("3 / (1 - x)", Interval(frac(1,2), 1))
        )
    )

# Given
# - a (sigma, tau, rho*) Region representing feasible LV*(\sigma, \tau) values 
# - a (sigma, tau, rho*) Region representing feasible LV*_{\zeta}(\sigma, \tau) values
# - a range of values of \sigma
# compute the best bound on A*(\sigma) using t0 = 2 and the bound 
# A*(s)(1 - s) \leq 
# max(sup_{2 \leq t < t0} LV*_{\zeta}(s, t)/t, sup_{t0 \leq t \leq 2t0} LV*(s, t)/t)
def approx_best_energy_bound(LV_region, sigma):
    LVs = LV_region.substitute({0: sigma})

    # Take tau0 = 2
    tau0 = 2
    N = 100
    taus = np.linspace(tau0, 3/2 * tau0, N)
    rhos = np.linspace(0, 10, N)
    LV_on_tau = []
    for tau in taus:
        LV_on_tau.append(max(rho for rho in rhos if LVs.contains([tau, rho])) / tau)
    return max(LV_on_tau)

