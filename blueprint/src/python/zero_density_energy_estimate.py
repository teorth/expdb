# code for bounding the function A*(\sigma) in the bound
#
# N*(\sigma - \delta, T) \ll T^{A*(\sigma)(1 - \sigma) + o(1)}
#
# for all \delta > 0, where N(\sigma, T) is the additive energy of 
# the imaginary parts of the zeroes \rho of the Riemann zeta-function 
# satisfying with real part \in [\sigma, 1] and imaginary part \in 
# [-T, T]. 

from fractions import Fraction as frac
from functions import Interval, RationalFunction as RF, SympyHelper
from hypotheses import Hypothesis, Hypothesis_Set
import numpy as np
from reference import Reference
from region import Region, Region_Type
import sympy

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
        return f"A*(x) \leq {self.expr} on {self.interval}"

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

# Numerically approximate 
# sup_{tau_lower \leq t <= tau_upper} LV*(s, t)/t for a fixed s, given a 2-dimensional
# region representing the feasible values of (t, LV*(s, t))
def approx_sup_LV_on_tau(LVER, sigma, tau_lower, tau_upper):
    taus = np.linspace(tau_lower, tau_upper, 500)
    rhos = np.linspace(0, 10, 100)
    rho1s = np.linspace(0, 10, 100)
    ss = np.linspace(0, 10, 100)

    _max = 0
    for tau in taus:
        for rho1 in rho1s:
            ratio = rho1 / tau
            if ratio < _max:
                continue
            if any(LVER.contains([sigma, tau, rho, rho1, s]) 
                   for rho in rhos 
                   for s in ss):
                _max = ratio
                        
    return _max

# Given
# - a (sigma, tau, rho*) Region representing feasible LV*(\sigma, \tau) values 
# - a (sigma, tau, rho*) Region representing feasible LV*_{\zeta}(\sigma, \tau) values
# - a range of values of \sigma
# compute the best bound on A*(\sigma) using t0 = 2 and the bound 
# A*(s)(1 - s) \leq 
# max(sup_{2 \leq t < t0} LV*_{\zeta}(s, t)/t, sup_{t0 \leq t \leq 2t0} LV*(s, t)/t)
def approx_best_energy_bound(hypotheses, sigma, tau0):

    # 1. A large value energy region is also a zeta large value energy region
    # 2. Large value energy regions have the raise to power hypothesis, which is 
    # a type of Large value energy region transform
    # 3. Zeta large value energy regions do not have the raise to power hypothesis
    # Therefore, the set of zeta large value energy regions can be obtained by 
    # first expanding the set of large value energy regions (by raising to a power)
    # then adding in the additional zeta large value energy regions later. 
    lvers = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")
        
    # Compute the large value energy bounding region
    # Use LVER transformations and use them to expand the set of LVERs
    transforms = hypotheses.list_hypotheses(hypothesis_type="Large value energy region transform")
    transformed_lvers = []
    for tf in transforms:
        transformed_lvers.extend(tf.data.transform(lver) for lver in lvers)
    lvers.extend(transformed_lvers)

    LVER = Region(Region_Type.INTERSECT, [lver.data.region for lver in lvers])
    sup1 = approx_sup_LV_on_tau(LVER, sigma, tau0, 2 * tau0)
    
    # Compute the zeta large value energy bounding region
    lvers.extend(hypotheses.list_hypotheses(hypothesis_type="Zeta large value energy region"))

    LVER_zeta = Region(Region_Type.INTERSECT, [lver.data.region for lver in lvers])
    sup2 = approx_sup_LV_on_tau(LVER_zeta, sigma, 2, tau0)

    return max(sup1, sup2)

# Given two Hypothesis objects of type "Large value energy region", representing 
# an estimate of the large value energy region and zeta large value energy region
# respectively, computes and returns the best bound on A*(\sigma) 
# 
# This function should give the same result as approx_best_energy_bound
def lver_to_energy_bound(LVER, LVER_zeta, sigma_interval, debug=False):
    
    if LVER is not None and (not isinstance(LVER, Hypothesis) or \
        LVER.hypothesis_type != "Large value energy region"):
        raise ValueError("Parameter LVER must be a Hypothesis of type 'Large value energy region'.")
    if LVER_zeta is not None and (not isinstance(LVER_zeta, Hypothesis) or \
        LVER_zeta.hypothesis_type != "Zeta large value energy region"):
        raise ValueError("Parameter LVER_zeta must be a Hypothesis of type 'Zeta large value energy region'.")
    
    fns = []
    deps = set()
    depcount = [0, 0] # The number of dependencies

    if LVER is not None:
        sup1 = compute_sup_LV_on_tau(LVER.data.region, sigma_interval)
        if debug:
            print("sup1")
            for s in sup1: print(s[0], "for x\in", s[1])
        fns.extend(list(sup1))
        deps.update(list(LVER.dependencies))
        depcount[0] = len(LVER.dependencies)
    
    if LVER_zeta is not None:
        sup2 = compute_sup_LV_on_tau(LVER_zeta.data.region, sigma_interval)
        if debug:
            print("sup2")
            for s in sup2: print(s[0], "for x\in", s[1])
        fns.extend(list(sup2))
        deps.update(list(LVER_zeta.dependencies))
        depcount[1] = len(LVER_zeta.dependencies)
    
    # Take maximum
    bounds, indexes = [(f[0], f[1]) for f in fns], [(f[2],) for f in fns]
    sup = RF.max(bounds, sigma_interval)
    print("A*(x)(1-x) \leq")
    for s in sup: print(s[0], "for x \in", s[1])
    
    return [
        derived_zero_density_energy_estimate(
            Zero_Density_Energy_Estimate.from_rational_func(s[0], s[1]),
            f"Follows from combining {depcount[0]} large value energy regions " + \
            "and {depcount[1]} zeta large value energy regions",
            deps
        )
        for s in sup
        ]

def compute_sup_LV_on_tau(LV_region, sigma_interval):
    
    # assume that LV_region is a union of polytopes and that the (sigma-tau) domain
    # is already correct
    polys = [r.child for r in LV_region.child]

    # each polytope is 3-dimensional. Find all edges and project them onto the 
    # sigma dimension. For those with a non-zero projection, compute rho / tau along 
    # the edge as a function of sigma 
    fns = []
    visited = set() # Keep track of the edges we have already visited to prevent duplication
    for p in polys:
        edges = p.get_edges()
        for edge in edges:
            # Vertices on either side of the edge
            v1 = tuple(edge[0])
            v2 = tuple(edge[1])
            if v1 + v2 in visited or v2 + v1 in visited: 
                continue
            visited.add(v1 + v2)

            (sigma1, tau1, rho1) = v1
            (sigma2, tau2, rho2) = v2

            # Skip edges with empty projection onto the sigma domain
            if sigma1 == sigma2: continue

            # tau as a (linear) function of sigma along this edge
            tau = RF([
                (tau1 - tau2) / (sigma1 - sigma2), 
                (sigma1 * tau2 - sigma2 * tau1) / (sigma1 - sigma2)
            ])
            # rho as a linear function of sigma along this edge
            rho = RF([
                (rho1 - rho2) / (sigma1 - sigma2), 
                (sigma1 * rho2 - sigma2 * rho1) / (sigma1 - sigma2)
            ])

            fns.append((rho.div(tau), Interval(min(sigma1, sigma2), max(sigma1, sigma2))))
    
    # Take the maximum of the functions
    return RF.max(fns, sigma_interval)

# Given a Hypothesis_Set, compute a piecewise function representing the best bound 
# on A*(\sigma)
# TODO: this function is quite similar to compute_best_density_estimate in zero_density_estimate.py
# Combine/merge them?
def compute_best_energy_bound(hypotheses: Hypothesis_Set) -> list:
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")
    
    hs = hypotheses.list_hypotheses(hypothesis_type="Zero density energy estimate")

    # Ensure bound is computed (i.e any RationalFunction objects that use lazy initialization
    # is parsed)
    for h in hs:
        h.data._ensure_bound_is_computed()

    minimum = RF.min([(h.data.bound, h.data.interval) for h in hs], Interval(frac(1,2), 1))

    # Pack into Hypothesis objects
    best_bound = []
    for (func, interval, ref) in minimum:
        zde = Zero_Density_Energy_Estimate.from_rational_func(func, interval)
        if ref < 0:
            best_bound.append(
                Hypothesis(
                    "Placeholder zero density estimate",
                    "Zero density estimate",
                    zde,
                    "Placeholder zero density estimate",
                    Reference.trivial(),
                )
            )
        else:
            h = hs[ref]
            best_bound.append(
                Hypothesis(h.name, h.hypothesis_type, zde, h.proof, h.reference)
            )
    return best_bound
