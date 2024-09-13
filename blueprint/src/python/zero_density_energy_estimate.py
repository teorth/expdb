# code for bounding the function A*(\sigma) in the bound
#
# N*(\sigma - \delta, T) \ll T^{A*(\sigma)(1 - \sigma) + o(1)}
#
# for all \delta > 0, where N(\sigma, T) is the additive energy of 
# the imaginary parts of the zeroes \rho of the Riemann zeta-function 
# satisfying with real part \in [\sigma, 1] and imaginary part \in 
# [-T, T]. 

import additive_energy as ad
import cdd
from constants import Constants
from fractions import Fraction as frac
from functions import Interval, RationalFunction as RF, SympyHelper
from hypotheses import Hypothesis, Hypothesis_Set
import numpy as np
from polytope import Polytope
from reference import Reference
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

# Numerically approximate 
# sup_{tau_lower \leq t <= tau_upper} LV*(s, t)/t for a fixed s, given a 2-dimensional
# region representing the feasible values of (t, LV*(s, t))
def approx_sup_LV_on_tau(tau_rho_region, tau_lower, tau_upper):
    taus = np.linspace(tau_lower, tau_upper, 500)

    # TODO: Project to find the largest feasible value of rho
    rhos = np.linspace(0, 10, 100)

    LV_on_tau = []
    for tau in taus:
        LV_on_tau.append(max(rho for rho in rhos if tau_rho_region.contains([tau, rho])) / tau)
    return max(LV_on_tau)

# Given
# - a (sigma, tau, rho*) Region representing feasible LV*(\sigma, \tau) values 
# - a (sigma, tau, rho*) Region representing feasible LV*_{\zeta}(\sigma, \tau) values
# - a range of values of \sigma
# compute the best bound on A*(\sigma) using t0 = 2 and the bound 
# A*(s)(1 - s) \leq 
# max(sup_{2 \leq t < t0} LV*_{\zeta}(s, t)/t, sup_{t0 \leq t \leq 2t0} LV*(s, t)/t)
def approx_best_energy_bound(LV_region, LVZ_region, sigma, tau0):

    LVs = LV_region.substitute({0: sigma})
    LVs.plot2d((tau0, 2 * tau0), (0, 20), resolution=500)
    #proj = LVs.project({1})
    #proj.simplify()
    #print("rho projection", proj)

    sup1 = approx_sup_LV_on_tau(LVs, tau0, 2 * tau0)
    
    LVZs = LVZ_region.substitute({0: sigma})
    sup2 = approx_sup_LV_on_tau(LVZs, 2, tau0)

    #print(sup1, sup2)

    return max(sup1, sup2)

# Given two Hypothesis objects of type "Large value energy region", representing 
# an estimate of the large value energy region and zeta large value energy region
# respectively, computes and returns the best bound on A*(\sigma) 
# 
# This function should give the same result as approx_best_energy_bound
def compute_best_energy_bound(LVER, LVER_zeta, sigma_interval, tau0):
    
    if not isinstance(LVER, Hypothesis) or \
        LVER.hypothesis_type != "Large value energy region":
        raise ValueError("Parameter LVER must be a Hypothesis of type 'Large value energy region'.")
    if not isinstance(LVER_zeta, Hypothesis) or \
        LVER_zeta.hypothesis_type != "Zeta large value energy region":
        raise ValueError("Parameter LVER_zeta must be a Hypothesis of type 'Zeta large value energy region'.")
        
    sup1 = compute_sup_LV_on_tau(LVER.data.region, sigma_interval, (tau0, 2 * tau0))
    print("sup1")
    for s in sup1: print(s[0], "for x\in", s[1])
    
    sup2 = compute_sup_LV_on_tau(LVER_zeta.data.region, sigma_interval, (2, tau0))
    print("sup2")
    for s in sup2: print(s[0], "for x\in", s[1])
    
    # Take maximum
    sup = max_RF(list(sup1) + list(sup2), Interval(sigma_interval[0], sigma_interval[1]))
    print("A*(x)(1-x) \leq")
    for s in sup: print(s[0], "for x\in", s[1])
    
    return [
        derived_zero_density_energy_estimate(
            Zero_Density_Energy_Estimate.from_rational_func(s[0], s[1]),
            f"Follows from combining {len(LVER.dependencies)} large value energy regions " + \
            "and {len(LVER_zeta.dependencies)} zeta large value energy regions",
            set(list(LVER.dependencies) + list(LVER_zeta.dependencies))
        )
        for s in sup
        ]



def compute_sup_LV_on_tau(LV_region, sigma_interval, tau_interval):
    
    # assume that LV_region is a union of polytopes 
    polys = [r.child for r in LV_region.child]

    # Crop the appropriate domain
    domain = Polytope.rect(
        sigma_interval, 
        tau_interval, 
        (0, Constants.LV_DEFAULT_UPPER_BOUND)
    )
    cropped_polys = [p.intersect(domain) for p in polys]
    polys = [p for p in cropped_polys if not p.is_empty(include_boundary=False)]

    # each polytope is 3-dimensional. Find all edges and project them onto the 
    # sigma dimension. For those with a non-zero projection, compute rho / tau along 
    # the edge as a function of sigma 
    fns = []  
    for p in polys:
        edges = p.get_edges()
        for edge in edges:
            # Vertices on either side of the edge
            (sigma1, tau1, rho1) = edge[0]
            (sigma2, tau2, rho2) = edge[1]

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
    return max_RF(fns, Interval(sigma_interval[0], sigma_interval[1]))

# TODO: this function is very similar to a function in zero_density_estimate.py - 
# merge/combine them?
def max_RF(fns, domain):
    # Then, compute their maximums as a piecewise RationalFunction
    # at this point we can just make use of the same code as best density estimate solver 
    x = RF.x
    sup = [(RF.parse("0"), domain)] # placeholder
    
    for (func1, in1) in fns:
        # For simplicity, work directly with sympy objects
        f1 = func1.num / func1.den
        new_sup = []
        for (func2, in2) in sup:
            f2 = func2.num / func2.den
            solns = sympy.solve(f1 - f2)
            crits = set(SympyHelper.to_frac(soln) for soln in solns if soln.is_real)
            crits.update([in1.x0, in1.x1])
            crits = set(c for c in crits if in2.contains(c))
            crits.update([in2.x0, in2.x1])

            crits = list(crits)
            crits.sort()
            for i in range(1, len(crits)):
                inter = Interval(crits[i - 1], crits[i])
                mid = inter.midpoint()
                
                if not in1.contains(mid):
                    new_sup.append((func2, inter))
                elif f2.subs(x, mid) >= f1.subs(x, mid):
                    new_sup.append((func2, inter))
                else:
                    new_sup.append((func1, inter))
        
        # Simplify
        sup = []
        i = 0
        while i < len(new_sup):
            (fi, inti) = new_sup[i]
            left = inti.x0
            right = inti.x1

            j = i + 1
            while j < len(new_sup):
                (fj, intj) = new_sup[j]
                if not (fi == fj and right == intj.x0):
                    break
                right = intj.x1
                j += 1
                
            sup.append((fi, Interval(left, right)))
            i = j

    return sup


