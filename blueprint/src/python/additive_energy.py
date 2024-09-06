# code implementation for large value additive energy region theorems

import copy 
from constants import Constants
import exponent_pair as ep
from fractions import Fraction as frac
from functions import *
from hypotheses import Hypothesis
from polytope import Polytope
from reference import Reference
from region import Region, Region_Type
from transform import Transform

class Large_Value_Energy_Region:
    
    # Construct a valid large value energy region represented as a union of 
    # 5-dimensional polytopes 
    # Parameters: region (Region) an object of type Region which represents the 
    # large value energy region
    def __init__(self, region):
        # Should this be a set instead of a list?
        if not isinstance(region, Region):
            raise ValueError("Parameter region must be of type Region.")
        self.region = region
    
    
    def __repr__(self):
        return str(self.region)
        
    def __copy__(self):
        return Large_Value_Energy_Region(copy.copy(self.region))
    
    # Static methods ---------------------------------------------------------

    # Computes a Large_Value_Energy_Region object representing the union of 
    # a list of polytopes, where each polytope is defined as the intersection of 
    # - a box, represented as a list of lists (a list of constraints)
    # - a halfplane, represented as a single list (a constraint)
    # The resulting Region object is represented as a DISJOINT_UNION, which 
    # greatly improves performance over a UNION representation
    # 
    # This method is useful for quickly initializing many commonly encountered 
    # regions in the study of large value energy regions, since they correspond to 
    # a region implied by a single max() function. 
    def union_of_halfplanes(halfplanes, box):

        # Once a halfplane has been added, include its complement in the list of 
        # neg_ineq, to add as a constraint to all remaining polytopes to be 
        # constructed. 
        neg_ineq = []
        polys = []
        for hp in halfplanes:
            polys.append(Region(Region_Type.POLYTOPE, Polytope(box + [hp] + neg_ineq)))
            neg_ineq.append([-x for x in hp])
        return Region.disjoint_union(polys)


    # The default bounds on the tuple (sigma, tau, rho, rho*, s). This is to ensure
    # that all large value regions are finite regions. 
    def default_constraints():
        limits = [
            (frac(1,2), frac(1)),
            (0, Constants.TAU_UPPER_LIMIT),
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND)
            ]
        
        dim = len(limits)
        bounds = []
        for i in range(dim):
            lim = limits[i]
            
            b = [-lim[0]] + ([0] * dim)
            b[i + 1] = 1
            bounds.append(b)   # x >= lim[0]
            
            b = [lim[1]] + ([0] * dim)
            b[i + 1] = -1
            bounds.append(b)   # x <= lim[1]
        return bounds

    # ------------------------------------------------------------------------
    
    # Returns whether the region contains a 5-dimensional point 
    def contains(self, point):
        if len(point) != 5: 
            raise ValueError("point must be 5-dimensional")
        return self.region.contains(point)
    
    # Raise this region to the kth power
    # (sigma, tau, rho, rho*, s) -> (sigma, tau / k, rho / k, rho* / k, s / k) 
    # TODO implement 
    def raise_to_power(self, k):
        if not isinstance(k, int) or k < 2:
            raise ValueError("Parameter k must be an integer and >= 2.")
        raise NotImplementedError("TODO") # TODO

def literature_large_value_energy_region(region, ref, params=""):
    return Hypothesis(
        f"{ref.author()} large value energy region" + params,
        "Large value energy region",
        Large_Value_Energy_Region(region),
        f"See [{ref.author()}, {ref.year()}]",
        ref,
    )

def derived_large_value_energy_region(data, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived large value energy region",
        "Large value energy region",
        data,
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound

# Returns a Hypothesis object representing raising the large value energy region 
# to a integer power k. 
def get_raise_to_power_hypothesis(k):
    # name, hypothesis_type, data, proof, reference
    name = f"Large value energy region raise to power hypothesis with k = {k}"
    
    def f(h):
        region = copy.copy(h.data)
        region.raise_to_power(k)
        return derived_large_value_energy_region(
            region, 
            f"Follows from raising {h} to the k = {k} power",
            {h}
            )
    
    return Hypothesis(
        name,
        "Large value energy region transform",
        Transform(name, f),
        "Classical",
        Reference.classical(),
        )

# Convert an exponent pair (k, l) to a large value energy region 
# Lemma 10.15 
# s \leq 1 + max(
# rho + 1,
# 5/3 rho + tau / 3
# (2 + 3k + 4l) / (1 + 2k + 2l) rho + (k + l) / (1 + 2k + 2l) tau
# )
def ep_to_lver(eph):

    k, l = eph.data.k, eph.data.l
    
    # Construct Polytope of (sigma, tau, rho, rho*, s)
    polys = []

    # default limits 
    rect = Large_Value_Energy_Region.default_constraints()

    # 2 + rho - s >= 0
    polys.append(Polytope(rect + [[2, 0, 0, 1, 0, -1]]))

    # 1 + 1/3 * tau + 5/3 * rho - s >= 0
    polys.append(Polytope(rect + [[1, 0, frac(1,3), frac(5,3), 0, -1]]))

    # 1 + (k + l) / (1 + 2k + 2l) tau + (2 + 3k + 4l) / (1 + 2k + 2l) rho - s >= 0
    polys.append(Polytope(rect + [
        [
            1,
            0,
            frac(k + l, 1 + 2 * k + 2 * l),
            frac(2 + 3 * k + 4 * l, 1 + 2 * k + 2 * l),
            0,
            -1
        ]]))

    region = Region.union([Region(Region_Type.POLYTOPE, p) for p in polys])
    return derived_large_value_energy_region(
        Large_Value_Energy_Region(region),
        f"Follows from {eph.data}",
        {eph})

# Given a set of hypotheses and the choices of (sigma, tau), compute 
# the best available bound on LV*(sigma, tau) numerically. 
def approx_sup_LV_star(hypotheses, sigma, tau):
    lvers = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")

    # Compute intersection 
    E = Region(Region_Type.INTERSECT, [lver.data.region for lver in lvers])

    print(E)

    DE = E.to_disjoint_union()
    print(DE)
    
    # sigma, tau, rho, rho*, s
    sup = 0
    argmax = []
    for rho in np.linspace(0, 10, 100):
        for s in np.linspace(0, 10, 100):
            for rho_star in np.linspace(0, 10, 100):
                if sup < rho_star and E.contains([sigma, tau, rho, rho_star, s]):
                    sup = rho_star
                    argmax = [sigma, tau, rho, rho_star, s]
    
    print(sup)
    print(argmax)
    return sup

