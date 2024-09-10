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
    # (sigma, tau, rho, rho*, s) -> (sigma, k * tau, k * rho, k * rho*, k * s) 
    # TODO implement 
    def raise_to_power(self, k):
        if not isinstance(k, int) or k < 2:
            raise ValueError("Parameter k must be an integer and >= 2.")
        raise NotImplementedError("TODO") # TODO

#############################################################################


# Computes a Region object representing the union of a list of polytopes, 
# where each polytope is defined as the intersection of 
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
    region = Large_Value_Energy_Region.union_of_halfplanes(
        [
            [2, 0, 0, 1, 0, -1],                    # 2 + rho - s >= 0
            [1, 0, frac(1,3), frac(5,3), 0, -1],    # 1 + 1/3 * tau + 5/3 * rho - s >= 0
            [   # 1 + (k + l) / (1 + 2k + 2l) tau + (2 + 3k + 4l) / (1 + 2k + 2l) rho - s >= 0
                1,
                0,
                frac(k + l, 1 + 2 * k + 2 * l),
                frac(2 + 3 * k + 4 * l, 1 + 2 * k + 2 * l),
                0,
                -1
            ]
        ],
        rect
    )
    return derived_large_value_energy_region(
        region,
        f"Follows from {eph.data}",
        {eph})

# Given a Hypothesis_Set, convert all Hypothesis objects of type "Large value estimate" into 
# large value energy regions and returns them
def lv_to_lver(hypotheses):
    lvs = hypotheses.list_hypotheses(hypothesis_type="Large value estimate")

    # A large value estimate currently is represented as a 2-dimensional affine function rho
    # <= f(sigma, tau). Convert this into a polytope representing the set of feasible rho
    # values in (sigma, tau, rho, rho*, s) space, with default limits on the unconstrained
    # variables rho* and s.
    hyps = []
    for lv in lvs:
        polys = []
        # Each piece is an affine function of (sigma, tau)
        for piece in lv.data.bound.pieces:
            # Express this piece as a polytope 
            # Lift (sigma, tau) -> (sigma, tau, rho, rho*, s)
            P = piece.domain.lift([
                0, 
                1,
                (0, Constants.LV_DEFAULT_UPPER_BOUND),
                (0, Constants.LV_DEFAULT_UPPER_BOUND),
                (0, Constants.LV_DEFAULT_UPPER_BOUND)
            ]).intersect(
                # rho <= f[0] + f[1] * sigma + f[2] * tau
                Polytope([
                    [piece.a[0], piece.a[1], piece.a[2], -1, 0, 0]
                ])
            )
            polys.append(Region(Region_Type.POLYTOPE, P))
        region = Region(Region_Type.DISJOINT_UNION, polys)
        hyps.append(
            derived_large_value_energy_region(
                Large_Value_Energy_Region(region),
                f"Follows from {lv}",
                {lv}
            )
        )
    return hyps


import random as rd
rd.seed(1007)

# Debugging method to check whether two regions agree
def sample_check(region1, region2, N=1000, info=None):
    ntrues = 0
    npassed = 0
    for i in range(N):
        x = (rd.uniform(1/2, 1), rd.uniform(0, 5), rd.uniform(0, 5), rd.uniform(0, 5), rd.uniform(0, 5))
        c1 = region1.contains(x)
        c2 = region2.contains(x)
        if c1 != c2:
            print(i, x)
            print(info)
            raise ValueError()
        else:
            npassed += 1
        if c1:
            ntrues += 1
    print(f"[Debug info] Checking regions equal. Passed: {npassed}/{N}", "Contained:", ntrues)

# Given a set of hypotheses, compute the best available bound on LV*(sigma, tau)
# as a polytope in R^3 with dimensions (sigma, tau, rho*)
def compute_LV_star(hypotheses, debug=True):
    lvers = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")
    
    # Compute intersection 
    E = Region(Region_Type.INTERSECT, [lver.data.region for lver in lvers])
    E1 = E.as_disjoint_union()

    # if debug: randomly sample some points, and test inclusion/exclusion 
    if debug:
        sample_check(E, E1, N=10000, info=lvers)
    
    # Project onto the (sigma, tau, rho*) dimension
    Eproj = E1.project({0, 1, 3})
    return Eproj




