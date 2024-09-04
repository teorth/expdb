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

    # Computes the union of n large value energy regions
    def union(*regions):
        return Large_Value_Energy_Region(Region.union([r.region for r in regions]))
        
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

def literature_large_value_energy_region(region, ref):
    return Hypothesis(
        f"{ref.author()} large value energy region",
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


def small_tau_relation():
    # divide into two cases: tau <= 1 and tau >= 1
    rect1 = Large_Value_Energy_Region.default_constraints()
    rect1.append([1, 0, -1, 0, 0, 0]) # tau <= 1

    rect2 =  Large_Value_Energy_Region.default_constraints()
    rect2.append([-1, 0, 1, 0, 0, 0]) # tau >= 1

    # No additional constraints on the second region with tau >= 1
    polys = [Polytope(rect2)]
    
    # Region 1 polytopes
    # 2 + rho - s >= 0
    polys.append(Polytope(rect1 + [[2, 0, 0, 1, 0, -1]]))
    # 1 + 2rho - s >= 0
    polys.append(Polytope(rect1 + [[1, 0, 0, 2, 0, -1]]))
    region = Region.union([Region(Region_Type.POLYTOPE, p) for p in polys])
    return region


# Implementation of "Simplified Heath-Brown relation" Corollary 10.19
# TODO: move this over to literature.py
def heath_brown_relation():
    # divide into two cases: Region 1 with tau < 3/2 and 
    # region 2 with tau > 3/2
    rect1 = Large_Value_Energy_Region.default_constraints()
    rect1.append([frac(3,2), 0, -1, 0, 0, 0]) # tau <= 3/2

    rect2 =  Large_Value_Energy_Region.default_constraints()
    rect2.append([-frac(3,2), 0, 1, 0, 0, 0]) # tau >= 3/2

    # No additional constraints on the second region with tau >= 3/2
    polys = [Polytope(rect2)]
    
    # Region 1 polytopes
    # 1 - 2sigma + 3rho - rho* >= 0
    polys.append(Polytope(rect1 + [[1, -2, 0, 3, -1, 0]]))
    # 4 - 4sigma + rho - rho* >= 0
    polys.append(Polytope(rect1 + [[4, -4, 0, 1, -1, 0]]))
    # 3/2 - 2sigma + 5/2rho - rho* >= 0
    polys.append(Polytope(rect1 + [[frac(3,2), -2, 0, frac(5,2), -1, 0]]))
    region = Region.union([Region(Region_Type.POLYTOPE, p) for p in polys])
    return region

# Implementation of second Heath-Brown relation (Lemma 10.20)
# TODO: move this over to literature.py
def heath_brown_relation2(k)
    if not isinstance(k,int):
        raise ValueError("Parameter k must be an integer.")
    if k < 2:
        raise ValueError("Parameter k must be >= 2.")
    
    rect = Large_Value_Energy_Region.default_constraints()
    polys = []
    # 2 - 2sigma - rho >= 0
    polys.append(Polytope(rect + [[2, -2, 0, -1, 0, 0]]))
    # 3k/4 - k sigma + tau/4 - rho + rho*/4 >= 0
    polys.append(Polytope(rect + [[frac(3*k,4), -k, frac(1,4), -1, frac(1,4), 0]]))
    # 3k/4 - k sigma - rho + rho*/4 >= 0
    polys.append(Polytope(rect + [[frac(3*k,4), -k, 0, -1, frac(1,4), 0]])) 
    
    region = Region.union([Region(Region_Type.POLYTOPE, p) for p in polys])
    return region

def guth_maynard_relation(k):
    if not isinstance(k, int):
        raise ValueError("Parameter k must be an integer.")
    if k < 1:
        raise ValueError("Parameter k must be >= 1.")

    rect = Large_Value_Energy_Region.default_constraints()
    polys = []
    # 2 - 2sigma - rho >= 0
    polys.append(Polytope(rect + [[2, -2, 0, -1, 0, 0]]))
    # 1 - 2sigma + S_1/3 - rho >= 0 
    # -7/3 - 2sigma - rho >= 0
    polys.append(Polytope(rect + [[-frac(7,3), -2, 0, -1, 0, 0]]))
    
    # TODO: complete the list of constraints
    raise NotImplementedError()

def guth_maynard_relation2():
    rect = Large_Value_Energy_Region.default_constraints()
    return Region(Region_Type.POLYTOPE, Polytope(
        rect + [[0, -2, 0, 1, -1, 1]] # -2sigma + rho - rho* - s >= 0
    ))

def guth_maynard_relation3():
    rect1 = Large_Value_Energy_Region.default_constraints()
    rect1.append([frac(4,3), 0, -1, 0, 0, 0]) # 1 <= tau <= 4/3
    rect1.append([-1, 0, 1, 0, 0, 0]) 

    rect2 = Large_Value_Energy_Region.default_constraints()
    rect2.append([-frac(4,3), 0, 1, 0, 0, 0]) # tau >= 4/3

    rect3 = Large_Value_Energy_Region.default_constraints()
    rect3.append([1, 0, -1, 0, 0, 0]) # tau <= 1
    
    # No additional constraints for second and third region
    polys = [Polytope(rect2), Polytope(rect3)]

    # Impose constraints for first region
    # 4 - 4sigma + rho - rho* >= 0
    polys.append(Polytope(rect + [[4, -4, 0, 1, -1, 0]]))
    # 1 - 2sigma + tau/4 + 21/8 rho - rho* >= 0 
    polys.append(Polytope(rect + [[1, -2, frac(1,4), frac(21,8), -1, 0]]))
    # 1 - 2sigma + 3rho - rho* >= 0
    polys.append(Polytope(rect + [[1, -2, 0, 3, -1, 0]]))

    region = Region.union([Region(Region_Type.POLYTOPE, p) for p in polys])
    return region
