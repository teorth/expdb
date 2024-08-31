# code implementation for large value additive energy region theorems

import copy 
from constants import Constants
from fractions import Fraction as frac
from functions import *
from hypotheses import Hypothesis
from polytope import Polytope
from reference import Reference
from transform import Transform

class Large_Value_Energy_Region:
    
    # Construct a valid large value energy region represented as a union of 
    # 5-dimensional polytopes 
    # Parameters:
    #   - polytopes (list of Polytope or Polytope) the set of polytopes whose union 
    #     represent the energy region. The polytopes may not necessarily be disjoint.
    def __init__(self, polytopes):
        # Should this be a set instead of a list?
        if not isinstance(polytopes, list) and \
            not isinstance(polytopes, Polytope):
            raise ValueError("Parameter polytopes must be a list of Polytope or Polytope")
            
        if isinstance(polytopes, Polytope):
            self.region = [polytopes]
        else:
            self.region = polytopes
    
    
    def __repr__(self):
        return "Union of " + str(self.region)
        
    def __copy__(self):
        return Large_Value_Energy_Region(copy.copy(self.region))
    
    # Static methods ---------------------------------------------------------

    # The default bounds on the tuple (sigma, tau, rho, rho*, s). This is to ensure
    # that all large value regions are finite regions. 
    def default_region():
        return Polytope.rect(
            (frac(1,2), frac(1)),
            (0, Constants.TAU_UPPER_LIMIT),
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND))


    # ------------------------------------------------------------------------
    
    # Returns whether the region contains a 5-dimensional point 
    def contains(self, point):
        if len(point) != 5: 
            raise ValueError("point must be 5-dimensional")
        # Slow method - simply loop through the polytopes 
        return any(poly.contains(point) for poly in self.region)
    
    
    # Computes the union of two large value regions (this object is modified )
    def union_with(self, other):
        if not isinstance(other, Large_Value_Energy_Region):
            raise ValueError("Parameter other must be of type Large_Value_Energy_Region")
        self.region.extend(other.region)


    # Raise this region to the kth power
    # (sigma, tau, rho, rho*, s) -> (sigma, tau / k, rho / k, rho* / k, s / k) 
    def raise_to_power(self, k):
        if not isinstance(k, int) or k < 2:
            raise ValueError("Parameter k must be an integer and >= 2.")
        new_regions = []
        for r in self.region:
            r_cpy = copy.copy(r)
            for i in range(1,5):
                r_cpy.scale(i, frac(1,k), [])
            new_regions.append(r_cpy)
        return Large_Value_Energy_Region(new_regions)


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




print(Large_Value_Energy_Region.default_region())







