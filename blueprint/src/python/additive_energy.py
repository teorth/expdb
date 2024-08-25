# code implementation for large value additive energy region theorems
from functions import *
from polytope import *
import copy 

class Large_Value_Energy_Region:
    
    # Construct a valid large value energy region represented as a union of 
    # 5-dimensional polytopes 
    # Parameters:
    #   - polytopes (list of Polytope) the set of polytopes whose union represent
    #               the energy region. The polytopes may not necessarily be disjoint.
    def __init__(self, polytopes):
        self.region = polytopes
    
    def __repr__(self):
        return "Union of " + str(self.region)
        
        
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
        
















