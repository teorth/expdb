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
        polys.append(
            Region(
                Region_Type.POLYTOPE, 
                Polytope(box + [hp] + neg_ineq, canonicalize=True)
            )
        )
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

def literature_zeta_large_value_energy_region(region, ref, params=""):
    return Hypothesis(
        f"{ref.author()} zeta large value energy region" + params,
        "Zeta large value energy region",
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

def derived_zeta_large_value_energy_region(data, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived zeta large value energy region",
        "Zeta large value energy region",
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
        region = h.data.region.scale_all([frac(1), frac(k), frac(k), frac(k), frac(k)])
        return derived_large_value_energy_region(
            Large_Value_Energy_Region(region), 
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


# Given a list of Hypothesis objects of type "Large value estimate" or "Zeta large value estimate", 
# convert them into large value energy regions and returns them as a list of Hypothesis
#
# If zeta is True, then "Zeta large value estimate" and "Zeta large value energy region" are
# considered instead 
def lv_to_lver(hypotheses, zeta=False):
    
    if zeta:
        lvs = hypotheses.list_hypotheses(hypothesis_type="Zeta large value estimate")
        constructor = derived_zeta_large_value_energy_region
    else:
        lvs = hypotheses.list_hypotheses(hypothesis_type="Large value estimate")
        constructor = derived_large_value_energy_region

    # A large value estimate currently is represented as a 2-dimensional affine function rho
    # <= f(sigma, tau). Convert this into a polytope representing the set of feasible rho
    # values in (sigma, tau, rho, rho*, s) space, with default limits on the unconstrained
    # variables rho* and s.
    hyps = []
    for lvh in lvs:
        polys = []
        # Each piece is an affine function of (sigma, tau)
        for piece in lvh.data.bound.pieces:
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
            constructor(
                Large_Value_Energy_Region(region),
                f"Follows from {lvh}",
                {lvh}
            )
        )
    return hyps


import random as rd
rd.seed(1007)

# Debugging method to check whether two regions agree
def sample_check(region1, region2, N=1000, dim=5, info=None):
    ntrues = 0
    npassed = 0
    for i in range(N):
        # Supports 3 and 5-dimensional tests
        if dim == 5:
            x = (rd.uniform(1/2, 1), rd.uniform(0, 5), rd.uniform(0, 5), rd.uniform(0, 5), rd.uniform(0, 5))
        elif dim == 3:
            x = (rd.uniform(1/2, 1), rd.uniform(0, 5), rd.uniform(0, 5))
        else:
            raise NotImplementedError()

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
# as a polytope in R^3 with dimensions (sigma, tau, rho*) for (sigma, tau) \in sigma_tau_domain 
# (represented as a Polytope)
# If zeta is true, the best available bound on LV*_\zeta(sigma, tau) is computed instead 
def compute_LV_star(hypotheses, sigma_tau_domain, debug=True, zeta=False):

    # 1. A large value energy region is also a zeta large value energy region
    # 2. Large value energy regions have the raise to power hypothesis, which is 
    # a type of Large value energy region transform
    # 3. Zeta large value energy regions do not have the raise to power hypothesis
    # Therefore, the set of zeta large value energy regions can be obtained by 
    # first expanding the set of large value energy regions (by raising to a power)
    # then adding in the additional zeta large value energy regions later. 
    lvers = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")
    
    # Use LVER transformations and use them to expand the set of LVERs
    transforms = hypotheses.list_hypotheses(hypothesis_type="Large value energy region transform")
    transformed_lvers = []
    for tf in transforms:
        transformed_lvers.extend(tf.data.transform(lver) for lver in lvers)
    lvers.extend(transformed_lvers)

    h1 = hypotheses.find_hypothesis(keywords="Heath-Brown large value energy region 2a")
    x = [frac(5, 6), 3, 1, frac(31, 12), 1000000]
    x2 = [frac(5,6), frac(3,2), frac(1,2), frac(31,24), 100]
    print("HERE", h1.data.region.contains(x), h1)
    print("This poly contains ", x2)
    for i in range(len(h1.data.region.child)):
        poly= h1.data.region.child[i]
        if poly.contains(x2):
            print(i, poly)
    for tr in transforms:
        th = tr.data.transform(h1)
        print("HERE", th.data.region.contains(x), tr)
        if "k = 2" in tr.name:
            for poly in th.data.region.child:
                if poly.contains(x):
                    print(poly)

    if zeta:
        # A large value energy region is also a zeta large value energy region
        # A zeta large value energy region has no raise to power hypothesis
        lvers.extend(hypotheses.list_hypotheses(hypothesis_type="Zeta large value energy region"))
        print(f"Found {len(lvers)} zeta large value energy regions")
    
    if debug:
        for lver in lvers:
            print(lver, ":", lver.proof)
        
        print("sigma-tau domain:")
        print(sigma_tau_domain)
    # Compute intersection over the domain
    domain = Region.from_polytope(
        sigma_tau_domain.lift([
            0, 
            1, 
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND)
        ])
    )

    # Compute the large value energy bounding region
    E = Region(Region_Type.INTERSECT, [domain] + [lver.data.region for lver in lvers])

    E1 = E.as_disjoint_union()

    # if debug: randomly sample some points, and test inclusion/exclusion 
    if debug:
        sample_check(E, E1, N=10000, dim=5, info=lvers)
    
    # -------------------------------------------------------------------------
    # Keep track of which hypotheses are required to generate the final region. 
    # 
    # The ultimate goal is that whenever a {Polytope} -> Polytope operation 
    # (e.g. union, intersection) takes place, we keep track of the minimal set 
    # of polytopes that are necessary for determining the resulting polytope. 
    # E.g. during the union of polytopes p1, p2, if p1 is a subset of p2 then 
    # the set of minimal polytopes is {p2}. 
    # 
    # Unfortunately this is likely to be computationally expensive, so instead 
    # we perform a crude match after E has already been computed. Those LVERs 
    # which share a (non-trivial) facet with E are included in its minimal 
    # dependency set.

    def identify_constraint(con):
        for lver in lvers:
            polys = [r.child for r in lver.data.region.child]
            for poly in polys:
                if con in set(poly.get_constraints()):
                    return lver
        return None
    
    # Iterate through the non-trivial facets of E1 to work out the min. dep. set 
    trivial_cons = set(domain.child.get_constraints())
    for r in E1.child:
        # Get non-trivial constraints of the current polytope
        cons = set(
            c for c in r.child.get_constraints() if c.inequality_type != Constraint.EQUALS
        )
        cons.discard(trivial_cons)
        # Find out where those constraints came from
        deps = set()
        for con in cons:
            h = identify_constraint(con)
            if h is not None:
                deps.add(h)
        r.dependencies = deps
    
    # -------------------------------------------------------------------------
    
    # Project onto the (sigma, tau, rho*) dimension
    Eproj = E1.project({0, 1, 3})
    
    if debug:
        cpy = copy.copy(Eproj)
        
    Eproj.simplify()

    if debug:
        print("Simplifying:", len(cpy.child), "->", len(Eproj.child), Eproj)
        sample_check(Eproj, cpy, N=10000, dim=3, info=lvers)

    if zeta:    
        return derived_zeta_large_value_energy_region(
            Large_Value_Energy_Region(Eproj), 
            f"Follows from taking the intersection of {len(deps)} zeta large value energy regions" + \
            " then projecting onto the (sigma, tau, rho*) dimensions",
            deps
        )
    else:
        return derived_large_value_energy_region(
            Large_Value_Energy_Region(Eproj), 
            f"Follows from taking the intersection of {len(deps)} large value energy regions" + \
            " then projecting onto the (sigma, tau, rho*) dimensions",
            set(deps)
        )




