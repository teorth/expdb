# code implementation for large value additive energy region theorems

import copy 
from constants import Constants
import exponent_pair as ep
from fractions import Fraction as frac
from functions import *
from hypotheses import Hypothesis, Hypothesis_Set
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

def expand_lver(hypotheses: Hypothesis_Set):
    """
    Add new hypothesis of type "Large value energy region" to a set of 
    hypotheses by applying transformations to existing large value energy
    regions. 

    Note: these transformations only apply to large value energy regions
    (and not zeta large value energy regions).
    """
    lvers = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")
    tfs = hypotheses.list_hypotheses(hypothesis_type="Large value energy region transform")
    for tf in tfs:
        hypotheses.add_hypotheses(tf.data.transform(lver) for lver in lvers)

# Convert an exponent pair (k, l) to a large value energy region
def ep_to_lver(eph):

    k, l = eph.data.k, eph.data.l
    
    # Terms of the first maximum in the form
    # [constant, sigma, tau, rho, rho^*, s] 
    first_max = [
        [1, 0, 0, 1, 0, 0],
        [0, 0, frac(1,3), frac(5,3), 0, 0],
        [0, 0, (k + l) / (1 + 2*k + 2*l), (2 + 3*k + 4*l) / (1 + 2*k + 2*l), 0, 0]
    ]
    # Terms of the second maximum
    second_max = [
        [1, 0, 0, 0, 1, 0],
        [0, 0, 0, 4, 0, 0],
        [0, 0, frac(1,2), 1, frac(3,4), 0]
    ]

    # Iterate through all combinations of first and second maximum
    constraints = []
    for max1 in first_max:
        for max2 in second_max:
            # additional terms
            add = [1, -2, 0, 0, -1, 0]
            constraint = [add[i] + frac(1,2) * (max1[i] + max2[i]) for i in range(6)] # constraint >= 0
            # Solve for \rho^* (this step is not strictly necessary - only for presentation)
            factor = -frac(1, constraint[4])
            constraint = [c * factor for c in constraint]
            constraints.append(constraint)

    region = Region.from_union_of_halfplanes(
        constraints, 
        Large_Value_Energy_Region.default_constraints()
    )
    return derived_large_value_energy_region(
        Large_Value_Energy_Region(region),
        f"Follows from {eph.data}",
        {eph})


# Given a Hypothesis_Set, find Hypothesis of type "Large value estimate" or "Zeta large value estimate", 
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
        # Lift (sigma, tau, rho) -> (sigma, tau, rho, rho*, s)
        region = lvh.data.region.lift([
            0, 
            1, 
            2, 
            (0, Constants.LV_DEFAULT_UPPER_BOUND),
            (0, Constants.LV_DEFAULT_UPPER_BOUND)
        ])
        hyps.append(
            constructor(
                Large_Value_Energy_Region(region),
                f"Follows from {lvh}",
                {lvh}
            )
        )
    return hyps

# Given a Hypothesis_Set, find all Hypothesis of type "Large value energy region", computes 
# their intersection as a (simplified) polytope, then projects onto the (sigma, tau, rho)
# domain, returning the result as Region in R^3
# TODO: when the large value estimates are converted to their Polytope representation, 
# change the return type of this function to Hypothesis. 
def lver_to_lv(hypotheses, zeta=False):

    if zeta:
        hyps = hypotheses.list_hypotheses(hypothesis_type="Zeta large value energy region")
    else:
        hyps = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")
    
    # Take the energy region (sigma, tau, rho, rho*, s)
    E = Region.intersect([h.data.region for h in hyps])
    E1 = E.as_disjoint_union()

    # Project onto the dimensions (sigma, tau, rho)
    proj = E1.project({0, 1, 2})
    proj.simplify()
    return proj


def compute_best_lver(
        hypotheses: Hypothesis_Set, 
        sigma_tau_domain: Region, 
        zeta: bool = False, 
        debug: bool = False
    ) -> Hypothesis:

    """
    Computes the best large value energy region (i.e. set of feasible tuples of 
    (sigma, tau, rho, rho, s)) implied by the hypotheses of a given hypothesis set. 

    Parameters
    ----------
    hypotheses: Hypothesis_Set
        The set of hypotheses to consider, e.g. of type "Large value energy region" 
        or "Large value energy region transform".
    sigma_tau_domain: Region
        The domain of (sigma, tau) to consider. 
    zeta: bool, optional
        If True, the zeta large value energy region will be computed instead (default
        is False).
    debug: bool, optional
        If True, additional debugging info will be logged to console (default is False).
    
    Returns
    -------
    Hypothesis
        A Hypothesis representing the computed large value energy region.
    """

    # 1. A large value energy region is also a zeta large value energy region
    # 2. Large value energy regions have the raise to power hypothesis, which is 
    # a type of Large value energy region transform
    # 3. Zeta large value energy regions do not have the raise to power hypothesis
    # Therefore, the set of zeta large value energy regions can be obtained by 
    # first expanding the set of large value energy regions (by raising to a power)
    # then adding in the additional zeta large value energy regions later. 
    lvers = hypotheses.list_hypotheses(hypothesis_type="Large value energy region")
    tfs = hypotheses.list_hypotheses(hypothesis_type="Large value energy region transform")
    transformed = []
    for tf in tfs: transformed.extend(tf.data.transform(lver) for lver in lvers)
    lvers.extend(transformed)

    if zeta:
        # A large value energy region is also a zeta large value energy region
        # A zeta large value energy region has no raise to power hypothesis
        lvers.extend(hypotheses.list_hypotheses(hypothesis_type="Zeta large value energy region"))
    
    # Compute intersection over the domain
    domain = sigma_tau_domain.lift([
                0, 
                1, 
                (0, Constants.LV_DEFAULT_UPPER_BOUND),
                (0, Constants.LV_DEFAULT_UPPER_BOUND),
                (0, Constants.LV_DEFAULT_UPPER_BOUND)
            ])

    # Compute the large value energy bounding region
    E = Region(Region_Type.INTERSECT, [domain] + [lver.data.region for lver in lvers])

    E1 = E.as_disjoint_union(verbose=debug)

    # Pack into Hypothesis object
    proof = f"Follows from taking the intersection of {len(lvers)} large value energy regions"
    # TODO Keep track of which hypotheses are required to generate the final region. 
    dependencies = set(lvers)
    if zeta:
        return derived_zeta_large_value_energy_region(
            Large_Value_Energy_Region(E1), proof, dependencies)
    else:
        return derived_large_value_energy_region(
            Large_Value_Energy_Region(E1), proof, dependencies)


# Given a set of hypotheses, compute the best available bound on LV*(sigma, tau)
# as a polytope in R^3 with dimensions (sigma, tau, rho*) for (sigma, tau) \in sigma_tau_domain 
# (represented as a Polytope)
# If zeta is true, the best available bound on LV*_\zeta(sigma, tau) is computed instead 
def compute_LV_star(hypotheses, sigma_tau_domain, debug=True, zeta=False):

    hyp = compute_best_lver(hypotheses, sigma_tau_domain, zeta=zeta, debug=debug)
    E = hyp.data.region
    deps = hyp.dependencies

    # Project onto the (sigma, tau, rho*) dimension
    Eproj = E.project({0, 1, 3})
    
    if Eproj is None: return None
    
    Eproj.simplify()

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




