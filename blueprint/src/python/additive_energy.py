# code implementation for large value additive energy region theorems

import copy 
from constants import Constants
from fractions import Fraction as frac
from functions import *
from hypotheses import Hypothesis, Hypothesis_Set
import large_values as lv
from reference import Reference
from region import Region, Region_Type
from transform import Transform
import zeta_large_values as zlv


class Additive_Energy_Estimate:
    
    """
    Class representing an additive energy estimate ρ* \\le LV*(σ, τ)
    """
    def __init__(self, region: Region):
        if not isinstance(region, Region):
            raise ValueError("Parameter region must be of type Region")
        self.region = region
        
    def __repr__(self):
        s = self.region.to_str(use_indentation=False, variables=["σ","τ","ρ*"])
        return f"(σ,τ,ρ*) in {s}"

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

def derived_additive_energy_estimate(
        region: Region, 
        proof: str, 
        deps: set[Hypothesis], 
        zeta: bool = False
    ) -> Hypothesis:

    """
    Construct a additive energy estimate Hypothesis object. 

    Parameters 
    ----------
    region : Region
        A subset of R^3 representing the set of feasible values 
        of (σ, τ, ρ*).
    proof : str
        The proof of this additive energy estimate as a string.
    deps : set of Hypothesis
        The set of hypothesis on which this result depends.
    zeta : bool, optional
        If True, the additive energy estimate is with respect to
        a zeta large value pattern. 
    
    Returns 
    -------
    Hypothesis
        The additive energy estimate as a Hypothesis. 
    """
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived additive energy estimate",
        "Additive energy estimate",
        Additive_Energy_Estimate(region),
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound

#############################################################################

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

def lv_to_lver(hypotheses: Hypothesis_Set, zeta: bool = False) -> list[Hypothesis]:

    """
    Converts all large value estimates to large value energy regions in a hypothesis
    set. 

    If zeta is True, then "Zeta large value estimate" and "Zeta large value energy 
    region" are considered instead.

    Parameters
    ----------
    hypotheses : Hypothesis_Set
        The set of hypotheses from which large value estimates are drawn.
    zeta : bool, optional
        Indicates whether zeta large value estimates should be used instead (default
        is False).

    Returns
    -------
    list of Hypothesis
        A list of hypotheses, each representing a large value energy region.
    """

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

def lver_to_lv(hypothesis: Hypothesis) -> Hypothesis:

    """
    Computes a large value region from a large value energy region by taking the 
    projection (σ, τ, ρ, ρ*, s) -> (σ, τ, ρ).

    Parameters
    ----------
    hypothesis: Hypothesis
        A hypothesis of type "Large value energy region" or "Zeta large value 
        energy region".
    
    Returns
    -------
    Hypothesis
        A hypothesis representing the computed large value region, of type "Large value 
        estimate".
    """

    if not isinstance(hypothesis, Hypothesis):
        raise ValueError("Parameter hypothesis must be of type Hypothesis")
    
    if hypothesis.hypothesis_type != "Large value energy region" and \
        hypothesis.hypothesis_type != "Zeta large value energy region":
        raise ValueError("Hypothesis must either be of type 'Large value energy region'" + \
                         " or 'Zeta large value energy region'")
    
    # Ensure that the region is represented as a union or disjoint union
    region = hypothesis.data.region.as_disjoint_union()

    # Project onto the dimensions (sigma, tau, rho)
    proj = region.project({0, 1, 2})
    proj.simplify()

    # Pack into Hypothesis object
    if hypothesis.hypothesis_type == "Large value energy region":
        return lv.derived_bound_LV(
            proj, 
            f"Follows from {hypothesis} and taking the projection (σ, τ, ρ, ρ*, s) -> (σ, τ, ρ).", 
            {hypothesis}
        )
    else:
        return zlv.derived_bound_zeta_LV(
            proj,
            f"Follows from {hypothesis} and taking the projection (σ, τ, ρ, ρ*, s) -> (σ, τ, ρ).", 
            {hypothesis} 
        )

def lver_to_energy(hypothesis: Hypothesis) -> Hypothesis:

    """
    Computes a additive energy region from a large value energy region by taking the 
    projection (σ, τ, ρ, ρ*, s) -> (σ, τ, ρ*). 

    Parameters
    ----------
    hypothesis: Hypothesis
        A hypothesis of type "Large value energy region" or "Zeta large value 
        energy region".
    
    Returns
    -------
    Hypothesis
        A hypothesis representing the computed additive energy region, of type "Large 
        value estimate".
    """

    if not isinstance(hypothesis, Hypothesis):
        raise ValueError("Parameter hypothesis must be of type Hypothesis")
    
    if hypothesis.hypothesis_type != "Large value energy region" and \
        hypothesis.hypothesis_type != "Zeta large value energy region":
        raise ValueError("Hypothesis must either be of type 'Large value energy region'" + \
                         " or 'Zeta large value energy region'")
    
    # Ensure that the region is represented as a union or disjoint union
    region = hypothesis.data.region.as_disjoint_union()

    # Project onto the dimensions (sigma, tau, rho*)
    proj = region.project({0, 1, 3})

    # Handle empty projections 
    if proj is None: return None
    
    proj.simplify()

    # Pack into Hypothesis object
    zeta = hypothesis.hypothesis_type == "Zeta large value energy region"
    return derived_additive_energy_estimate(
        proj, 
        f"Follows from {hypothesis} and taking the projection (σ, τ, ρ, ρ*, s) -> (σ, τ, ρ*).", 
        {hypothesis}
    )

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




