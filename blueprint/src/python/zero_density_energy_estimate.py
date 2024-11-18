# code for bounding the function A*(\sigma) in the estimate
#
# N*(\sigma - \delta, T) \ll T^{A*(\sigma)(1 - \sigma) + o(1)}
#
# for all \delta > 0, where N(\sigma, T) is the additive energy of 
# the imaginary parts of the zeroes \rho of the Riemann zeta-function 
# satisfying with real part \in [\sigma, 1] and imaginary part \in 
# [-T, T]. 

import additive_energy as ad
from fractions import Fraction as frac
from functions import Affine, Interval, RationalFunction as RF
from hypotheses import Hypothesis, Hypothesis_Set
from polytope import Polytope
from reference import Reference
from region import Region


class Zero_Density_Energy_Estimate:

    """
    Class representing a zero-density energy estimate 

    For now, this is identical to the Zero_Density_Estimate class 
    TODO: set up an inheritance structure for these classes
    """
    
    def __init__(self, expr: str, interval: Interval):
        """
        Parameters
        ----------
        expr : str
            An expression that represents a function of x.
        interval : Interval 
            The Interval object representing the range of validity of the
            bound.
        """
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

#############################################################################
# Useful methods to construct Hypothesis instances of type "Zero density 
# energy estimate"

def literature_zero_density_energy_estimate(
        estimate: str, 
        interval: Interval,
        ref: Reference, 
        params: str = ""
    ) -> Hypothesis:

    """
    Returns a Hypothesis object representing a zero density energy theorem from 
    the literature

    Parameters
    ----------
    estimate: str
        A univarate function of x representing a bound on A*(x)
    interval: Interval
        The range of sigma on which the bound is valid. 
    ref: Reference
        The literature reference for this bound.
    params: str, optional
        Additional parameters to include in the name of this Hypothesis object. 
        (Default is "").  
    
    Returns
    -------
    Hypothesis
        A Hypothesis object representing a zero density energy estimate from 
        the literature. 
    """

    return Hypothesis(
            f"{ref.author()} ({ref.year()}) zero density energy estimate" + params,
            "Zero density energy estimate",
            Zero_Density_Energy_Estimate(estimate, interval),
            f"See [{ref.author()}, {ref.year()}]",
            ref,
        )

def derived_zero_density_energy_estimate(
        data: Zero_Density_Energy_Estimate, 
        proof: str, 
        deps: set[Hypothesis]
    ) -> Hypothesis:

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

def trivial_zero_density_energy_estimate(
        data: Zero_Density_Energy_Estimate) -> Hypothesis:
    bound = Hypothesis(
        "Trivial zero density energy estimate",
        "Zero density energy estimate",
        data,
        "Trivial",
        Reference.trivial(),
    )
    return bound

#############################################################################
# Methods for converting between zero density energy theorems and theorems 
# on other mathematical objects.

def add_trivial_zero_density_energy_estimates(hypotheses: Hypothesis_Set):

    """
    Computes the "trivial" additive energy estimates implied by zero density 
    estimates from a set of hypothesis, and add them to the hypothesis set.
    """

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

def lver_to_energy_bound(
        hypotheses: Hypothesis_Set,
        tau0: Affine,
        debug: bool = False
    ) -> list[Hypothesis]:
    
    """
    Given a set of hypotheses, compute the best large value energy region
    then computes the best bound on the additive energy A^*(\\sigma).

    Parameters
    ----------
    hypotheses : Hypothesis_Set
        The set of hypothesis from which to construct the large value 
        energy region.
    debug : bool, optional
        If True, additional debugging information will be logged to console
        (default is False).
    """
    
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set.")
    if not isinstance(tau0, Affine):
        raise ValueError("Parameter tau0 must be of type Affine.")
    
    fns = []
    deps = set()

    # domain representing tau0 <= tau <= 2 tau0
    sigma_interval = tau0.domain
    LVER_domain = Region.from_polytope(
        Polytope([
            [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
            [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
            [-tau0.c, -tau0.m, 1],       # tau >= tau0 = m sigma + c
            [2 * tau0.c, 2 * tau0.m, -1] # tau <= 2 tau0 = 2 m sigma + 2 c
        ])
    )
    lver = ad.compute_best_lver(hypotheses, LVER_domain, zeta=False, debug=debug)
    energy = ad.lver_to_energy(lver)
    if energy is not None:
        sup1 = compute_sup_LV_on_tau(energy.data.region, sigma_interval)
        fns.extend(list(sup1))
        deps.add(energy)
    
    # domain representing 2 <= tau <= tau0
    LVER_zeta_domain = Region.from_polytope(
        Polytope([
            [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
            [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
            [-2, 0, 1],                  # tau0 >= 2
            [tau0.c, tau0.m, -1],        # tau <= tau0 = m sigma + c
        ])
    )
    zlver = ad.compute_best_lver(hypotheses, LVER_zeta_domain, zeta=True, debug=debug)
    energy_zeta = ad.lver_to_energy(zlver)
    if energy_zeta is not None:
        sup2 = compute_sup_LV_on_tau(energy_zeta.data.region, sigma_interval)
        fns.extend(list(sup2))
        deps.add(energy_zeta)
    
    # Take maximum
    bounds = [(f[0], f[1]) for f in fns]
    sup = RF.max(bounds, sigma_interval)
    return [
        derived_zero_density_energy_estimate(
            Zero_Density_Energy_Estimate.from_rational_func(s[0], s[1]),
            f"Follows from combining {len(lver.dependencies)} large value energy regions " + \
            f"and {len(zlver.dependencies)} zeta large value energy regions",
            deps
        )
        for s in sup
        ]

def compute_sup_LV_on_tau(
        LV_region: Region, 
        sigma_interval: Interval
    ) -> list[tuple[RF, Interval, int]]:
    
    """
    Given a set of (sigma, tau, rho) values, compute the supremum of (rho / tau)
    as sigma ranges in an interval. 

    The result is expressed as a function of sigma.

    Parameters
    ----------
    LV_region : Region
        The region of feasible (sigma, tau, rho) values. This region must be a
        union or disjoint union of Polytope objects.
    sigma_interval : Interval
        The range of sigma values to consider. 

    Returns
    -------
    list of tuple of RationalFunction, Interval and int
        A list of piecewise-defined univariate functions of sigma. 
    """
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

def compute_best_energy_bound(hypotheses: Hypothesis_Set) -> list:

    """
    Given a Hypothesis_Set, compute a piecewise function representing the best bound 
    on A*(\\sigma)

    TODO: this function is quite similar to compute_best_density_estimate in zero_density_estimate.py
    Combine/merge them?
    """

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
