# code for bounding the function A(\sigma) in the bound
#
# N(\sigma, T) \ll T^{A(\sigma)(1 - \sigma) + o(1)}
#
# where N(\sigma, T) is the number of zeroes of the Riemann zeta-function
# in the rectangle \sigma \le \Re z \le 1 and |\Im z| \le T.

from constants import Constants
import copy
import exponent_pair as ep
from fractions import Fraction as frac
from functions import Affine, Interval, RationalFunction as RF
from hypotheses import *
import large_values as lv
import math
import matplotlib.pyplot as plt
import numpy as np
from polytope import Polytope
from reference import Reference
from region import Region
import scipy
import zeta_large_values as zlv

class Zero_Density_Estimate:

    """
    Object representing an estimate on A(\\sigma) represented as a RationalFunction
    on an Interval.

    The life cycle of this object is:
    1) Creation with the parameters expr (string) and interval (Interval)
    2) As needed, lazily parse expr to create bound (RationalFunction) that handles
    evaluations
    3) Once bound is set, expr is only used for stringify methods

    When creating derived zero-density estimates, occasionally it is more convenient
    to initialise with the bound object directly. In that case, use the static
    from_rational_func function to create the an instance which sets expr to the
    default __str__ representation of RationalFunction, which may or may not coincide
    with the original expr used to generate the bound object in the first place.
    This should cause any circular reference problems (only potential display
    inconsistencies) since, once the bound object is initialised, the expr object
    is not used for computation.

    In the future, we will probably move to a canonical model with a guaranteed
    one-to-one correspondence between expr and bound. At present this is challenging
    since for each rational function there are multiple possible valid representations.

    """

    def __init__(self, expr, interval):

        """
        Constructor

        Parameters
        ----------
        expr : str
            A function of x representing the zero-density bound A(sigma) <= f(sigma).
        interval : Interval 
            The range of validity of the zero-density estimate.
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

    def from_rational_func(rf: RF, interval: Interval) -> 'Zero_Density_Estimate':
        if not isinstance(rf, RF):
            raise ValueError("Parameter rf must be of type RationalFunction")

        zde = Zero_Density_Estimate(str(rf), interval)
        zde.bound = rf
        return zde

###############################################################################

def derived_zero_density_estimate(data, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived zero density estimate",
        "Zero density estimate",
        data,
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound


def add_zero_density(hypotheses, estimate, interval, ref, params=""):
    """
    Add a zero-density estimate to the hypothesis set

    Parameters
    ----------
    hypotheses: Hypothesis
        The hypothesis set to add to. 
    estimate: str
        A function of x representing a bound on A(x).
    interval - the domain of \sigma for which the estimate holds
    """
    hypotheses.add_hypotheses(
        Hypothesis(
            f"{ref.author()} ({ref.year()}) zero density estimate" + params,
            "Zero density estimate",
            Zero_Density_Estimate(estimate, interval),
            f"See [{ref.author()}, {ref.year()}]",
            ref,
        )
    )

###############################################################################

def compute_large_value_region(
        hypotheses : Hypothesis_Set,
        domain : Region,
        zeta = False, 
        debug = False) -> Hypothesis:

    """
    Computes the feasible region of large value tuples (\\sigma, \\tau, \\rho) 
    implied by a hypothesis set, returning a Hypothesis object representing the 
    computed region.

    Parameters
    ----------
    hypotheses : Hypothesis_Set
        The set of hypothesis to assume.
    domain : Region
        A 3-dimensional region representing the range of (sigma, tau, rho) 
        values to consider. 
    zeta : bool, optional
        If True, a zeta large value region will be computed instead (default is 
        False).
    debug : bool, optional 
        If True, debugging information will be logged on the console (default 
        is False).

    Returns
    -------
    Hypothesis
        A hypothesis of type "Large value estimate" or "Zeta large value estimate"
        representing the smallest large value region obtainable from the set of 
        hypothesis. 
    """
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")
    if not isinstance(domain, Region):
        raise ValueError("Parameter domain must be of type Region.")
    
    # Get a list of all large value estimates from the hypothesis set
    hyps = hypotheses.list_hypotheses(hypothesis_type="Large value estimate")

    # Use LVR transformations and use them to expand the set of LVRs
    tfs = hypotheses.list_hypotheses(hypothesis_type="Large value estimate transform")
    transformed_lvs = []
    for tf in tfs:
        transformed_lvs.extend(tf.data.transform(h) for h in hyps)
    hyps.extend(transformed_lvs)

    if zeta:
        hyps.extend(hypotheses.list_hypotheses("Zeta large value estimate"))

    # Clip each hypothesis large value region to domain. Benefits of doing this 
    # (as opposed to taking an intersection with the domain region in the main 
    # intersection step) are:
    # 1) There tends to be many duplicate hypotheses after clipping - so we can 
    # identify them early (in a shallow fashion) and remove them
    # 2) We may copy region objects so that any labelling does not affect the 
    # Polytopes in the original hypotheses
    regions = []
    uniq_hash = set()
    for i in range(len(hyps)):
        r = Region.intersect([domain, hyps[i].data.region])
        r = r.as_disjoint_union(verbose=False, track_dependencies=False)
        r.simplify()
        # Simple hash to quickly remove most duplicates
        hsh = str(r) 
        if hsh in uniq_hash: 
            continue
        uniq_hash.add(hsh)
        # Create a copy of the current region object, then set label
        r = copy.copy(r)
        r.set_label(i)
        regions.append(r)

    # Compute the large value energy bounding region
    # LV_region = Region(Region_Type.INTERSECT, [domain] + [h.data.region for h in hyps])
    LV_region = Region.intersect(regions)
    LV_region = LV_region.as_disjoint_union(verbose=debug, track_dependencies=True)

    # Calculate dependency set using labels 
    deps = set()
    for r in LV_region.child:
        deps = deps | r.child.dependencies 
    deps = set(hyps[d] for d in deps if d >= 0)

    # Construct hypothesis object and return
    if zeta:
        return zlv.derived_bound_zeta_LV(
            LV_region, f"Follows from {len(deps)} zeta large value estimates", deps
        )
    else:
        return lv.derived_bound_LV(
            LV_region, f"Follows from {len(deps)} large value estimates", deps
        )

def compute_sup_rho_on_tau(polys: list, sigma_interval: Interval) -> list:

    """
    Given a set of feasible (sigma, tau, rho) tuples (represented as the union of 
    a list of Polytope objects in R^3), compute the supremum of rho / tau as a 
    function of sigma, as sigma varies in sigma_interval. 

    Parameters
    ----------
    polys : list of Polytope 
        The convex polytopes whose union represents the set of feasible 
        (sigma, tau, rho) values. 
    sigma_interval : Interval
        The range of values of sigma to consider. 

    Returns
    -------
    list of (RationalFunction, Interval) tuples 
        The supremum of rho / tau represented as a piecewise-defined function 
        of sigma.
    """

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
            if v1 + v2 in visited or v2 + v1 in visited: continue
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

def lv_zlv_to_zd(
        hypotheses : Hypothesis_Set, 
        sigma_interval : Interval, 
        tau0 : numbers.Number = frac(3)
    ) -> list[Hypothesis]:

    """
    Compute the best zero density estimate implied by a set of hypotheses 
    containing large value estimates and zeta large value estimates. This 
    function implements Corollary 11.7 in the web blueprint. Specifically, 
    this function computes the maximum of the two functions 

    \\sup_{tau \\in [\\tau_0, 2\\tau_0]} LV(sigma, tau) / tau

    \\sup_{tau \\in [2, \\tau_0]} LV_{\\zeta}(sigma, tau) / tau

    given a choice of \\tau_0, where LV(sigma, tau) and LV_{\\zeta}(sigma, tau)
    represent respectively the best large value estimate and zeta large 
    value estimate obtainable from hypotheses present in the given hypothesis
    set.

    Compared with the function lv_zlv_to_zd2(), there is less reliance on 
    making a good choice of \\tau_0 since the method uses raise-to-power 
    hypotheses to ensure that as long as \\tau_0 is taken to be sufficiently 
    large, the zero density estimate can be obtained. However the use of 
    raise-to-power hypotheses means this function is typically slower. 

    Parameters
    ----------
    hypotheses : Hypothesis_Set
        A set of hypothesis containing large value estimates and zeta 
        large value estimates.
    sigma_interval : Interval
        The range of sigma values to compute the maximum function on.
    tau0 : Number, optional
        The value of \\tau_0 to use in the supremum computation (default is 3).

    Returns 
    -------
    list of Hypothesis
        The derived zero density estimates. 
    """

    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set.")
    if not isinstance(sigma_interval, Interval):
        raise ValueError("Parameter sigma_interval must be of type Interval.")
    
    sigma_lim = (sigma_interval.x0, sigma_interval.x1)
    rho_lim = (0, Constants.LV_DEFAULT_UPPER_BOUND)

    # Get large value bounds in the interval tau0 <= tau <= 2 * tau0
    lvr = compute_large_value_region(
        hypotheses,
        Region.from_polytope(
            Polytope.rect(sigma_lim, (tau0, 2 * tau0), rho_lim)
        ),
        zeta = False
    )
    sup1 = compute_sup_rho_on_tau([r.child for r in lvr.data.region.child], sigma_interval)

    # Get zeta large value bounds in the interval 2 <= tau <= tau0
    hypotheses.add_hypotheses(zlv.compute_large_value_estimate(hypotheses))
    zlvr = compute_large_value_region(
        hypotheses,
        Region.from_polytope(
            Polytope.rect(sigma_lim, (frac(2), tau0), rho_lim)
        ),
        zeta = True
    )
    sup2 = compute_sup_rho_on_tau([r.child for r in zlvr.data.region.child], sigma_interval)

    bounds = [(s[0], s[1]) for s in sup1]
    bounds.extend((s[0], s[1]) for s in sup2)
    sup = RF.max(bounds, sigma_interval, track_dependencies=False)

    # pack into Hypothesis
    hyps = []
    for s in sup:
        proof = f'Follows from computed large value estimates (for σ in {sigma_interval}, ' + \
                f'τ in [{tau0},{2 * tau0}]) and zeta large value estimates (for σ in {sigma_interval}, τ in [{2},{tau0}])'
        deps = {lvr, zlvr}
        hyps.append(derived_zero_density_estimate(
            Zero_Density_Estimate.from_rational_func(s[0], s[1]),
            proof,
            deps
        ))
    return hyps


def lv_zlv_to_zd2(
        hypotheses:Hypothesis_Set, 
        sigma_interval:Interval, 
        tau0:Affine
    ) -> list[Hypothesis]:

    """
    Tries to prove the zero-density estimate using Corollary 11.8 by specifying the 
    value of tau0 to use and the range of sigma to consider. 

    Parameters
    ----------
    hypotheses : Hypothesis_Set
        The set of hypotheses to assume. 
    tau0 : RationalFunction 
        The choice of tau0 as a function of sigma 
    sigma_interval : Interval 
        The range of sigma values to consider
    """

    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set.")
    if not isinstance(sigma_interval, Interval):
        raise ValueError("Parameter sigma_interval must be of type Interval.")
    if not isinstance(tau0, Affine):
        raise ValueError("Parameter tau0 must be of type Affine.")
    
    sigma_rho_constraints = [
        [-sigma_interval.x0, 1, 0, 0],              # sigma >= x0
        [sigma_interval.x1, -1, 0, 0],              # sigma <= x1
        [0, 0, 0, 1],                               # rho >= 0
        [Constants.LV_DEFAULT_UPPER_BOUND, 0, 0, -1] # rho <= large number
    ]

    # Get large value bounds in the interval 2/3 tau0 <= tau <= tau0
    lvr = compute_large_value_region(
        hypotheses,
        Region.from_polytope(
            Polytope([
                [-frac(2,3) * tau0.c, -frac(2,3) * tau0.m, 1, 0], # tau >= 2/3 * tau0
                [tau0.c, tau0.m, -1, 0] # tau <= tau0
            ] + sigma_rho_constraints)
        ),
        zeta = False
    )
    sup1 = compute_sup_rho_on_tau([r.child for r in lvr.data.region.child], sigma_interval)

    # Get zeta large value bounds in the interval 2 <= tau <= 4/3 tau0
    hypotheses.add_hypotheses(zlv.compute_large_value_estimate(hypotheses))
    zlvr = compute_large_value_region(
        hypotheses,
        Region.from_polytope(
            Polytope([
                [-2, 0, 1, 0], # tau >= 2
                [frac(4,3) * tau0.c, frac(4,3) * tau0.m, 1, 0] # tau <= 4/3 * tau0
            ])
        ),
        zeta = True
    )
    sup2 = compute_sup_rho_on_tau([r.child for r in zlvr.data.region.child], sigma_interval)

    bounds = [(s[0], s[1]) for s in sup1]
    bounds.extend((s[0], s[1]) for s in sup2)
    sup = RF.max(bounds, sigma_interval, track_dependencies=False)

    # pack into Hypothesis
    hyps = []
    for s in sup:
        proof = "Follows from computed large value estimates and zeta large " + \
                f"value estimates with τ0 = {tau0.to_str('σ')})"
        deps = {lvr, zlvr}
        hyps.append(derived_zero_density_estimate(
            Zero_Density_Estimate.from_rational_func(s[0], s[1]),
            proof,
            deps
        ))
    return hyps

# Computes the zero-density estimate obtained from
#
# A(s) \leq 3m / ((3m - 2) s + 2 - m)
#
# for s \geq max {
#   (9m^2 - 4m + 2) / (12m^2 - 6m + 2)
#   (3m^2(1 + 2k + 2l) - (4k + 2l) m + 2k + 2l) / (4m^2(1 + 2k + 2l) - (6k + 4l) m + 2k + 2l)
# }
def ivic_ep_to_zd(exp_pairs, m=2):

    # Search only among the vertices of H
    dep = None
    sigma0 = 1
    for eph in exp_pairs:
        (k, l) = eph.data.k, eph.data.l
        v = (3 * m * m * (1 + 2 * k + 2 * l) - (4 * k + 2 * l) * m + 2 * k + 2 * l) / (
            4 * m * m * (1 + 2 * k + 2 * l) - (6 * k + 4 * l) * m + 2 * k + 2 * l
        )
        if v < sigma0:
            sigma0 = v
            dep = eph

    sigma0 = max(sigma0, frac(9 * m * m - 4 * m + 2, 12 * m * m - 6 * m + 2))
    sigma0 = min(sigma0, frac(6 * m * m - 5 * m + 2, 8 * m * m - 7 * m + 2))

    zde = Zero_Density_Estimate(f"{3*m}/({3*m-2}x + {2-m})", Interval(sigma0, 1))
    return derived_zero_density_estimate(
        zde, f"Follows from {dep.data}", {dep}
    )

def approx_bourgain_ep_to_zd(exp_pairs):

    sigmas = np.linspace(1/2, 1, 1000)

    points = [[p.data.k, p.data.l] for p in exp_pairs]
    conv = scipy.spatial.ConvexHull(np.array(points))
    vertices = [points[v] for v in conv.vertices]
    poly = Polytope.from_V_rep(vertices)

    for v in vertices:
        print(v[0], v[1])

    for s in sigmas:
        R1 = poly.intersect(Polytope([
            [0, 1, 0],
            [frac(11,85), -1, 0],
            [-frac(3,5), 0, 1],
            [1, 0, -1],
            [-13, 20, 15],
            [2 * s - 1, 2 * s, -1]
        ]))

        R2 = poly.intersect(Polytope([
            [-frac(11,85), 1, 0],
            [frac(1,5), -1, 0],
            [-frac(3,5), 0, 1],
            [1, 0, -1],
            [-13, 20, 15],
            [2 * s - 1, 2 * s, -1],
            [11 - 22 * s, 170 * s - 144, 11]
        ]))

        # iterate through the vertices of R1
        Abound = float('inf')
        argmin = (0, 1)
        if not R1.is_empty(include_boundary=False):
            verts = R1.get_vertices()
            for v in verts:
                (k, l) = v
                if 2 * (1 + k) * s - 1 - l <= 0: continue
                A = 4 * k / (2 * (1 + k) * s - 1 - l)
                if A < Abound:
                    Abound = A
                    argmin = (k, l)

        if not R2.is_empty(include_boundary=False):
            verts = R2.get_vertices()
            for v in verts:
                (k, l) = v
                if 2 * (1 + k) * s - 1 - l <= 0: continue
                A = 4 * k / (2 * (1 + k) * s - 1 - l)
                if A < Abound:
                    Abound = A
                    argmin = (k, l)

        print(s, argmin, Abound)


# Computes the zero-density estimate
#
# A(s) \leq 4k/(2(1 + k)s - 1 - l)   (s > s0)
#
# for any exponent pair (k, l) and number s0 jointly satisfying
#   - k < 11/85
#   - 11/85 < k < 1/5
#   - s0 = (144k - 11l - 11)/(170k - 22)
def bourgain_ep_to_zd():
    
    # Temporary:
    # For now, manually enter the exponent pairs found from the above numerical optimisation
    eps = [
        (frac(11,85), frac(59,85)),
        (frac(391, 4595), frac(3461, 4595)),
        (frac(2779, 38033), frac(58699, 76066)),
        (frac(1101653, 15854002), frac(12327829, 15854002)),
        (frac(1959, 47230), frac(3975, 4723)),
        (frac(1175779, 38456886), frac(16690288, 19228443)),
        (frac(89, 3478), frac(15327, 17390)),
        (frac(1, 100), frac(14, 15))
        ]
    
    bounds = []
    for (k, l) in eps:
        h = ep.derived_exp_pair(k, l, f"See proof of ({k}, {l})", set())
        if k <= frac(1, 5) and l >= frac(3, 5) and 15 * l + 20 * k >= 13:
            if k <= frac(11, 85):
                s0 = max(frac(1, 2), (l + 1) / (2 * (k + 1)))
                if s0 >= 1:
                    continue
                bounds.append(
                    (RF([4 * k], [2 * (1 + k), -1 - l]), Interval(s0, 1), h)
                )
            else:
                s0 = max(
                    frac(1, 2),
                    (l + 1) / (2 * (k + 1)),
                    frac(144 * k - 11 * l - 11, 170 * k - 22),
                )
                if s0 >= 1:
                    continue
                bounds.append(
                    (RF([4 * k], [2 * (1 + k), -1 - l]), Interval(s0, 1), h)
                )

    crits = set()
    for i in range(len(bounds)):
        for j in range(i):
            (b1, int1, h1) = bounds[i]
            (b2, int2, h1) = bounds[j]
            crits.update(b1.intersections(b2, int1.intersect(int2)))
            crits.update([int1.x0, int1.x1, int2.x0, int2.x1])

    crits = list(crits)
    crits.sort()

    soln = []
    for i in range(1, len(crits)):
        s1 = crits[i - 1]
        s2 = crits[i]
        if s2 == s1 or s1 < frac(1, 2):
            continue

        interval = Interval(s1, s2)
        s = interval.midpoint()  # the test point

        # Iterate through the faces, collecting t-coordinates where the vertical
        # line \sigma = s intersects with a face
        sup = float("inf")
        argmax = None
        intersecting_faces = [b for b in bounds if b[1].contains(s)]
        for b in intersecting_faces:
            q = b[0].at(s)
            if q < sup:
                sup = q
                argmax = b
        soln.append((argmax[0], interval, argmax[2]))

    # Simplify ranges by merging neighbouring intervals
    for i in range(len(soln) - 1, 0, -1):
        s1 = soln[i - 1]
        s2 = soln[i]
        if s1[0] == s2[0] and s1[1].x1 == s2[1].x0:
            # Merge
            s1[1].x1 = s2[1].x1
            soln.pop(i)

    # for s in soln:
    #     print(s[0], "for x \\in", s[1], s[1].contains(frac(15, 16)), s[2])
    #     s[2].recursively_list_proofs()
    return [derived_zero_density_estimate(
                Zero_Density_Estimate.from_rational_func(s[0], s[1]),
                f"Follows from applying [Bourgain, 1995] and taking {s[2]}", 
                {s[2]}
                ) 
            for s in soln]


# Compute all zero-density estimates derived from exponent pairs
def ep_to_zd(hypotheses):

    # TODO: this routine is used quite often, should we roll it into a separate method?
    hypotheses.add_hypotheses(
        ep.compute_exp_pairs(hypotheses, search_depth=5, prune=True)
    )
    hypotheses.add_hypotheses(ep.exponent_pairs_to_beta_bounds(hypotheses))
    hypotheses.add_hypotheses(ep.compute_best_beta_bounds(hypotheses))
    ephs = ep.beta_bounds_to_exponent_pairs(hypotheses)

    return bourgain_ep_to_zd(ephs) + [ivic_ep_to_zd(ephs, m=2)]


# Given a list of tuples (RationalFunction, interval), 
# returns a simplified list of tuples representing the same piecewise defined 
# function         
def simplify(pieces):
    # Simplify
    simplified_pieces = []
    i = 0
    while i < len(pieces):
        (fi, inti) = pieces[i]
        left = inti.x0
        right = inti.x1

        j = i + 1
        while j < len(pieces):
            (fj, intj) = pieces[j]
            if not (fi == fj and right == intj.x0):
                break

            right = intj.x1
            j += 1

        simplified_pieces.append((fi, Interval(left, right)))
        i = j
    return simplified_pieces
    
# Compute the set of density estimates implied by Pintz theorem using subdivision
def compute_pintz_density_estimate_subdiv():
    N = 20
    bounds = []
    for k in range(4, N):
        k_int = Interval(frac(1, k*(k+1)), frac(1,k*(k-1))) 
        # Bound 1: 4/(((k - 1) k - 4) x + (2 - k) k + 4)
        b1 = RF([4], [(k - 1) * k - 4, (2 - k) * k + 4])
        #b1 = RF([4], [(k - 1) * k, (2 - k) * k])
        for l in range(3, N):
            l_int = Interval(frac(1, 2*l*(l+1)), frac(1,2*l*(l-1)))
            # Bound 2: 3/((2 l - 2) l x + (3 - 2 l) l)
            b2 = RF([3], [(2 * l - 2) * l, (3 - 2 * l) * l])
            
            interval = k_int.intersect(l_int)
            sigma_interval = Interval(1 - interval.x1, 1 - interval.x0, 
                                      interval.include_upper, interval.include_lower)
            if not interval.is_empty():
                bounds.append((k, l, sigma_interval, [b1, b2]))
    
    # Compute the maximum 
    zdes = []
    for b in bounds:
        interval = b[2]
        (b1, b2) = b[3]
        zdes.extend(
            RF.max(
                [(b1, interval), (b2, interval)], 
                interval, 
                track_dependencies=False
            )
        )
    
    zdes = simplify(zdes)
    return zdes
    
# Aggregate the zero-density estimates in the Hypothesis_Set and returns a piecewise
# function that represents the best zero-density estimate in each subinterval of [1/2, 1]
# Note that this function does not compute zero-density estimates from other
# Hypothesis objects - it only aggregates existing Hypothesis objects of type
# "Zero density estimate".
def best_zero_density_estimate(hypotheses, verbose=False):
    hs = hypotheses.list_hypotheses(hypothesis_type="Zero density estimate")

    # Ensure bound is computed (i.e any RationalFunction objects that use lazy initialization
    # is parsed)
    for h in hs:
        h.data._ensure_bound_is_computed()

    minimum = RF.min([(h.data.bound, h.data.interval) for h in hs], Interval(frac(1,2), 1))

    # Pack into Hypothesis objects
    best_bound = []
    for (func, interval, ref) in minimum:
        zde = Zero_Density_Estimate.from_rational_func(func, interval)
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
    
    if verbose:
        for b in best_bound:
            if b.data.interval.x0 < Constants.ZERO_DENSITY_SIGMA_LIMIT:
                print(f"{b.data}  {b.proof}")

        # plot zero-density estimate
        N = 500
        xs = []
        computed_zdt = []
        literature_zdt = []

        for i in range(N):
            sigma = 1 / 2 + 1 / 2 * i / N
            xs.append(sigma)
            min_ = 1000000
            for h in hs:
                if h.data.interval.contains(sigma):
                    q = h.data.at(sigma)
                    if q.is_real and math.isfinite(q) and min_ > q:
                        min_ = q
            A1 = min_
            A2 = next((b.data.at(sigma) for b in best_bound if b.data.interval.contains(sigma)), 0)
            literature_zdt.append(A1)
            computed_zdt.append(A2)

        plt.figure(dpi=1200)
        plt.xlabel(r"$\sigma$")
        plt.ylabel(r"$A(\sigma)$")
        plt.plot(xs, computed_zdt, linewidth=0.5, label="Computed zero-density estimate")
        plt.plot(xs, literature_zdt, linewidth=0.5, label="Literature zero-density estimate")
        plt.title("Best zero density estimate")
        plt.legend(loc="lower left")
        plt.show()

    return best_bound

def optimize_pintz_zero_density(hypotheses):

    # compute best beta bounds
    hypotheses.add_hypotheses(
        ep.compute_exp_pairs(hypotheses, search_depth=5, prune=True)
    )
    hypotheses.add_hypotheses(ep.exponent_pairs_to_beta_bounds(hypotheses))
    beta_hyps = ep.compute_best_beta_bounds(hypotheses)

    for b in beta_hyps:
        (m, c) = b.data.bound.m, b.data.bound.c
        a_interval = b.data.bound.domain

        # solve t0 = t0(sigma) such that t0 beta (1/t0) = sigma
        # as a linear function of sigma
        # Since t beta (1/t) \leq m + ct, t0 = -m/c + sigma / c
        # Range is t1 <= t <= t2 i.e. m + c t1 <= sigma <= m + c t2
        if c > 0:

            # exclude this edge case for now
            if a_interval.x0 == 0: continue

            sigma_domain = Interval(
                m + c / a_interval.x1, m + c / a_interval.x0,
                a_interval.include_upper, a_interval.include_lower)
            tau0 = Affine(m/c, -1/c, sigma_domain)
            print('tau_0 ', tau0)
            print(m, c, a_interval)
        elif c == 0:
            tau0 = None
            print(m, a_interval)
        else:
            raise NotImplementedError()
