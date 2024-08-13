# code for bounding the function A(\sigma) in the bound
#
# N(\sigma, T) \ll T^{A(\sigma)(1 - \sigma) + o(1)}
#
# where N(\sigma, T) is the number of zeroes of the Riemann zeta-function
# in the rectangle \sigma \le \Re z \le 1 and |\Im z| \le T.

from constants import Constants, Proof_Optimization_Method
import bound_beta as bb
import exponent_pair as ep
from fractions import Fraction as frac
from functions import Affine2, Interval, Piecewise, Polytope, RationalFunction as RF
from hypotheses import *
import large_values as lv
import matplotlib.pyplot as plt
from reference import Reference
import sympy
import scipy.optimize
import zeta_large_values as zlv


import time


###############################################################################
# Object representing an estimate on A(\sigma) represented as a piecewise affine
# function. 
#
# The life cycle of this object is:
# 1) Creation with expr (string), interval (Interval)
# 2) As needed, lazily parse expr to create bound (RationalFunction) that handles 
# evaluations
# 3) Once bound is set, expr is only used for stringify methods
#
# When creating derived zero-density estimates, occasionally it is more convenient
# to initialise with the bound object directly. In that case, use the static 
# from_rational_func function to create the an instance which sets expr to the 
# default __str__ representation of RationalFunction, which may or may not coincide
# with the original expr used to generate the bound object in the first place. 
# This should cause any circular reference problems (only potential display 
# inconsistencies) since, once the bound object is initialised, the expr object 
# is not used for computation. 
# 
# In the future, we will probably move to a canonical model with a guaranteed 
# one-to-one correspondence between expr and bound. At present this is challenging 
# since for each rational function there are multiple possible valid representations. 
class Zero_Density_Estimate:

    # parameters:
    #   - expr: an expression that represents a function of x
    #   - interval: the Interval object representing the range of validity of the
    #               bound.
    def __init__(self, expr, interval):
        if not isinstance(expr, str):
            raise ValueError("Parameter expr must be of type string")
        if not isinstance(interval, Interval):
            raise ValueError("Parameter interval must be of type Interval")
        self.expr = expr
        self.interval = interval

        # Do not compute this yet
        self.bound = None

    def __repr__(self):
        return f"A(x)(1-x) \leq {self.expr} on {self.interval}"

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


# Adds a zero-density estimate to the hypothesis set
# estimate - bound on A(\sigma)
# interval - the domain of \sigma for which the estimate holds
def add_zero_density(hypotheses, estimate, interval, ref, params=""):
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

# Given a list of large value estimates or zeta large value estimates, compute
#
# sup_{t \in [tau_lower, tau_upper]} LV(s, t) / t
#
# for s in sigma_interval, returning the result as a piecewise function.
# Returns: list of (RationalFunction, Interval) tuples.
def compute_sup_LV_on_tau(hypotheses, sigma_interval, tau_lower, tau_upper):

    # Keep track of which hypothesis each piece came from
    lookup = []
    pieces = []
    for h in hypotheses:
        for p in h.data.bound.pieces:
            lookup.append(h)
            pieces.append(p)

    # debugging only - for extending Heath-Brown's estimates
    # fn = Piecewise(pieces)
    # print('debugging:', tau_lower, tau_upper)
    # for f in fn.pieces:
    #     print(f)
    # fn.plot_domain(xlim=(7/8, 1), ylim=(tau_lower, tau_upper), title='Debugging')

    # Critical points are partition the interval sigma_interval into subintervals
    # s_i, with the property that
    # d/dt(f(s, t)/t) is one-signed for all s \in s_i and (s, t) \in D,
    # where D is a polytope (in this case it is the domain of a piecewise defined function).
    crits = set([sigma_interval.x0, sigma_interval.x1])

    # Also, cache all the faces of the polytopes, which are 4-tuples of
    # (domain, function, interval, estimate), where
    #
    # 1. interval: the sigma interval on which the facet is defined
    # 2. domain: the equation of the facet as a function of sigma (i.e. constraint)
    # 3. function: the value of the function on the domain (i.e. a)
    # 4. hypothesis: the hypothesis object from which the bound is derived
    faces = []
    for pi in range(len(pieces)):
        p = pieces[pi]
        # Compute vertices and faces of the domain
        verts = p.domain.get_vertices()
        incidences = [list(x) for x in p.domain.polyhedron.get_input_incidence()]
        constraints = p.domain.get_constraints()

        crits.update([v[0] for v in verts if sigma_interval.contains(v[0])])
        # d/dt(f(s,t)/t) = -(A + Bs)/t^2, which vanishes for s = -A/B
        if p.a[1] != 0:
            s = -p.a[0] / p.a[1]
            if sigma_interval.contains(s):
                crits.add(s)

        # Iterate through the edges
        for i in range(len(constraints)):
            c = constraints[i].coefficients
            vs = [verts[j] for j in incidences[i]]
            sub = Interval(min(v[0] for v in vs), max(v[0] for v in vs), True, True)

            if sub.length() > 0:
                faces.append((RF([-c[1], -c[0]], [c[2]]), p, sub, lookup[pi]))

    return max_RF(crits, faces)

# Bounds are 4-tuples of
# (domain, function, interval, estimate), where
#
# 1. interval: the sigma interval on which the facet is defined
# 2. domain: the equation of the facet as a function of sigma (i.e. constraint)
# 3. function: the value of the function on the domain (i.e. a)
# 4. hypothesis: the hypothesis object from which the bound is derived
def max_RF(crits, faces):

    # Iterate through each subinterval of sigma
    crits = list(crits)
    crits.sort()

    soln = []
    for i in range(1, len(crits)):
        s1 = crits[i - 1]
        s2 = crits[i]
        interval = Interval(s1, s2)
        if interval.length() == 0:
            continue

        s = interval.midpoint()  # the test point

        # Iterate through the faces, collecting t-coordinates where the vertical
        # line \sigma = s intersects with a face
        sup = float("-inf")
        argmax = None
        intersecting_faces = [f for f in faces if f[2].contains(s)]
        dependencies = {}
        for constraint, func, _, hyp in intersecting_faces:
            t = constraint.at(s)
            q = func.at([s, t]) / t
            if q > sup:
                sup = q
                # objective function (a[0] + a[1]s + a[2]t)/t when t = f(s)
                argmax = (
                    RF([func.a[1], func.a[0]]).div(constraint).add(func.a[2]),
                    hyp,
                )
            dependencies[str(hyp.data)] = hyp

        soln.append((argmax[0], interval, list(dependencies.values())))

    # Simplify ranges by merging neighbouring intervals
    for i in range(len(soln) - 1, 0, -1):
        s1 = soln[i - 1]
        s2 = soln[i]
        if s1[0] == s2[0] and s1[1].x1 == s2[1].x0:
            # Merge
            s1[1].x1 = s2[1].x1
            soln.pop(i)
    return soln


# Computes the maximum of
#
# \sup_{t \in [\tau_0, 2\tau_0]} LV(s, t) / t,
# \sup_{t \in [2, \tau_0]} LV_{\zeta}(s, t) / t
#
# for all s \in sigma_interval, as a piecewise RationalFunction. This function
# should return the same result as approximate_best_zero_density_estimate.
def lv_zlv_to_zd(hypotheses, sigma_interval, tau0=frac(3), debug=False):

    s_lim = (sigma_interval.x0, sigma_interval.x1)
    
    if debug:
        start_time = time.time()
        
    # Get large value bounds
    hyps = lv.best_large_value_estimate(
        hypotheses, Polytope.rect(s_lim, (tau0, 2 * tau0))
    )
    
    if debug:
        print(time.time() - start_time, "s")
        start_time = time.time()
        print(f"computing sup LV with {len(hyps)} estimates")

    sup1 = compute_sup_LV_on_tau(
        hyps, sigma_interval, tau0, 2 * tau0
    )

    if debug:
        print(time.time() - start_time, "s")
        start_time = time.time()

    # Get zeta large value bounds
    hyps = zlv.best_large_value_estimate(
        hypotheses, Polytope.rect(s_lim, (frac(2), tau0))
    )
    
    if debug:
        print(time.time() - start_time, "s")
        start_time = time.time()
        print(f"computing sup LVZ with {len(hyps)} estimates")

    sup2 = compute_sup_LV_on_tau(
        hyps, sigma_interval, frac(2), tau0
    )
    
    if debug:
        print(time.time() - start_time, "s")

    # Compute the maximum as a piecewise function
    crits = set(s[1].x0 for s in sup1)
    crits.update(s[1].x1 for s in sup1)
    crits.update(s[1].x0 for s in sup2)
    crits.update(s[1].x1 for s in sup2)

    for func1, int1, _ in sup1:
        for func2, int2, _ in sup2:
            crits.update(func1.intersections(func2, int1.intersect(int2)))

    crits = list(crits)
    crits.sort()
    soln = []
    for i in range(1, len(crits)):
        s1 = crits[i - 1]
        s2 = crits[i]
        interval = Interval(s1, s2)
        if interval.length() == 0:
            continue

        s = interval.midpoint()  # the test point

        f1 = next(f for f in sup1 if f[1].contains(s))
        f2 = next(f for f in sup2 if f[1].contains(s))
        f = f1[0] if f1[0].at(s) > f2[0].at(s) else f2[0]
        soln.append((f, interval, [f1[2], f2[2]]))

    # Simplify ranges by merging neighbouring intervals
    for i in range(len(soln) - 1, 0, -1):
        s1 = soln[i - 1]
        s2 = soln[i]
        if s1[0] == s2[0] and s1[1].x1 == s2[1].x0:
            # Merge
            s1[1].x1 = s2[1].x1
            soln.pop(i)
    
    # pack into Hypothesis
    hyps = []
    for s in soln:
        proof = f'Follows from {len(s[2][0])} large value estimates and {len(s[2][1])} zeta large value estimates'
        deps = s[2][0]          # the LV dependencies
        deps.extend(s[2][1])    # the LVZ dependencies
        hyps.append(derived_zero_density_estimate(
            Zero_Density_Estimate.from_rational_func(s[0], s[1]), 
            proof, 
            deps
        ))
    return hyps

# Tries to prove the zero-density estimate 
# A(sigma) \leq Abound (sigma in sigma_interval)
# using Corollary 11.8
def prove_density_estimate(hypothesis, Abound, sigma_interval):

    # Try to prove that LV(s, t) / t \leq 3(1 - s)/tau0 = (1 - s) * Abound in the range 
    # 2/3 tau0 \leq t \leq tau0
    lv.prove_LV_on_tau_bound(hypothesis, Abound.mul(RF([-1, 1])), sigma_interval, (RF([2]), Abound))
    zlv.prove_LV_on_tau_bound(hypothesis, Abound.mul(RF([-1, 1])), sigma_interval, (RF([2]), Abound))

    pass


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
    print(zde, dep)
    return derived_zero_density_estimate(
        zde, f"Follows from the exponent pair {dep.data}", {dep}
    )


# Computes the zero-density estimate
#
# A(s) \leq 4k/(2(1 + k)s - 1 - l)   (s > s0)
#
# for any exponent pair (k, l) and number s0 jointly satisfying
#   - k < 11/85
#   - 11/85 < k < 1/5
#   - s0 = (144k - 11l - 11)/(170k - 22)
def bourgain_ep_to_zd(exp_pairs):

    bounds = []
    for eph in exp_pairs:
        (k, l) = eph.data.k, eph.data.l
        if k <= frac(1, 5) and l >= frac(3, 5) and 15 * l + 20 * k >= 13:
            if k <= frac(11, 85):
                s0 = max(frac(1, 2), (l + 1) / (2 * (k + 1)))
                if s0 >= 1:
                    continue
                bounds.append(
                    (RF([4 * k], [2 * (1 + k), -1 - l]), Interval(s0, 1), eph)
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
                    (RF([4 * k], [2 * (1 + k), -1 - l]), Interval(s0, 1), eph)
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

    for s in soln:
        print(s[0], "for x \\in", s[1], s[1].contains(frac(15, 16)), s[2])
        s[2].recursively_list_proofs()
    return soln


# Compute all zero-density estimates derived from exponent pairs
def ep_to_zd(hypotheses):

    # TODO: this routine is used quite often, should we roll it into a separate method?
    hypotheses.add_hypotheses(
        ep.compute_exp_pairs(hypotheses, search_depth=5, prune=True)
    )
    hypotheses.add_hypotheses(ep.exponent_pairs_to_beta_bounds(hypotheses))
    hypotheses.add_hypotheses(ep.compute_best_beta_bounds(hypotheses))
    ephs = ep.beta_bounds_to_exponent_pairs(hypotheses)

    zdts = [bourgain_ep_to_zd(ephs)]
    zdts.append(ivic_ep_to_zd(ephs, m=2))

    return zdts

def approx_optimise_bourgain_zero_density_estimate():
    
    sigma_range = (frac(25,32), 1)
    tau_range = (1, 3)

    TOL = 1e-6
    N = 1000
    max_A = 0

    regions = {}
    for i in range(N + 1):
        # for each sigma, compute the tau limits 
        sigma = sigma_range[0] + (sigma_range[1] - sigma_range[0]) * i / N
        #tau_lower = min(frac(3,2), (24 * sigma - 18) / (2 * sigma - 1))
        tau_lower = min(frac(3,2), 48 * (8 * sigma ** 2 - 10 * sigma + 3) / (25 * sigma - 17))

        for j in range(N + 1):
            tau = tau_range[0] + (tau_range[1] - tau_range[0]) * j / N

            # In this region, we cannot guarantee that rho \leq 1 using Jutila's k = 3 
            # estimate - for now just ignore the region
            if max(2 - 2 * sigma, tau + 18 - 24 * sigma) > min(1, 4 - 2 * tau): 
                continue
            
            # Treating sigma, tau as fixed, solve the LP
            # min (a1, a2)
            #       max {f1(a1,a2), f2(a1,a2), ... , f5(a1,a2)} 
            # s.t. 
            #       a1 >= 0, 
            #       a2 >= 0
            # 
            # which is equivalent to the LP 
            # 
            # min (a1, a2, M) 
            #       M
            # s.t. 
            #       M >= f1(a1,a2), ..., M >= f5(a1,a2),
            #       a1 >= 0, 
            #       a2 >= 0

            # objective function is u(a1, a2, M) := M
            obj_func = [0, 0, 1]
            # Constraints of the form Ax \leq b
            Ab = [
                ([0, 1, -1], -2 + 2 * sigma),               # f1(s, t) = a2 + (2 - 2s) < M
                ([1, frac(1,2), -1], -2 + 2 * sigma),       # f2(s, t) = a1 + a2/2 + (2 - 2s) < M
                ([0, -1, -1], -(2 * tau + 4 - 8 * sigma)),  # f3(s, t) = -a2 + (2t + 4 - 8s) < M
                ([-2, 0, -1], -(tau + 12 - 16 * sigma)),    # f4(s, t) = -2a1 + (t + 12 - 16s) < M
                ([4, 0, -1], -(2 + max(1, 2 * tau - 2) - 4 * sigma)),# f5(s, t) = 4a1 + (2 + max(1, 2t - 2) - 4s) < M
                ([-1, 0, 0], 0),                             # a1 >= 0
                ([0, -1, 0], 0)                             # a2 >= 0
            ]
            A = [r[0] for r in Ab]
            b = [r[1] for r in Ab] 

            res = scipy.optimize.linprog(obj_func, A_ub=A, b_ub=b)
            if not res.success:
                print('Warning: linprog did not succeed')

            # check against Bourgain's choice of (sigma, tau)
            if tau > 4 * (1 + sigma) / 5:
                bourgain_alpha1 = tau / 8 - (9 * sigma - 7) / 2
                bourgain_alpha2 = 5 * tau / 4 - (1 + sigma)
            else:
                bourgain_alpha1 = tau / 3 - 2 / 3 * (7 * sigma - 5)
                bourgain_alpha2 = 0
            
            # assert abs(res.x[0] - bourgain_alpha1) < TOL and abs(res.x[1] - bourgain_alpha2) < TOL


            # Save for the plot later 
            equations = []
            if res.slack[0] < TOL:
                equations.append("a2 + 2 - 2s = M")
            if res.slack[1] < TOL:
                equations.append("a1 + a2/2 + 2 - 2s = M")
            if res.slack[2] < TOL:
                equations.append("-a2 + 2t + 4 - 8s = M")
            if res.slack[3] < TOL:
                equations.append("-2a1 + t + 12 - 16s = M")
            if res.slack[4] < TOL:
                equations.append("4a1 + 2 + max(1, 2t - 2) - 4s = M")
            if res.slack[5] < TOL:
                equations.append("a1 = 0")
            if res.slack[6] < TOL:
                equations.append("a2 = 0")
            
            key = str(equations)
            if key not in regions:
                regions[key] = []
            regions[key].append((sigma, tau))

    for key in regions:
        if len(regions[key]) > 0:
            xs = [v[0] for v in regions[key]]
            ys = [v[1] for v in regions[key]]
            plt.plot(xs, ys, label=key)
    plt.legend(loc="lower right")
    plt.xlabel("sigma")
    plt.ylabel("tau")
    plt.title("Bourgain's LV(sigma, tau) estimate: equations governing the optimal choices of alpha_1, alpha_2")
    plt.show()


def optimise_bourgain_zero_density_estimate():

    # The domain of definition is (sigma, tau, alpha1, alpha2) jointly satisfying
    # 25/32 <= sigma <= 11/14
    # 1 <= tau <= 3/2
    # alpha1 >= 0
    # alpha2 >= 0

    ALPHA_MAX = 1000
    constraints = [
        [-frac(25, 32), 1, 0, 0, 0],  # \sigma >= 25/32
        [frac(11, 14), -1, 0, 0, 0],  # \sigma <= 11/14
        [-1, 0, 1, 0, 0],  # \tau >= 1
        [frac(3, 2), 0, -1, 0, 0],  # \tau <= 3/2
        [0, 0, 0, 1, 0],  # \alpha_1 >= 0
        [0, 0, 0, 0, 1],  # \alpha_2 >= 0
        [ALPHA_MAX, 0, 0, -1, 0],  # \alpha_1 <= ALPHA_MAX
        [ALPHA_MAX, 0, 0, 0, -1],  # \alpha_2 <= ALPHA_MAX
    ]

    domain = Polytope(constraints)

    # Objective function (bounds on rho)
    bounds = [
        [2, -2, 0, 0, 1],  # \alpha_2 + 2 - 2\sigma
        [2, -2, 0, 1, frac(1, 2)],  # \alpha_1 + \alpha_2 / 2 + 2 - 2\sigma
        [4, -8, 2, 0, -1],  # -\alpha_2 + 2\tau + 4 - 8\sigma
        [12, -16, 1, 2, 0],  # 12 - 16\sigma + \tau + 2\alpha_1
        [3, -4, 0, 4, 0],  # 3 - 4\sigma + 4\alpha_1
    ]

    # Construct affine objects to represent the bounds
    fns = [Affine2(b, domain) for b in bounds]

    # Compute the pairwise intersections
    intersections = []
    for i in range(len(bounds)):
        for j in range(0, i):
            intersections.append(fns[i].intersection(fns[j]))

    # Iterate over the power set of intersection lines
    regions = []
    for i in range(2 ** len(intersections)):
        b = bin(i)[2:].zfill(len(intersections))  # binary representation
        s = list(constraints)  # shallow copy
        for j in range(len(b)):
            if b[j] == "1":
                s.append(intersections[j].constraint.coefficients)
            else:
                s.append([-x for x in intersections[j].constraint.coefficients])
        p = Polytope(s, linear=False, canonicalize=True)
        if not p.is_empty(include_boundary=False):
            regions.append(p)

    # For each region, compute the centroid and find the maximal function
    pieces = []
    for r in regions:
        center = r.get_centroid()
        index = max(range(len(bounds)), key=lambda i: fns[i].at(center))
        argmax = fns[index]
        pieces.append(Affine2(argmax.a, r))

    # Since the domains are in 4-dimensional space, facets are 3-dimensional and
    # ridges are 2-dimensional. Every ridge is the intersection of two facets
    # however since the domain is not necessarily a simplex, not every combination
    # of two facets is a ridge. We need an algorithm for computing all ridges
    # of a 4-polytope.
    """
    p = pieces[0]
    verts = p.domain.get_vertices()
    incidences = [list(x) for x in p.domain.polyhedron.get_input_incidence()]
    constraints = p.domain.get_constraints()

    for i in range(len(constraints)):
        vs = [verts[j] for j in incidences[i]]
        print('Vertices of the facet', constraints[i], 'are:')
        for v in vs:
            print('\t', [str(x) for x in v])
    """
    # test

    N = 30
    for i in range(N):
        sigma = frac(25, 32) + (frac(11, 14) - frac(25, 32)) * i / N
        for j in range(N):
            tau = 1 + (frac(3, 2) - 1) * j / N

            if tau < (24 * sigma - 18) / (2 * sigma - 1):
                continue

            min_bound = 1000000000
            argmin = None
            for k in range(N):
                alpha1 = frac(k, 2 * N)  # assume alpha1 is in [0, 1/2]

                for l in range(N):
                    alpha2 = frac(l, 2 * N)  # assume alpha2 is in [0, 1/2]

                    f = max(fn.at([sigma, tau, alpha1, alpha2]) for fn in fns)

                    if f < min_bound:
                        min_bound = f
                        argmin = [alpha1, alpha2, f]

            # Bourgain's choice - this was to check the max argument
            bp = []
            if tau > 4 / 5 * (1 + sigma):
                bp = [tau / 8 - (9 * sigma - 7) / 2, 5 * tau / 4 - (1 + sigma)]
            else:
                bp = [tau / 3 - 2 / 3 * (7 * sigma - 5), 0]

            bf = max(fn.at([sigma, tau, bp[0], bp[1]]) for fn in fns)
            print(
                sigma,
                tau,
                tau > 4 / 5 * (1 + sigma),
                argmin,
                float(min_bound),
                bf,
                float(bf) / (1 - sigma) / tau,
            )
            if bf < min_bound:
                print("our choices:", argmin[0:2], "gives", min_bound)
                print([float(fn.at([sigma, tau, argmin[0], argmin[1]])) for fn in fns])
                print("bourgain choices", bp, "gives", bf)
                print([fn.at([sigma, tau, bp[0], bp[1]]) for fn in fns])

                raise ValueError()
