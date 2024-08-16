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
from functions import Affine2, Interval, Piecewise, Polytope, RationalFunction as RF, SympyHelper
from hypotheses import *
import large_values as lv
import matplotlib.pyplot as plt
import numpy as np
from reference import Reference
import sympy
import scipy
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
    # fn.plot_domain(xlim=(0.75, 0.82), ylim=(tau_lower, tau_upper), title='Debugging')

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

    # print('sup 1')
    # for p in sup1:
    #     print(p[0], p[1])
    # for p in sup1:
    #     if p[1].contains(0.78):
    #         for q in p[2]:
    #             q.recursively_list_proofs(1)
    
    # temp = []
    # for p in sup1:
    #     if p[1].contains(0.78):
    #         for q in p[2]:
    #             temp.extend(q.data.bound.pieces)
    # temp = Piecewise(temp)
    # temp.plot_domain(xlim=(0.75, 0.82), ylim=(tau0, 2 * tau0), resolution=1000)
    
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
    
    # print('sup 2')
    # for p in sup2:
    #     print(p[0], p[1])
    # for p in sup2:
    #     if p[1].contains(0.78):
    #         for q in p[2]:
    #             q.recursively_list_proofs(1)
    # temp = []
    # for p in sup2:
    #     if p[1].contains(0.78):
    #         for q in p[2]:
    #             temp.extend(q.data.bound.pieces)
    # temp = Piecewise(temp)
    # temp.plot_domain(xlim=(0.75, 0.82), ylim=(2, tau0), resolution=1000)
    
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

    raise NotImplementedError()


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

    # zdts = [bourgain_ep_to_zd(ephs)]
    approx_bourgain_ep_to_zd(ephs)
    # zdts.append(ivic_ep_to_zd(ephs, m=2))

    return zdts


# Aggregate the zero-density estimates in the Hypothesis_Set and returns a piecewise 
# function that represents the best zero-density estimate in each subinterval of [1/2, 1]
# Note that this function does not compute zero-density estimates from other 
# Hypothesis objects - it only aggregates existing Hypothesis objects of type 
# "Zero density estimate". 
def best_zero_density_estimate(hypotheses, verbose=False):
    hs = hypotheses.list_hypotheses(hypothesis_type="Zero density estimate")

    # Ensure bound is computed
    for h in hs:
        h.data._ensure_bound_is_computed()
    
    # Start with default bound 
    default = Zero_Density_Estimate("10000000", Interval(frac(1,2), 1))
    default._ensure_bound_is_computed()
    best_bound = [
        Hypothesis(
            "Placeholder zero density estimate",
            "Zero density estimate",
            default,
            "Placeholder zero density estimate",
            Reference.trivial(),
        )
    ]
    x = RF.x

    for h1 in hs:
        # For simplicity, work directly with sympy objects
        f1 = h1.data.bound.num / h1.data.bound.den
        in1 = h1.data.interval 

        new_best_bound = []
        for h2 in best_bound:
            f2 = h2.data.bound.num / h2.data.bound.den
            in2 = h2.data.interval
            solns = sympy.solve(f1 - f2)
            crits = set(SympyHelper.to_frac(soln) for soln in solns if soln.is_real)
            crits.update([in1.x0, in1.x1])
            crits = set(c for c in crits if in2.contains(c))
            crits.update([in2.x0, in2.x1])

            crits = list(crits)
            crits.sort()
            for i in range(1, len(crits)):
                inter = Interval(crits[i - 1], crits[i])
                mid = inter.midpoint()
                if not in1.contains(mid) or f2.subs(x, mid) < f1.subs(x, mid):
                    # f2 is the better bound 
                    zde = Zero_Density_Estimate(str(f2), inter)
                    h = h2
                else:
                    zde = Zero_Density_Estimate(str(f1), inter)
                    h = h1

                zde._ensure_bound_is_computed()
                # name, hypothesis_type, data, proof, reference
                new_best_bound.append(
                    Hypothesis(h.name, h.hypothesis_type, zde, h.proof, h.reference)
                )

        # Simplify
        best_bound = []
        i = 0
        while i < len(new_best_bound):
            bi = new_best_bound[i]
            (fi, inti) = (bi.data.bound, bi.data.interval)
            left = inti.x0
            right = inti.x1

            j = i + 1
            while j < len(new_best_bound):
                bj = new_best_bound[j]
                (fj, intj) = (bj.data.bound, bj.data.interval)
                if not (fi == fj and right == intj.x0 and \
                        bi.proof == bj.proof):
                    break

                right = intj.x1
                j += 1

            bi.data.interval = Interval(left, right)
            best_bound.append(bi)
            i = j

    if verbose:
        print("A(x) \leq ")
        for b in best_bound:
            if b.data.interval.x0 < Constants.ZERO_DENSITY_SIGMA_LIMIT:
                print(f"\t{b.data}  {b.proof}")

        # plot zero-density estimate 
        N = 500
        xs = []
        computed_zdt = []
        literature_zdt = []

        for i in range(N):
            sigma = 1 / 2 + 1 / 2 * i / N
            xs.append(sigma)
            A1 = min(h.data.at(sigma) for h in hs if h.data.interval.contains(sigma))
            A2 = next((b.data.at(sigma) for b in best_bound if b.data.interval.contains(sigma)), 0)
            literature_zdt.append(A1)
            computed_zdt.append(A2)

        plt.figure(dpi=1200)
        plt.xlabel("σ")
        plt.ylabel("A(σ)")
        plt.plot(xs, computed_zdt, linewidth=0.5, label="Computed zero-density estimate")
        plt.plot(xs, literature_zdt, linewidth=0.5, label="Literature zero-density estimate")
        plt.title("Best zero density estimate")
        plt.legend(loc="lower left")
        plt.show()
    
    return best_bound

