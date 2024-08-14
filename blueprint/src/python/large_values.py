# code for representing a an upper bound on how frequently a Dirichlet polynomial
# of length N can be large. For

from constants import *
import copy
from fractions import Fraction as frac
from functions import *
from hypotheses import *
import numpy as np
import matplotlib.pyplot as plt
from polytope import *
from reference import *
import itertools
import time


###############################################################################
# Object representing a large value bound for a certain domain on (\sigma, \tau)
class Large_Value_Estimate:

    # parameters:
    #   - bound: (Piecewise type) the piecewise affine function (of 2 variables)
    #           representing the bound on LV(\sigma, \tau)
    def __init__(self, bound):
        if not isinstance(bound, Piecewise):
            raise ValueError("bound must be of type Piecewise")
        self.bound = bound

    def __repr__(self):
        return "LV(x, y) \leq " + str(self.bound)


class Large_Value_Estimate_Transform:

    def __init__(self, transform):
        self.transform = transform

    def __repr__(self):
        return "Raising to a power"


###############################################################################


# Compute the maximum of a list of bounds represented in coefficient vector form
# Returns the result as a piecewise function
# bounds (list of list) vectors a representing a function f(x) = a^Tx with the first
#                   coefficient denoting the constant term
# domain (Polytope object) the (sigma, tau) domain of definition
def max_of(bounds, domain=None):

    # The standard domain of definition is \sigma \in [1/2, 1], \tau \geq 0
    if domain is None:
        domain_ineq = [
            [-frac(1, 2), 1, 0],  # \sigma >= 1/2
            [1, -1, 0],  # \sigma <= 1
            [0, 0, 1],  # \tau >= 0
            [Constants.TAU_UPPER_LIMIT, 0, -1],  # \tau <= large number
        ]
        domain = Polytope(domain_ineq)

    # Construct affine objects to represent the bounds
    fns = [Affine2(b, domain) for b in bounds]

    # Handle edge cases
    if len(fns) == 1:
        return Piecewise(fns)

    # Compute the pairwise intersections
    intersections = []
    for i in range(len(bounds)):
        for j in range(0, i):
            intersections.append(fns[i].intersection(fns[j]))

    # Iterate over the power set of intersection lines
    regions = []
    for i in range(2 ** len(intersections)):
        b = bin(i)[2:].zfill(len(intersections))  # binary representation
        s = []
        for j in range(len(b)):
            if b[j] == "1":
                s.append(intersections[j].constraint.coefficients)
            else:
                s.append([-x for x in intersections[j].constraint.coefficients])
        p = Polytope(s).intersect(domain)
        if not p.is_empty(include_boundary=False):
            regions.append(p)

    # For each region, compute the centroid and find the maximal function
    pieces = []
    for r in regions:
        center = r.get_centroid()
        index = max(range(len(bounds)), key=lambda i: fns[i].at(center))
        argmax = fns[index]
        pieces.append(Affine2(argmax.a, r))

    # Construct the piecewise function and simplify
    bound = Piecewise(pieces)
    bound.simplify()

    return bound


def literature_bound_LV_max(bounds, ref, params=""):
    piecewise = max_of(bounds)
    return Hypothesis(
        f"{ref.author()} large value estimate" + params,
        "Large value estimate",
        Large_Value_Estimate(piecewise),
        f"See [{ref.author()}, {ref.year()}]",
        ref,
    )


def derived_bound_LV(bound, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived large value estimate",
        "Large value estimate",
        Large_Value_Estimate(bound),
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound


def classical_LV_estimate(bounds):
    return Hypothesis(
        "Classical large value estimate",
        "Large value estimate",
        Large_Value_Estimate(max_of(bounds)),
        "Classical",
        Reference.classical(),
    )


def conjectured_LV_estimate(bounds, name):
    return Hypothesis(
        name,
        "Large value estimate",
        max_of(bounds),
        f"Conjecture",
        Reference.conjectured(),
    )

montgomery_conjecture = conjectured_LV_estimate([[2, -2, 0]], "Montgomery conjecture")

def get_optimized_bourgain_lv_estimate(ref):
    pieces = [
        # For now - assume that we can't say anthing about LV estimates
        # for tau < 1
        Affine2(
            [Constants.LV_DEFAULT_UPPER_BOUND, 0, 0],
            Polytope([
                [-frac(1,2), 1, 0],  # \sigma >= 1/2
                [1, -1, 0],  # \sigma <= 1
                [0, 0, 1],  # \tau >= 0
                [1, 0, -1],  # \tau <= 1
            ])
        ),
        Affine2(
            [frac(16,3), -frac(20,3), frac(1,3)],
            Polytope([
                [10, -14, 1], # 10 - 14s + t >= 0
                [-1, 0, 1], # -1 + t >= 0
                [4, 4, -5], # 4/5 + 4/5s - t >= 0
                [-11, 16, -1]
            ])
        ),
        Affine2(
            [5, -7, frac(3,4)],
            Polytope([
                [8, -8, -1],
                [-16, 20, frac(1,3)],
                [-6, 10, -frac(7,6)],
                [-4, -4, 5]
            ])
        ),
        Affine2(
            [3, -5, 1],
            Polytope([
                [-8, 8, 1],
                [2, -6, 2],
                [-10, 14, -frac(2,3)],
                [6, -2, -2]
            ])
        ),
        Affine2(
            [0, -4, 2],
            Polytope([
                [-6, 2, 2],
                [-12, 12, 1],
                [Constants.TAU_UPPER_LIMIT, 0, -1],
                [1, -1, 0]
            ])
        ),
        Affine2(
            [8, -12, frac(4,3)],
            Polytope([
                [15, -21, 1],
                [12, -12, -1],
                [-frac(3,2), 0, 1],
                [6, -10, frac(7,6)]
            ])
        ),
        Affine2(
            [2, -2, 0],
            Polytope([
                [1, -1, 0],
                [-10, 14, -1],
                [-1, 0, 1],
                [-2, 6, -2]
            ])
        ),
        Affine2(
            [9, -12, frac(2,3)],
            Polytope([
                [frac(3,2), 0, -1],
                [-frac(1,2), 1, 0],
                [-1, 0, 1],
                [11, -16, 1],
                [16, -20, -frac(1,3)]
            ])
        )
    ]
    return Hypothesis(
        f"{ref.author()} optimized large value estimate",
        "Large value estimate",
        Large_Value_Estimate(Piecewise(pieces)),
        f"See [{ref.author()}, {ref.year()}]",
        ref,
    )




###############################################################################

# Classical estimates

# L^2 mean value theorem: LV(s, t) \leq max(2 - 2s, 1 - 2s + t)
large_value_estimate_L2 = classical_LV_estimate([[2, -2, 0], [1, -2, 1]])


# Raising to a power (large value estimate transform theorem)
def raise_to_power_hypothesis(k):

    # Convert into fraction
    k = frac(k)

    # LV(s, kt) \leq kLV(s, t)
    def transform(hypothesis):
        # Scale while re-imposing upper limit on tau
        new_bound = hypothesis.data.bound.scale(
            2, k, [[Constants.TAU_UPPER_LIMIT, 0, -1]]
        )
        for piece in new_bound.pieces:
            piece.a = [k * ai for ai in piece.a]
            piece.a[2] /= k
        return derived_bound_LV(
            new_bound,
            f"Follows from raising {hypothesis} to the k = {k} power",
            {hypothesis},
        )

    return Hypothesis(
        f"Raising to a power large value estimate transform with k = {k}",
        "Large value estimate transform",
        Large_Value_Estimate_Transform(transform),
        "Classical",
        Reference.classical(),
    )


###############################################################################

# Debugging functions
def covers(estimate, xlim, ylim):
    N = 100
    for xi in range(1, N):
        for yi in range(1, N):
            x = xlim[0] + (xlim[1] - xlim[0]) * xi / N
            y = ylim[0] + (ylim[1] - ylim[0]) * yi / N
            if estimate.at([x, y]) is None:
                for p in estimate.pieces:
                    print(p)
                raise ValueError(f"{x}, {y}")


###############################################################################

# Given a list of Piecewise objects, compute their minimum over a given domain
# Returns result as a list of hypotheses, created using the constructor function
def piecewise_min(estimates, domain, constructor):

    # Compute bounds and crop domains (taking care not the alter the original estimates
    # objects)
    bounds = [e.data.bound for e in estimates]
    for i in range(len(bounds)):
        bounds[i] = copy.copy(bounds[i])
        bounds[i].crop(domain)

    # No combination required - just crop a copy of the pieces to the required domain
    if len(estimates) == 1:
        parent = estimates[0]
        return [constructor(bounds[0], f"Follows from {parent.name}", {parent})]

    # Temporarily set the 'label' field of all estimates, and create a lookup table
    lookup = {}
    for i in range(len(bounds)):
        b = bounds[i]
        lookup[i] = estimates[i]
        for piece in b.pieces:
            piece.label = i

    # Iterate through the list of estimates, taking pairwise minimum each time
    best_est = None
    simplify_every = 5
    for i in range(len(bounds)):
        b = bounds[i]
        if best_est is None:
            best_est = b
        else:
            best_est = best_est.min_with(b)
        if i % simplify_every == 0:
            best_est.simplify()

    if len(bounds) % simplify_every != 0:
        best_est.simplify()

    # Set hypothesis objects with dependencies and proofs
    hyps = []
    if best_est is not None:
        for e in best_est.pieces:
            parent = lookup[e.label]
            hyps.append(
                constructor(Piecewise([e]), f"Follows from {parent.name}", {parent})
            )

    # Remove 'label' field
    for b in bounds:
        for piece in b.pieces:
            piece.label = None

    return hyps


# Returns the best estimate on LV(\sigma, \tau) by combining all hypotheses of
# type 'Large value estimate' in the hypothesis set.
# Warning: this method is not thread safe
# Parmeters:
#   - hypotheses: (Hypothesis_Set object) the set of hypothesis to assume
#   - domain: (Polytope object) [Optional]the domain on which to compute the large value estimate
#               must be of dimension 2.
# Returns: (list of Hypothesis objects)
def best_large_value_estimate(hypotheses, domain=None):
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("hypotheses must be of type Hypothesis_Set")
    if domain is not None and not isinstance(domain, Polytope):
        raise ValueError("domain must be of type Polytope")

    # Default domain (\sigma, \tau) \in [1/2, 1] x [0, TAU_UPPER_LIMIT]
    if domain is None:
        domain = Polytope.rect((frac(1, 2), frac(1)), (0, Constants.TAU_UPPER_LIMIT))

    lves = hypotheses.list_hypotheses(hypothesis_type="Large value estimate")
    lv_transforms = hypotheses.list_hypotheses(
        hypothesis_type="Large value estimate transform"
    )

    # Generate set of LV estimates (original + transformed)
    lv_estimates = list(lves)
    for tr_hyp in lv_transforms:
        lv_estimates.extend([tr_hyp.data.transform(lve) for lve in lves])

    return piecewise_min(lv_estimates, domain, derived_bound_LV)


# Optimise Bourgain's large value estimate by choosing the best value of \alpha_1, \alpha_2
# in each subregion of (\sigma, \tau)
def optimize_bourgain_large_value_estimate():
    # Variables are (in order)
    # [a1, a2, sigma, tau, M, constant]
    # Regions over (sigma, tau) are defined by solving a system of 3 equations
    # to obtain a1 = f(sigma, tau), a2 = g(sigma, tau) where f, g are linear.
    # The equations are given by
    eqns = [
        [0, 1, -2, 0, -1, 2],            # a2 + 2 - 2s = M
        [1, frac(1,2), -2, 0, -1, 2],    # a1 + a2/2 + 2 - 2s = M
        [0, -1, -8, 2, -1, 4],           # -a2 + 2t + 4 - 8s = M
        [-2, 0, -16, 1, -1, 12],         # -2a1 + t + 12 - 16s = M
        [4, 0, -4, 0, -1, 3],            # 4a1 + 3 - 4s = M
        [4, 0, -4, 2, -1, 0],            # 4a1 + 2t - 4s = M
        [1, 0, 0, 0, 0, 0],              # a1 = 0
        [0, 1, 0, 0, 0, 0],              # a2 = 0
    ]

    domain = Polytope.rect((frac(1,2), frac(1)), (frac(1), frac(3)))

    # Sympy solvers give empty solution sets for these systems - so instead we
    # (temporarily) use sympy's matrix rref computer
    # a1, a2, s, t, M = sympy.symbols("a1, a2, s, t, M")
    # var = [a1, a2, s, t, M]
    hypotheses = []
    for c in itertools.combinations(eqns, 3):
        mat = sympy.Matrix([eq for eq in c])
        (rows, cols) = (len(c), len(c[0]))
        res = mat.rref() # Compute reduced row-echelon form
        mat = [[SympyHelper.to_frac(x) for x in res[0].row(i)] for i in range(rows)] # unpack

        # Take advantage of the reduced row-echelon form
        if mat[0][0] == 0 or mat[1][1] == 0:
            print('skipping')
            for r in mat: print(' '.join(str(v) for v in r))

            continue

        # Scale if required (it shouldn't be required, but just in case)
        if mat[0][0] != 1:
            d = mat[0][0]
            mat[0] = [mat[0][i] / d for i in range(cols)]
        if mat[1][1] != 1:
            d = mat[1][1]
            mat[1] = [mat[1][i] / d for i in range(cols)]

        # Eliminate the M variable
        if mat[0][4] != 0:
            if mat[2][4] == 0:
                # Unsolvable for M
                print('skipping 2')
                for r in mat: print(' '.join(str(v) for v in r))
                continue
            else:
                r = mat[0][4] / mat[2][4]
                mat[0] = [mat[0][i] - mat[2][i] * r for i in range(cols)]
        if mat[1][4] != 0:
            if mat[2][4] == 0:
                print('skipping 3')
                for r in mat: print(' '.join(str(v) for v in r))
                continue
            else:
                r = mat[1][4] / mat[2][4]
                mat[1] = [mat[1][i] - mat[2][i] * r for i in range(cols)]

        # Given the definitions of a1, a2, substitute into each of the
        # 6 functions under the max, computing their maximum in the domain
        a1_defn = [-mat[0][5], -mat[0][2], -mat[0][3]] # C + A s + B t
        a2_defn = [-mat[1][5], -mat[1][2], -mat[1][3]] # C + A s + B t

        # This is the region where a1 >= 0, a2 >= 0
        region = Polytope([a1_defn, a2_defn]).intersect(domain)

        if not region.is_empty(include_boundary=False):
            # Hack to ensure uniqueness (using the fact that tuples are hashable)
            fns = set()
            for i in range(6):
                fn = [eqns[i][5], eqns[i][2], eqns[i][3]]
                fn = tuple(fn[j] + eqns[i][0] * a1_defn[j] + eqns[i][1] * a2_defn[j] for j in range(len(fn)))
                fns.add(fn)
            # Compute maximum
            start_time = time.time()
            lst = [list(fn) for fn in fns]
            func = max_of(lst, region)
            neg_regions = domain.set_minus(region)
        else:
            func = Piecewise([])
            neg_regions = [domain]

        for reg in neg_regions:
            func.pieces.append(Affine2([10000000, 0, 0], reg))

        if not func.check((1/2, 1), (1, 3)):
            print("original list")
            for p in lst:
                print(' '.join(str(pi) for pi in p))

            print("maxed")
            for p in func.pieces:
                print(p)

            print("region")
            print(region)

            print("recreation")
            print(max_of(lst, region).check((1/2, 1), (1,3)))



        a1_proof = "a1 = " + Affine2.to_string(a1_defn, "st")
        a2_proof = "a2 = " + Affine2.to_string(a2_defn, "st")
        hypotheses.append(derived_bound_LV(func, f"Follows from taking {a1_proof} and {a2_proof}", {}))

    # Ensure that hypotheses objects are complete 
    print('Checking individual hypotheses')
    for h in hypotheses:
        print('Checking:')
        for p in h.data.bound.pieces:
            print('\t', p)
        h.data.bound.check(xlim=(frac(25,32), 1), ylim=(1, 3))
    
    print('computing piecewise min of ', len(hypotheses))
    best_lv_estimate = piecewise_min(hypotheses, domain, derived_bound_LV)

    pieces = []
    for h in best_lv_estimate:
        for p in h.data.bound.pieces:
            pieces.append(p)

    fn = Piecewise(pieces)
    fn.plot_domain(xlim=(1/2, 1), ylim=(1, 3), title='Before simplifying')

    fn.simplify(5)

    # debugging
    print('-----------------------------------------------------------------------------')
    fn = Piecewise(pieces)
    for f in fn.pieces:
        print(f)
    fn.plot_domain(xlim=(1/2, 1), ylim=(1, 3), title='Debugging')

    return best_lv_estimate


# Tries to prove the bound LV(s, t) / t \leq f(s) on the specified domain defined by
# s \in sigma_range
# t \in tau_range(s)
def prove_LV_on_tau_bound(hypotheses, f, sigma_range, tau_range):
    pass


# Given a large-value estimate as a Hypothesis, apply Huxley subdivison (see Basic
# properties (ii) of Large value estimates section) to obtain a better large
# value estimate.
#
# If the large value estimate is unchanged, returns None. Otherwise, the new
# large value estimate is returned as a Hypothesis
def apply_huxley_subdivision(hypothesis):
    if not isinstance(hypothesis, Hypothesis):
        raise ValueError('Parameter hypothesis must be of type Hypothesis')
    if hypothesis.hypothesis_type != 'Large value estimate':
        raise ValueError('Parameter hypothesis must be a Hypothesis of type "Large value estimate"')

    # Iterate through the pieces, extracting the facets (which are just lines)
    # then compute their projection onto the \sigma space
    for p in lv_hypothesis.data.bound.pieces:
        p = pieces[k]
        # Compute vertices and facets of the domain
        verts = p.domain.get_vertices()
        incidences = [list(x) for x in p.domain.polyhedron.get_input_incidence()]
        constraints = p.domain.get_constraints()

        # Iterate through the edges
        for i in range(len(constraints)):
            c = constraints[i].coefficients
            vs = [verts[j] for j in incidences[i]]
            sub = Interval(min(v[0] for v in vs), max(v[0] for v in vs), True, True)
            if sub.length() > 0:
                faces.append((RF([-c[1], -c[0]], [c[2]]), p, sub, lookup[k]))

    raise NotImplementedError()
