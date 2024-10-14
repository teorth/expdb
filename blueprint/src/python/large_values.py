# code for representing a an upper bound on how frequently a Dirichlet polynomial
# of length N can be large. For

from constants import *
import copy
from fractions import Fraction as frac
from functions import *
from hypotheses import *
import itertools
from polytope import *
from reference import *
from region import Region
import scipy
import time


###############################################################################
"""
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
        return "LV(σ,τ) \leq " + str(self.bound)
"""

# Object representing a large value estimate LV(σ,τ)
class Large_Value_Estimate:
    # Parameters: region (Region) a region in R^3 representing the set of 
    # feasible (sigma, tau, rho) values 
    def __init__(self, region):
        if not isinstance(region, Region):
            raise ValueError("Parameter region must be of type region")
        self.region = region

    def __repr__(self):
        return f"(σ,τ,ρ) in {self.region}"


class Large_Value_Estimate_Transform:

    def __init__(self, transform):
        self.transform = transform

    def __repr__(self):
        return "Raising to a power"

# Object representing a set in R^3 containing feasible (sigma, tau, rho*) values 
class Large_Value_Energy_Estimate:
    
    def __init__(self, region):
        if not isinstance(region, Region):
            raise ValueError("Parameter region must be of type Region")
        self.region = region
        
    def __repr__(self):
        return f"(σ,τ,ρ*) in {self.region}"

###############################################################################
"""
# Compute the maximum of a list of bounds represented in coefficient vector form
# Returns the result as a piecewise function
# bounds (list of list) vectors a representing a function f(x) = a^Tx with the first
#                   coefficient denoting the constant term
# domain (Polytope object) the (sigma, tau) domain of definition
# 
# TODO: since we have implemented the Region class this method may be replaced 
# by e.g. union_of_halfplanes method in additive_energy. 
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
"""

def convert_bounds_to_region(bounds):
    # The default domain of definition in (\sigma, \tau, \rho) space
    domain = Polytope.rect(
        (frac(1,2), frac(1)),                       # 1/2 \le \sigma \le 1
        (frac(0), Constants.TAU_UPPER_LIMIT),       # 0 \le \tau \le TAU_UPPER_LIMIT
        (frac(0), Constants.LV_DEFAULT_UPPER_BOUND),# 0 \le \rho \le RHO_UPPER_LIMIT
    )
    # Convert bounds from 'bounds', of the form 
    # \rho \le f(\sigma, \tau) 
    # into a region in R^3 of the form 
    # f(\sigma, \tau, \rho) \ge 0.
    return Region.from_union_of_halfplanes(
        [b + [-1] for b in bounds],
        domain
    )

def literature_bound_LV_max(bounds, ref, params=""):
    return Hypothesis(
        f"{ref.author()} large value estimate" + params,
        "Large value estimate",
        Large_Value_Estimate(convert_bounds_to_region(bounds)),
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
        Large_Value_Estimate(convert_bounds_to_region(bounds)),
        "Classical",
        Reference.classical(),
    )

def conjectured_LV_estimate(bounds, name):
    return Hypothesis(
        name,
        "Large value estimate",
        Large_Value_Estimate(convert_bounds_to_region(bounds)),
        "Conjecture",
        Reference.conjectured(),
    )

###############################################################################

# Classical estimates

# L^2 mean value theorem: LV(s, t) \leq max(2 - 2s, 1 - 2s + t)
large_value_estimate_L2 = classical_LV_estimate([[2, -2, 0], [1, -2, 1]])

montgomery_conjecture = conjectured_LV_estimate([[2, -2, 0]], "Montgomery conjecture")

###############################################################################

def get_optimized_bourgain_lv_estimate(ref):
    region = Region.as_disjoint_union([
        # For now - assume that we can't say anthing about LV estimates
        # for tau < 1
        Polytope.rect(
            (frac(1,2), frac(1)),                       # 1/2 \le \sigma \le 1
            (frac(0), frac(1)),                         # 0 \le \tau \le 1
            (frac(0), Constants.LV_DEFAULT_UPPER_BOUND) # 0 \le \rho
        ),
        Polytope([
            [frac(16,3), -frac(20,3), frac(1,3), -1],
            [10, -14, 1, 0],            # 10 - 14s + t >= 0
            [-1, 0, 1, 0],              # t >= 1
            [4, 4, -5, 0],              # 4/5 + 4/5s - t >= 0
            [-11, 16, -1, 0]            # t <= 16s - 11
        ]),
        Polytope([
            [5, -7, frac(3,4), -1],     # p <= 5 - 7s + 3s/4
            [8, -8, -1, 0],             # t <= 8 - 8s
            [-16, 20, frac(1,3), 0],    # t >= 3(20s - 16)
            [-6, 10, -frac(7,6), 0],    # -6 + 10s - 7t/6 >= 0
            [-4, -4, 5, 0]              # -4 - 4s + 5t >= 0
        ]),
        Polytope([
            [3, -5, 1, -1],             # p <= 3 - 5s + t
            [-8, 8, 1, 0],
            [2, -6, 2, 0],
            [-10, 14, -frac(2,3), 0],
            [6, -2, -2, 0]
        ]),
        Polytope([
            [0, -4, 2, -1],             # p <= -4s + 2t
            [-6, 2, 2, 0],
            [-12, 12, 1, 0],
            [Constants.TAU_UPPER_LIMIT, 0, -1, 0],
            [1, -1, 0, 0]
        ]),
        Polytope([
            [8, -12, frac(4,3), 0],     # p <= 8 - 12s + 4t/3
            [15, -21, 1, 0],
            [12, -12, -1, 0],
            [-frac(3,2), 0, 1, 0],
            [6, -10, frac(7,6), 0]
        ]),
        Polytope([
            [2, -2, 0, -1],             # p <= 2 - 2s
            [1, -1, 0, 0],
            [-10, 14, -1, 0],
            [-1, 0, 1, 0],
            [-2, 6, -2, 0]
        ]),
        Polytope([
            [9, -12, frac(2,3), -1],    # p <= 9 - 12s + 2t/3
            [frac(3,2), 0, -1, 0],
            [-frac(1,2), 1, 0, 0],
            [-1, 0, 1, 0],
            [11, -16, 1, 0],
            [16, -20, -frac(1,3), 0]
        ])
    ])
    
    return Hypothesis(
        f"{ref.author()} optimized large value estimate",
        "Large value estimate",
        Large_Value_Estimate(region),
        f"See [{ref.author()}, {ref.year()}]",
        ref,
    )


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

    # Compute bounds and crop domains (taking care not the alter the original 
    # estimates objects)
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
            if mat[2][4] == 0: # Unsolvable for M
                continue
            else:
                r = mat[0][4] / mat[2][4]
                mat[0] = [mat[0][i] - mat[2][i] * r for i in range(cols)]
        if mat[1][4] != 0:
            if mat[2][4] == 0:
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

        a1_proof = "a1 = " + Affine2.to_string(a1_defn, "st")
        a2_proof = "a2 = " + Affine2.to_string(a2_defn, "st")
        hypotheses.append(derived_bound_LV(func, f"Follows from taking {a1_proof} and {a2_proof}", {}))

    return piecewise_min(hypotheses, domain, derived_bound_LV)

# Check the literature Hypothesis object against a linear program computing the numerical
# solution for fixed (sigma, tau)
def check_bourgain_large_value_estimate(hypothesis):
    
    sigma_range = (frac(1,2), 1)
    tau_range = (frac("1.001"), 5) # do not include 1

    TOL = 1e-6
    N = 1000
    max_A = 0
    for i in range(N + 1):
        # for each sigma, compute the tau limits 
        sigma = sigma_range[0] + (sigma_range[1] - sigma_range[0]) * i / N
        #tau_lower = min(frac(3,2), (24 * sigma - 18) / (2 * sigma - 1))
        # tau_lower = min(frac(3,2), 48 * (8 * sigma ** 2 - 10 * sigma + 3) / (25 * sigma - 17))

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

            assert abs(res.x[2] - hypothesis.data.bound.at([sigma, tau])) < TOL

    print("Check complete")


# Tries to prove the bound LV(s, t) / t \leq f(s) on the specified domain defined by
# s \in sigma_range
# t \in tau_range(s)
def prove_LV_on_tau_bound(hypotheses, f, sigma_range, tau_range):
    pass


# Given a large-value estimate as a Hypothesis, apply Huxley subdivison (see Basic
# properties (ii) of Large value estimates section) to obtain a better large
# value estimate. The domain of the transformed hypothesis will be the same as 
# the original hypothesis. 
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
    raise NotImplementedError()
