# code for representing a an upper bound on how frequently a Dirichlet polynomial
# of length N can be large. For

from constants import *
import copy
from fractions import Fraction as frac
from functions import *
from helpers.str_helper import Str_Helper
from hypotheses import *
import itertools
from polytope import *
from reference import *
from region import Region
import scipy
import time


###############################################################################

class Large_Value_Estimate:

    """
    Class representing a large value estimate ρ \\le LV(σ, τ).

    Internally a large value estimate is stored as a 3-dimensional Region 
    representing the set of feasible tuples of (σ, τ, ρ). 
    """

    def __init__(self, region:Region, repr:str=None):
        
        """
        Constructs a large value estimate.

        Parameters
        ----------
        region : Region
            The 3-dimensional region representing a set of feasible (σ, τ, ρ)
            values 
        repr : str, optional
            The string representation for this object. Defaults to None. 
        """

        if not isinstance(region, Region):
            raise ValueError("Parameter region must be of type Region")
        
        self.region = region
        self.repr = repr

    def __repr__(self):
        # If there is a specified string representation for this object, 
        # return it
        if self.repr is not None:
            return self.repr
        
        # Otherwise: use specially formatted Region object
        s = self.region.to_str(use_indentation=False, variables="στρ")
        return f"(σ,τ,ρ) in {s}"


class Large_Value_Estimate_Transform:

    """
    Class representing a large value estimate transform 
    """

    def __init__(self, transform):
        self.transform = transform

    def __repr__(self):
        return "Raising to a power"


class Large_Value_Energy_Estimate:
    
    """
    Class representing a additive energy estimate ρ* \\le LV*(σ, τ)
    """
    def __init__(self, region):
        if not isinstance(region, Region):
            raise ValueError("Parameter region must be of type Region")
        self.region = region
        
    def __repr__(self):
        return f"(σ,τ,ρ*) in {self.region}"

###############################################################################

def convert_bounds_to_region(bounds):
    # The default domain of definition in (\sigma, \tau, \rho) space
    box = [
        [-frac(1,2), 1, 0, 0],                      # sigma >= 1/2
        [1, -1, 0, 0],                              # sigma <= 1
        [0, 0, 1, 0],                               # tau >= 0
        [Constants.TAU_UPPER_LIMIT, 0, -1, 0],      # tau <= TAU_UPPER_LIMIT
        [0, 0, 0, 1],                               # rho >= 0
        [Constants.LV_DEFAULT_UPPER_BOUND, 0, 0, -1]# rho <= RHO_UPPER_LIMIT
    ]
    # Convert bounds from 'bounds', of the form 
    # \rho \le f(\sigma, \tau) 
    # into a region in R^3 of the form 
    # f(\sigma, \tau, \rho) \ge 0.
    return Region.from_union_of_halfplanes(
        [b + [-1] for b in bounds],
        box
    )

def literature_bound_LV_max(bounds, ref, params=""):
    # Compute a string representation of bounds
    repr = f"ρ <= max({', '.join(Str_Helper.format(b, ['σ','τ']) for b in bounds)})"

    return Hypothesis(
        f"{ref.author()} large value estimate" + params,
        "Large value estimate",
        Large_Value_Estimate(convert_bounds_to_region(bounds), repr=repr),
        f"See [{ref.author()}, {ref.year()}]",
        ref,
    )

def derived_bound_LV(region, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived large value estimate",
        "Large value estimate",
        Large_Value_Estimate(region),
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound

def classical_LV_estimate(bounds):
    # Compute a string representation of bounds
    repr = f"ρ <= max({', '.join(Str_Helper.format(b, ['σ','τ']) for b in bounds)})"

    return Hypothesis(
        "Classical large value estimate",
        "Large value estimate",
        Large_Value_Estimate(convert_bounds_to_region(bounds), repr=repr),
        "Classical",
        Reference.classical(),
    )

def conjectured_LV_estimate(bounds, name):
    # Compute a string representation of bounds
    repr = f"ρ <= max({', '.join(Str_Helper.format(b, ['σ','τ']) for b in bounds)})"
    
    return Hypothesis(
        name,
        "Large value estimate",
        Large_Value_Estimate(convert_bounds_to_region(bounds), repr=repr),
        "Conjecture",
        Reference.conjectured(),
    )

###############################################################################

# Classical estimates

# L^2 mean value theorem: LV(s, t) \leq max(2 - 2s, 1 - 2s + t)
large_value_estimate_L2 = classical_LV_estimate([[2, -2, 0], [1, -2, 1]])

# Montgometry's conjecture
montgomery_conjecture = conjectured_LV_estimate([[2, -2, 0]], "Montgomery conjecture")

###############################################################################

# Raising to a power (large value estimate transform theorem)
def raise_to_power_hypothesis(k):

    # Convert into fraction
    k = frac(k)

    # LV(s, kt) \leq kLV(s, t)
    def transform(hypothesis):
        # Scale while re-imposing upper limit on tau
        new_region = hypothesis.data.region.scale_all([1, k, k])
        return derived_bound_LV(
            new_region,
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

def optimize_bourgain_large_value_estimate():
    """
    Finds the optimal choice of alpha_1, alpha_2 in Bourgain's large values theorem
    (Corollary 10.30 in the web blueprint) by solving a series of linear programs. 

    The optimal choice of alpha_1, alpha_2 is a piecewise-defined function of sigma
    and tau. 
    """

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

    box = Polytope.rect(
        (frac(1,2), frac(1)), 
        (frac(1), frac(3))
        (0, Constants.LV_DEFAULT_UPPER_BOUND)
    )

    # Sympy solvers give empty solution sets for these systems - so instead we
    # (temporarily) use sympy's matrix rref computer
    # a1, a2, s, t, M = sympy.symbols("a1, a2, s, t, M")
    # var = [a1, a2, s, t, M]
    regions = []
    dependencies = []

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
        domain = Polytope([a1_defn, a2_defn]).lift(
            0, 1, (0, Constants.LV_DEFAULT_UPPER_BOUND)
        )

        if not domain.is_empty(include_boundary=False):
            # Hack to ensure uniqueness (using the fact that tuples are hashable)
            fns = set()
            for i in range(6):
                fn = [eqns[i][5], eqns[i][2], eqns[i][3]]
                fn = tuple(fn[j] + eqns[i][0] * a1_defn[j] + eqns[i][1] * a2_defn[j] for j in range(len(fn)))
                fns.add(fn)

            # Compute feasible (sigma, tau, rho) region
            union = [
                domain.intersect(Polytope([list(fn) + [-1]])) # Region where f(sigma, tau) - rho >= 0
                for fn in fns
            ]

            # Deal with the unbounded region
            neg = box.set_minus(domain)
            if not neg.is_empty(include_boundary=False):
                union.extend(neg)

            regions.append(Region.union(union))
        else:
            regions.append(box)

        a1_proof = "a1 = " + Affine2.to_string(a1_defn, "st")
        a2_proof = "a2 = " + Affine2.to_string(a2_defn, "st")
        dependencies.append(f"Follows from taking {a1_proof} and {a2_proof}")

    # Compute intersection
    R = Region.intersect(regions)
    R = R.to_disjoint_union()

    # TODO: implement dependency tracking

    return R

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
