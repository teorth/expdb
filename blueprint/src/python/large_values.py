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
def max_of(bounds):

    # The standard domain of definition is \sigma \in [1/2, 1], \tau \geq 0
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
        s = list(domain_ineq)  # shallow copy
        for j in range(len(b)):
            if b[j] == "1":
                s.append(intersections[j].constraint.coefficients)
            else:
                s.append([-x for x in intersections[j].constraint.coefficients])
        p = Polytope(s, canonicalize=True)
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


# Try to prove the given large value estimate (hypothesis) by assuming the hypotheses
# in hypothesis_set
def prove_large_value_estimate(hypothesis, hypothesis_set):
    if not isinstance(hypothesis, Hypothesis):
        raise ValueError("Parameter hypothesis must be of type Hypothesis")
    if not isinstance(hypothesis_set, Hypothesis_Set):
        raise ValueError("Parameter hypothesis_set must be of type Hypothesis_Set")
