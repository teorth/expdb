# code for representing a bound on the exponential sum
# \sum_{n \asymp N}e(TF(n/N)) \ll N^{\beta(\alpha) + o(1)}
# where $\alpha = \log N/\log T.

from constants import *
from hypotheses import *

# from mpmath import fraction as frac
from fractions import Fraction as frac
import numpy as np
from scipy.spatial import ConvexHull
from functions import *
import matplotlib.pyplot as plt
from reference import *


# Represents a single linear bound on the \beta(\alpha) function
class Bound_beta:
    # Creates a beta bound object.
    # Parameters:
    #   - bound:  (Affine object) a bound on \beta(\alpha) that is a affine function
    def __init__(self, bound):
        self.bound = bound

    def __repr__(self):
        return f"\\beta(x) \\leq {self.bound}"


class Bound_Beta_Transform:
    def __init__(self, name, transform):
        self.name = name
        self.transform = transform

    def __repr__(self):
        return self.name

###############################################################################


# Use this constructor to create a beta bound and add it to a set of hypotheses
def add_beta_bound(hypotheses, bounds, ref):
    for b in bounds:
        author = ref.entries["author"]
        year = ref.entries["year"]
        hypotheses.add_hypothesis(
            Hypothesis(
                f"{author} ({year}) bound on \\beta on {b.domain}",
                "Upper bound on beta",
                Bound_beta(b),
                f"See [{author}, {year}]",
                ref,
            )
        )


def derived_bound_beta(bound, proof, dependencies):
    year = Reference.max_year(tuple(d.reference for d in dependencies))
    bound = Hypothesis(
        "Derived bound on \\beta",
        "Upper bound on beta",
        Bound_beta(bound),
        proof,
        Reference.derived(year),
    )
    bound.dependencies = dependencies
    return bound


# use this constructor to create an upper bound on beta that is so classical it does not require citation
def classical_bound_beta(bound):
    return Hypothesis(
        f"Classical bound on \\beta",
        "Upper bound on beta",
        Bound_beta(bound),
        f"Classical",
        Reference.classical(),
    )


trivial_beta_bound_1 = classical_bound_beta(Affine(1, 0, Interval(0, 1, True, True)))
trivial_beta_bound_2 = classical_bound_beta(
    Affine(1, -1, Interval(1, Constants.ALPHA_UPPER_LIMIT, False, True))
)

# TODO: create a method to implement the B-process for beta bounds


# Given a set of hypotheses, compute the bounds on \beta(\alpha) implied by the
# assumed exponent pairs. Returns a list of Hypothesis objects representing derived
# bounds on \beta(\alpha).
#
# The exponent pair (k, l) implies \beta(\alpha) \leq (l - k)\alpha + k for
# 0 \leq \alpha \leq 1/2.
#
# Parameters:
#   - hypothesis_set: (Hypothesis_Set object)
# Returns:
#   - a list of Hypothesis objects, each representing a derived beta bound
# TODO: combine all beta bounds into a single beta bound (i.e. remove redundancies)
def exponent_pairs_to_beta_bounds(hypothesis_set):
    if not isinstance(hypothesis_set, Hypothesis_Set):
        raise "hypothesis_set must be of type Hypothesis_Set"

    hypotheses = []
    ephs = hypothesis_set.list_hypotheses(hypothesis_type="Exponent pair")
    transforms = hypothesis_set.list_hypotheses(hypothesis_type="Exponent pair to beta bound transform")

    domain = Interval(0, frac(1, 2), True, True)
    for eph in ephs:
        k = eph.data.k
        l = eph.data.l
        hypotheses.append(
            derived_bound_beta(
                Affine(l - k, k, domain),
                f'Follows from "{eph.name}"',
                {eph},
            )
        )

        # Append all beta bounds obtained via transformations
        for tr in transforms:
            hypotheses.extend(tr.data.transform(eph))

    return hypotheses


# Computes the best bound on \beta(\alpha) given a set of hypothesis, over a specified
# domain (domain is None by default, in which case [0, 1/2] is assumed). Returns
# the set of best bounds as a list of Hypothesis
def compute_best_beta_bounds(hypothesis_set, domain=None):
    bounds = hypothesis_set.list_hypotheses("Upper bound on beta")

    if domain is None:
        domain = Interval(0, frac(1, 2), include_lower=True, include_upper=True)

    # handle edge cases
    if len(bounds) < 2:
        return bounds

    # Create temporary deep copy of bounds to work with
    bds = [b.data.bound.deep_copy() for b in bounds]

    # Attach working id
    for i in range(len(bds)):
        bds[i].label = i

    # Iterate through all bounds
    best_bounds = [bds[0]]
    for i in range(1, len(bds)):
        best_bounds = bds[i].min_with(best_bounds, domain)

    # create temp dictionary for efficient lookup
    derived_bounds = []
    for b in best_bounds:
        depend = bounds[b.label]
        b.label = None  # label is no longer used - clear
        # if the bound was already in the hypotheses set, don't create a new object
        if b == depend.data.bound:
            derived_bounds.append(depend)
        else:
            derived_bounds.append(
                derived_bound_beta(b, f"Follows from {depend}", {depend})
            )
    return derived_bounds


# Displays a beta bound (a Hypothesis object), both in console and as a plot
def display_beta_bounds(hypotheses):

    # print to console
    print("\t\\beta(x) \\leq ")
    for h in hypotheses:
        print("\t\t", h.data.bound, "\t", h.name, "\t depends on", h.dependencies)

    # display in plot
    for h in hypotheses:
        p = h.data.bound
        plt.plot(
            (p.domain.x0, p.domain.x1),
            (p.at(p.domain.x0, True), p.at(p.domain.x1, True)),
            color="black",
        )
    plt.xlim((-0.05, 0.55))
    plt.ylim((-0.05, 0.5))
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta(\alpha)$")
    plt.title(r"Best bound on $\beta(\alpha)$")
    plt.show()
