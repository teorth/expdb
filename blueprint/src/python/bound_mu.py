# Code for handing upper bounds on $\mu(\sigma)$, defined as the best exponent in the bound
#
# $$ |\zeta(\sigma+it)| \ll (1+|t|)^{\mu(\sigma)+o(1)}.$$
#
# No support for lower bounds is implemented, since conjecturally (Lindelof hypothesis) the trivial lower bound is the truth.

from hypotheses import *
from exponent_pair import *

# from mpmath import fraction as frac
from fractions import Fraction as frac
import numpy as np
from scipy.spatial import ConvexHull
from reference import *


########################################################################################
# A Bound_mu is simply an upper bound on mu(sigma) at a particular choice of sigma.  It is raw data, containing no metadata;
# to upgrade the bound to a hypothesis containing metadata as well, use one of the constructors below.


class Bound_mu:
    def __init__(self, sigma, mu):
        self.sigma = sigma
        self.mu = mu

    def __repr__(self):
        return f"\\mu({self.sigma}) \\leq {self.mu}"

    def __eq__(self, other):
        if isinstance(other, Bound_mu):
            return (self.sigma, self.mu) == (other.sigma, other.mu)
        return NotImplemented


# use this constructor to create an upper bound on mu that is so classical it does not require citation
def classical_bound_mu(sigma, mu):
    return Hypothesis(
        f"Classical bound on \\mu({sigma})",
        "Upper bound on mu",
        Bound_mu(sigma, mu),
        f"Classical",
        Reference.classical(),
    )


# use this constructor to create an upper bound on mu from the literature
def literature_bound_mu(sigma, mu, ref):
    return Hypothesis(
        f"{ref.author()} bound on \\mu({sigma})",
        "Upper bound on mu",
        Bound_mu(sigma, mu),
        f"See [{ref.author()}, {ref.year()}]",
        ref,
    )


# use this constructor to create a conjectural upper bound hypothesis
def conjecture_bound_mu(sigma, mu, name):
    return Hypothesis(
        name,
        "Upper bound on mu",
        Bound_mu(sigma, mu),
        f"Conjecture",
        Reference.conjectured(),
    )


# use this constructor to create an upper bound hypothesis derived from other hypotheses
def derived_bound_mu(sigma, mu, proof, dependencies):
    year = Reference.max_year(tuple(d.reference for d in dependencies))
    bound = Hypothesis(
        f"Derived bound on \\mu({sigma})",
        "Upper bound on mu",
        Bound_mu(sigma, mu),
        proof,
        Reference.derived(year),
    )
    bound.dependencies = dependencies
    return bound


# this trivial bound will be automatically added to the optimization routine; no need to explicitly add them to the list of hypotheses
bound_mu_trivial = classical_bound_mu(1, 0)

Lindelof_hypothesis = conjecture_bound_mu(frac(1, 2), 0, "The Lindelof hypothesis")

# TODO: a method to clean up a list of mu bounds by removing any bound in the convex hull of existing bounds (including the trivial ones)


########################################################################################
# Methods for computing upper bounds on mu from a given list of hypotheses


# The trivial upper bound
def trivial_mu_bound(sigma):
    if sigma <= 0:
        return frac(1, 2) - sigma
    if sigma >= 1:
        return 0
    return frac(1, 2) - sigma / 2


def apply_trivial_mu_bound(sigma):
    return derived_bound_mu(
        sigma,
        trivial_mu_bound(sigma),
        f"This is the convexity bound on the Riemann zeta function",
        set(),
    )


# Apply the functional equation to a bound on mu to obtain another bound on mu
def apply_functional_equation(bound):
    new_sigma = 1 - bound.data.sigma
    new_mu = bound.data.mu + bound.data.sigma - frac(1, 2)
    return derived_bound_mu(
        new_sigma,
        new_mu,
        f'Follows from "{bound.name}" and the functional equation',
        {bound},
    )


# Apply the convexity of mu to two bounds on mu to obtain another bound on mu
def apply_mu_convexity(b1, b2, theta):
    if theta < 0 or theta > 1:
        print(f"FAILED: the parameter theta={theta} is not between 0 and 1.")
    mu = (1 - theta) * b1.data.mu + theta * b2.data.mu
    sigma = (1 - theta) * b1.data.sigma + theta * b2.data.sigma
    return derived_bound_mu(
        sigma,
        mu,
        f'This bound follows from convexity and the bounds "{b1}", "{b2}"',
        {b1, b2},
    )


# convert an exponent pair hypothesis to a bound on mu
def obtain_mu_bound_from_exponent_pair(bound):
    return derived_bound_mu(
        bound.data.l - bound.data.k,
        bound.data.k,
        f"This bound follows from the exponent pair ({bound.data.k}, {bound.data.l})",
        {bound},
    )


# Calculate the set of bounds on mu(sigma) implied by a hypothesis set
def get_bounds(hypothesis_set):
    # create the list of available bounds
    bounds = hypothesis_set.list_hypotheses("Upper bound on mu")
    bounds.append(bound_mu_trivial)

    # add bounds derived from exponent pairs. If (k, l) is an exponent pair,
    # then \mu(l - k) \leq k (see e.g. Ivic 1980)
    exp_pairs = compute_exp_pairs(hypothesis_set)
    bounds.extend(obtain_mu_bound_from_exponent_pair(p) for p in exp_pairs)

    # add bounds coming from the functional equation
    for bound in bounds:
        # mpmath does not currently support comparisons between Fraction and
        # the mpf object, use this for now to determine whether sigma > 1/2
        if 2 * bound.data.sigma > 1:
            bounds.append(apply_functional_equation(bound))
    return bounds


# Compute the 2-d convex hull containing all points (sigma, mu) implied
# by mu bounds in 'bounds' list and saves it in the 'hypothesis_set' object.
def compute_convex_hull(bounds, hypothesis_set):
    conv = ConvexHull(np.array([[b.data.sigma, b.data.mu] for b in bounds]))
    vertices = [bounds[v] for v in conv.vertices]
    hypothesis_set.data["convex_hull"] = vertices
    hypothesis_set.data_valid = True


# The best bound for mu(sigma) given a set of hypotheses, returned in the format
# of a new hypothesis
def best_mu_bound(sigma, hypothesis_set):
    # deal with trivial cases
    if sigma <= 0 or sigma >= 1:
        return apply_trivial_mu_bound(sigma)

    # Find the best bound on \mu(sigma).  ans is the current best candidate
    ans = apply_trivial_mu_bound(sigma)

    # Computed convex hull is stored within `hypothesis_list` the first time
    # to avoid repeatedly computing it for multiple values of sigma.
    if not hypothesis_set.data_valid or "convex_hull" not in hypothesis_set.data:
        # Generate set of bounds on mu implied by hypothesis set
        bounds = get_bounds(hypothesis_set)
        # The convex hull package requires us to provide at least 3 points, while some
        # hypothesis sets may only contain 2 points, so we include a placeholder point as
        # a third point to account for this edge case.
        pts = [b for b in bounds] + [conjecture_bound_mu(1 / 2, 10, "Placeholder")]
        compute_convex_hull(pts, hypothesis_set)

    # The vertices are guaranteed to be in counterclockwise order for 2D hulls
    verts = hypothesis_set.data["convex_hull"]
    for i in range(len(verts)):
        b1 = verts[i]
        b2 = verts[(i + 1) % len(verts)]
        if b1.data.sigma <= sigma and sigma <= b2.data.sigma:
            # Handle edge case
            if b1.data.sigma == b2.data.sigma:
                if b1.data.mu < b2.data.mu:
                    bound = b1
                else:
                    bound = b2
            else:
                theta = frac(sigma - b1.data.sigma, b2.data.sigma - b1.data.sigma)
                bound = apply_mu_convexity(b1, b2, theta)
            if bound.data.mu < ans.data.mu:
                ans = bound
    return ans


# TODO: supply not just a human-readable proof in the above method, but also a Lean-compilable proof.


# Computes the best bound on mu(sigma) in the range [sigma0, sigma1), as a piecewise
# linear function that the given hypotheses imply.
def best_mu_bound_piecewise(sigma_interval, hypothesis_set):
    sigma0, sigma1 = sigma_interval.x0, sigma_interval.x1

    if sigma0 > sigma1:
        raise ValueError("sigma0 must be <= sigma1")

    # deal with edge pieces
    if sigma0 < 0:
        return [("-x + 1/2", sigma0, 0)] + best_mu_bound_piecewise(
            Interval(0, sigma1, include_lower=True, include_upper=sigma_interval.include_upper),
            hypothesis_set
        )
    if sigma1 > 1:
        return best_mu_bound_piecewise(
            Interval(sigma0, 1, include_lower=sigma_interval.include_lower, include_upper=True),
            hypothesis_set) + [("0", 1, sigma1)]

    # Computed convex hull is stored within `hypothesis_set` the first time
    # to avoid repeatedly computing it for multiple values of sigma.
    if not hypothesis_set.data_valid or "convex_hull" not in hypothesis_set.data:
        bounds = get_bounds(hypothesis_set)
        compute_convex_hull(bounds, hypothesis_set)

    # The vertices are guaranteed to be in counterclockwise order for 2D hulls
    verts = hypothesis_set.data["convex_hull"]
    mu_bounds = []
    for i in range(len(verts)):
        b1 = verts[i].data
        b2 = verts[(i + 1) % len(verts)].data
        interval = Interval(max(sigma0, b1.sigma), min(sigma1, b2.sigma))
        if b1.sigma < b2.sigma and not interval.is_empty():
            mu_bounds.append(
                Affine(
                    (b2.mu - b1.mu) / (b2.sigma - b1.sigma),
                    (b2.sigma * b1.mu - b1.sigma * b2.mu) / (b2.sigma - b1.sigma),
                    Interval(max(sigma0, b1.sigma), min(sigma1, b2.sigma))
                )
            )
    return mu_bounds


# This method attempts to prove the bound mu(sigma) <= mu (or better) from the set of hypotheses.
# If successful, it adds the bound to the set; otherwise, it returns an error.
def prove_mu_bound(sigma, mu, hypothesis_list):
    bound = best_mu_bound(sigma, hypothesis_list)
    if bound.data.mu <= mu:
        # add the new hypothesis without requiring recomputation of cached data
        hypothesis_list.add_hypothesis(bound, invalidate_data=False)
        print(f"Proved {bound.description}.  {bound.proof}.")
        return bound
    else:
        print(
            f"FAILED to prove \\mu({sigma}) <= {mu}.  The best bound found on \\mu({sigma}) was {bound.data.mu}."
        )
        # should throw some exception here
        return False
