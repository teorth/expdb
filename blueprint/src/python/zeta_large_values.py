# This class handles bounds on the function LV_{\zeta}(\sigma, \tau)

import copy
from hypotheses import *
from polytope import *
from functions import *
from constants import *
import bound_beta as bbeta
import large_values as lv
from reference import Reference
from region import Region, Region_Type

###############################################################################


# Returns the trivial zeta large value estimates:
# LV_zeta(s, t) \leq t for 1/2 \leq \s \leq 1, t \geq 0
def get_trivial_zlv():
    # Trivial estimate
    domain = Polytope(
        [
            [-frac(1, 2), 1, 0],  # \sigma > 1/2
            [1, -1, 0],  # \sigma <= 1
            [-2, 0, 1],  # \tau >= 2
            [Constants.TAU_UPPER_LIMIT, 0, -1],  # \tau <= some large number
        ]
    )
    trivial = Hypothesis(
        "Trivial zeta large value estimate",
        "Zeta large value estimate",
        lv.Large_Value_Estimate(Piecewise([Affine2([0, 0, 1], domain)])),
        "Trivial",
        Reference.trivial(),
    )
    return trivial


def literature_bound_ZLV_Max(bounds, ref, params=""):
    hyp = lv.literature_bound_LV_max(bounds, ref, params)
    return Hypothesis(
        f"{ref.author()} ({ref.year()}) large value estimate{params}",
        "Zeta large value estimate",
        hyp.data,
        f"See [{ref.author(), ref.year()}]",
        ref,
    )


def derived_bound_zeta_LV(bound, proof, deps):
    year = Reference.max_year(tuple(d.reference for d in deps))
    bound = Hypothesis(
        "Derived zeta large value estimate",
        "Zeta large value estimate",
        lv.Large_Value_Estimate(bound),
        proof,
        Reference.derived(year),
    )
    bound.dependencies = deps
    return bound


###############################################################################


# Compute the set of best zeta large values estimates using the hypotheses
def compute_large_value_estimate(hypotheses):
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("hypotheses must be of type Hypothesis_Set")

    # Generate set of LV estimates (original + transformed)
    lv_estimates = list(
        hypotheses.list_hypotheses(hypothesis_type="Zeta large value estimate")
    )
    lv_estimates.extend(lv_to_zlv(hypotheses))
    lv_estimates.extend(beta_to_zlv(hypotheses))
    # lv_estimates.extend(mu_to_zlv(hypotheses))

    return lv_estimates


# Use existing large value estimates in the hypothesis set to derive zeta large
# value estimates
def lv_to_zlv(hypotheses):
    lvs = hypotheses.list_hypotheses("Large value estimate")
    lv_transforms = hypotheses.list_hypotheses("Large value estimate transform")

    hyps = list(lvs)
    for tr_hyp in lv_transforms:
        hyps.extend([tr_hyp.data.transform(lve) for lve in lvs])

    # As a convention we restrict the domain to tau >= 2, while large value estimates
    # are restricted to tau >= 0
    domain = Region.from_polytope(Polytope([[-2, 0, 1, 0]]))
    zlvs = []
    for lve in hyps:
        region = Region.intersect([domain, lve.data.region])
        zlvs.append(derived_bound_zeta_LV(region, f"Follows from {lve.name}", {lve}))
    return zlvs


# Use bounds on beta to obtain new zeta large value estimates
# For t \geq 0, LV(s, t) = 0 whenever s \geq t \beta(1/t)
# i.e. given a affine bound \beta(\alpha) \leq A + B\alpha (a0 < \alpha < a1),
# we have
# LV(s, t) = 0 for s \geq A t + B (1/a1 < t < 1/a0)
def beta_to_zlv(hypotheses):

    # Compute the best beta_estimates
    beta_estimates = bbeta.compute_best_beta_bounds(hypotheses)

    zlv_estimates = []
    for be in beta_estimates:
        A = be.data.bound.c
        B = be.data.bound.m
        alpha0 = be.data.bound.domain.x0
        alpha1 = be.data.bound.domain.x1

        tau_upper = (
            Constants.TAU_UPPER_LIMIT
            if alpha0 == 0
            else min(Constants.TAU_UPPER_LIMIT, 1 / alpha0)
        )
        tau_lower = (
            Constants.TAU_UPPER_LIMIT if alpha1 == 0 else max(frac(2, 1), 1 / alpha1)
        )

        # In order for the min_with function to work we require that we have a
        # piecewise function defined on the whole domain. Only the first of
        # these pieces if the bound, the other 3 pieces complete the domain with
        # the trivial bound LVZ(s, t) \leq t
        polys = [
            Polytope([
                [0, 0, 0, -1],          # \rho <= 0 
                [0, 0, 0, 1],           # \rho >= 0
                [1, -1, 0, 0],          # \sigma <= 1
                [-frac(1, 2), 1, 0, 0], # \sigma >= 1/2
                [-tau_lower, 0, 1, 0],  # \tau >= max(2, 1/a1)
                [tau_upper, 0, -1, 0],  # \tau <= min(some large number, 1/a0)
                [-B, 1, -A, 0],         # \sigma >= A\tau + B
            ]),
            # Trivial bound elsewhere -----------------------------------------
            Polytope([
                [0, 0, 1, -1],          # \rho <= \tau
                [0, 0, 0, 1],           # \rho >= 0
                [1, -1, 0, 0],          # \sigma <= 1
                [-frac(1, 2), 1, 0, 0], # \sigma >= 1/2
                [-tau_lower, 0, 1, 0],  # \tau >= max(2, 1/a1)
                [tau_upper, 0, -1, 0],  # \tau <= min(some large number, 1/a0)
                [B, -1, A, 0],          # \sigma <= A\tau + B
            ]),
            Polytope([
                [0, 0, 1, -1],          # \rho <= \tau
                [0, 0, 0, 1],           # \rho >= 0
                [1, -1, 0, 0],          # \sigma <= 1
                [-frac(1, 2), 1, 0, 0], # \sigma >= 1/2
                [tau_lower, 0, -1, 0],  # \tau <= max(2, 1/a1)
                [-2, 0, 1, 0],          # \tau >= 2
            ]),
            Polytope([
                [0, 0, 1, -1],          # \rho <= \tau
                [0, 0, 0, 1],           # \rho >= 0
                [1, -1, 0, 0],          # \sigma <= 1
                [-frac(1, 2), 1, 0, 0], # \sigma >= 1/2
                [-tau_upper, 0, 1, 0],  # \tau >= min(some large number, 1/a0)
                [
                    Constants.TAU_UPPER_LIMIT,
                    0,
                    -1,
                    0
                ],  # \tau <= TAU_UPPER_LIMIT
            ])
        ]
        region = Region(
            Region_Type.DISJOINT_UNION, 
            [Region.from_polytope(p) for p in polys]
        )
        zlv_estimates.append(derived_bound_zeta_LV(region, f"Follows from {be.name}", {be}))

    return zlv_estimates


# Use bounds on the function \mu(\sigma) to obtain new zeta large value estimates
def mu_to_zlv(hypotheses):

    mu_bounds = hypotheses.list_hypotheses("Upper bound on mu")

    # TODO: compute the best mu bounds first?
    zlv_estimates = []
    for mb in mu_bounds:
        sigma0 = mb.data.sigma
        mu0 = mb.data.mu
        polys = [
            Polytope([
                [0, 0, 0, 1],           # \rho >= 0
                [0, 0, 0, -1],          # \rho <= 0
                [1, -1, 0, 0],          # \sigma <= 1
                [-frac(1, 2), 1, 0, 0], # \sigma > 1/2
                [-2, 0, 1, 0],          # \tau >= 2
                [
                    Constants.TAU_UPPER_LIMIT,
                    0,
                    -1,
                    0
                ],  # \tau <= min(some large number, 1/a0)
                [-sigma0, 1, -mu0, 0],  # \sigma >= sigma0 + \tau * mu0
            ]),
            # Trivial bound elsewhere -----------------------------------------
            Polytope([
                [0, 0, 1, -1],          # \rho <= \tau
                [0, 0, 0, 1],           # \rho >= 0
                [1, -1, 0, 0],          # \sigma <= 1
                [-frac(1, 2), 1, 0, 0], # \sigma > 1/2
                [-2, 0, 1, 0],          # \tau >= 2
                [
                    Constants.TAU_UPPER_LIMIT,
                    0,
                    -1,
                    0
                ],  # \tau <= min(some large number, 1/a0)
                [sigma0, -1, mu0, 0],  # \sigma <= sigma0 + \tau * mu0
            ])
        ]
        region = Region(
            Region_Type.DISJOINT_UNION, 
            [Region.from_polytope(p) for p in polys]
        )
        zlv_estimates.append(derived_bound_zeta_LV(region, f"Follows from {mb.name}", {mb}))

    return zlv_estimates
