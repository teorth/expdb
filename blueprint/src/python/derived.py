# This file contains examples of bounds that can be derived from other bounds,
# including those from literature.py

from literature import *
import exponent_pair as ep
import zero_density_estimate as zd
import prime_gap as pg
from typing import Optional

######################################################################################
# Derivations of exponent pairs

def van_der_corput_pair(k: int) -> Hypothesis:
    """
    Prove the classical van der Corput exponent pair (1/(2^k-2), 1 - (k-1)/(2^k-2)),
    returning the proof as a Hypothesis object.
    """
    if k < 2:
        raise ValueError("k must be at least 2.")
    A_transform = literature.find_hypothesis(keywords="van der Corput A transform")
    B_transform = literature.find_hypothesis(keywords="van der Corput B transform")
    exp_pair = B_transform.data.transform(trivial_exp_pair)
    for _ in range(k - 2):
        exp_pair = A_transform.data.transform(exp_pair)
    print(f"The van der Corput pair for k = {k} is {exp_pair.desc()}")
    return exp_pair

def prove_hardy_littlewood_mu_bound() -> Hypothesis:
    """
    Prove the Hardy-Littlewood bound \\mu(1/2) \\leq 1/6 using the van der Corput
    exponent pair (1/6, 2/3), returning the proof as a Hypothesis object.
    """
    HL_bound = literature.find_hypothesis(data=Bound_mu(frac(1, 2), frac(1, 6)))
    A_transform = literature.find_hypothesis(keywords="van der Corput A transform")
    B_transform = literature.find_hypothesis(keywords="van der Corput B transform")
    print(f"We will reprove {HL_bound.desc()}.")
    B_exp_pair = B_transform.data.transform(trivial_exp_pair)
    print(f"We have {B_exp_pair.desc_with_proof()}")
    AB_exp_pair = A_transform.data.transform(B_exp_pair)
    print(f"This implies {AB_exp_pair.desc_with_proof()}")
    mu_bound = exponent_pair_to_mu_bound(AB_exp_pair)
    print(f"This implies {mu_bound.desc_with_proof()}")
    return mu_bound

def prove_exponent_pair_with_hypotheses(hypotheses:Hypothesis_Set, k:frac, l:frac):
    """
    Tries to find a proof of the exponent pair (k, l) by assuming a specific
    set of hypotheses.

    Parameters
    ----------
    hypotheses : Hypothesis_Set
        The set of hypotheses to assume.
    k : Fraction
        The first element of the exponent pair (k, l).
    l : Fraction
        The second element of the exponent pair (k, l).

    Returns
    -------
    Hypothesis or None
        If a proof is found, then the exponent pair is returned as a Hypothesis.
        Otherwise, None is returned.
    """

    # Expand the set of exponent pairs using transforms
    hypotheses.add_hypotheses(compute_exp_pairs(hypotheses, search_depth=1))

    # Perform 1 iteration of exponent pair -> beta bounds -> exponent pair
    hypotheses.add_hypotheses(exponent_pairs_to_beta_bounds(hypotheses))
    hypotheses.add_hypotheses(compute_best_beta_bounds(hypotheses))
    new_exp_pairs = beta_bounds_to_exponent_pairs(hypotheses)

    return next((h for h in new_exp_pairs if h.data.k == k and h.data.l == l), None)

def prove_exponent_pair(
        k:frac,
        l:frac,
        simplify_deps:bool = False,
        verbose:bool = False
    ) -> Hypothesis:

    """
    Tries to find a proof of the exponent pair (k, l) using results from the
    literature. For the analogous method given a specific set of hypotheses,
    see prove_exponent_pair_with_hypotheses().

    Parameters
    ----------
    k : Fraction
        The first element of the exponent pair (k, l).
    l : Fraction
        The second element of the exponent pair (k, l).
    simplify_deps : bool, optional
        If True, then the proof will be simplified by removing redundant dependencies.
        This is similar to applying ep.find_best_proof() with
        method = Proof_Optimization_Method.COMPLEXITY.
    verbose : bool, optional
        If True, then additional info will be logged to console. Default is False.

    Returns
    -------
    Hypothesis or None
        If a proof is found, then the exponent pair is returned as a Hypothesis.
        Otherwise, None is returned.
    """

    hypotheses = Hypothesis_Set()

    hypotheses.add_hypothesis(trivial_exp_pair)

    hyp_types = ["Upper bound on beta", "Exponent pair", "Exponent pair transform"]
    for ht in hyp_types:
        hypotheses.add_hypotheses(literature.list_hypotheses(hypothesis_type=ht))

    eph = prove_exponent_pair_with_hypotheses(Hypothesis_Set(hypotheses), k, l)

    # If simplification required, greedily remove dependencies
    if simplify_deps:
        for i in range(100):
            cpy = Hypothesis_Set(hypotheses)
            # Randomly remove an element from the hypothesis set, if the reduced set of
            # hypotheses is still sufficient, then update the set of hypotheses
            hyp = cpy.hypotheses.pop()
            eph1 = prove_exponent_pair_with_hypotheses(Hypothesis_Set(cpy), k, l)
            if eph1 is not None:
                hypotheses = cpy
                eph = eph1

    # Log to console
    if verbose:
        if eph is None:
            print(f"Failed to prove the exponent pair ({k}, {l}).")
        else:
            print(f"Proof of the ({k}, {l}) exponent pair:")
            eph.recursively_list_proofs()

    return eph

def prove_heathbrown_exponent_pairs():

    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(trivial_exp_pair)
    hypotheses.add_hypotheses(
        l
        for l in literature.list_hypotheses(hypothesis_type="Upper bound on beta")
        if l.reference.author() == "Heath-Brown"
    )

    for h in hypotheses.list_hypotheses():
        print(h.desc_with_proof())

    print("-------------------------------------------------------------------")
    best_beta_bounds = compute_best_beta_bounds(hypotheses)
    for h in best_beta_bounds:
        print(h.desc_with_proof())
    hypotheses.add_hypotheses(best_beta_bounds)
    new_exp_pairs = beta_bounds_to_exponent_pairs(hypotheses)
    for ep in new_exp_pairs:
        print(ep)

def best_proof_of_exponent_pair(
        k: frac,
        l: frac,
        proof_method: int=Proof_Optimization_Method.DATE,
        verbose: bool=True
    ) -> Optional[Hypothesis]:

    """
    Finds the best proof of the exponent pair (k, l) according to some measure
    of "goodness" of a proof. Examples of such measures include the date of the
    latest dependency of the proof, or the "complexity" of a proof measured by
    the total number of recursive dependencies.

    Parameters
    ----------
    k: frac
        The first element of the exponent pair (k, l).
    l: frac
        The second element of the exponent pair (k, l).
    proof_method: int, optional
        The measure of the goodness of a proof. Possible values are given in the
        class Proof_Optimization_Method, e.g. Proof_Optimization_Method.DATE,
        Proof_Optimization_Method.COMPLEXITY. Default is Proof_Optimization_Method.DATE.
    verbose: bool, optional
        If True, additional debugging info will be logged to console. Default is True.

    Returns
    -------
    Hypothesis or None
        A hypothesis object representing the exponent pair (k, l) if a proof
        is found. Otherwise, None is returned.

    """
    hypothesis = copy.copy(literature)
    hypothesis.add_hypothesis(ep.trivial_exp_pair)

    hyp = ep.find_best_proof(
        k, l, hypothesis, proof_method
    )
    if verbose:
        print()
        if hyp is not None:
            print(f'Found proof of ({k}, {l}) with complexity = {hyp.proof_complexity()} and date = {hyp.proof_date()}:')
            hyp.recursively_list_proofs()
        else:
            print(f'Failed to prove the exponent pair ({k}, {l}).')
    return hyp

def prove_exponent_pairs():
    # prove_heathbrown_exponent_pairs()
    prove_exponent_pair(frac(89,1282), frac(997,1282), simplify_deps=False, verbose=True)
    prove_exponent_pair(frac(652397,9713986), frac(7599781,9713986), simplify_deps=False, verbose=True)
    prove_exponent_pair(frac(10769,351096), frac(609317,702192), simplify_deps=False, verbose=True)
    prove_exponent_pair(frac(89,3478), frac(15327,17390), simplify_deps=False, verbose=True)

    best_proof_of_exponent_pair(frac(1, 6), frac(2, 3))
    best_proof_of_exponent_pair(frac(13, 31), frac(16, 31))
    best_proof_of_exponent_pair(frac(4, 11), frac(6, 11))
    best_proof_of_exponent_pair(frac(2, 7), frac(4, 7))
    best_proof_of_exponent_pair(frac(5, 24), frac(15, 24))
    best_proof_of_exponent_pair(frac(4, 18), frac(11, 18))
    best_proof_of_exponent_pair(frac(3, 40), frac(31, 40), Proof_Optimization_Method.DATE)
    #best_proof_of_exponent_pair(frac(3, 40), frac(31, 40), Proof_Optimization_Method.COMPLEXITY)

######################################################################################
# Derivations of bounds on zeta growth rates in the critical strip (mu bounds)

def compute_best_mu_bound(
        sigma_interval: Interval = Interval(frac(1,2), frac(99,100))
    ) -> list[Affine]:
    """
    Compute the best-known mu bounds from known exponent pairs and bounds on the
    beta function.

    Parameters
    ----------
    sigma_interval: Interval
        The values of sigma on which the bound on mu(sigma) is computed.
    """

    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(trivial_exp_pair)

    assume = [
        "Upper bound on beta",
        "Exponent pair",
        "Exponent pair transform",
        "Exponent pair to beta bound transform"
    ]
    for hyp_type in assume:
        hyps = literature.list_hypotheses(hypothesis_type=hyp_type)
        print(f"Assuming {len(hyps)} hypotheses of type {hyp_type}.")
        hypotheses.add_hypotheses(hyps)

    # Expand the hypothesis set using a few iterations of
    # exp pairs <-> beta bounds
    hypotheses.add_hypotheses(compute_exp_pairs(hypotheses, search_depth=1))
    hypotheses.add_hypotheses(exponent_pairs_to_beta_bounds(hypotheses))
    hypotheses.add_hypotheses(compute_best_beta_bounds(hypotheses))
    hypotheses.add_hypotheses(beta_bounds_to_exponent_pairs(hypotheses))

    print("The following bounds on mu were obtained in the range [1/2, 99/100]")
    mbs = best_mu_bound_piecewise(sigma_interval, hypotheses)
    for b in mbs:
        print(f"\\mu(x) \\leq {b}. {b.domain.x1} = {float(b.domain.x1)}")
    return mbs


######################################################################################
# Derivations of large value estimates

def prove_bourgain_large_values_theorem():
    """
    Proves the optimized form of Bourgain (2000) large values theorem by optimizing the
    choice of alpha_1, alpha_2 as functions of sigma and tau.
    """
    lv_hyps = lv.optimize_bourgain_large_value_estimate()
    pieces = []
    for lvh in lv_hyps:
        pieces.extend(lvh.data.bound.pieces)
    bound = Piecewise(pieces)
    bound.plot_domain((1/2, 1), (1, 3))

    for lvh in lv_hyps:
        print(lvh.data, lvh.proof)

def prove_guth_maynard_large_values_theorem():
    """
    Prove the Guth--Maynard large values theorem (Theorem 1.1 in their paper) using the
    Large Value Energy Region estimates due to Guth--Maynard
    """
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value energy region 1 with k = 2"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value energy region 2"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value energy region 3"))

    # Compute the feasible region of (sigma, tau, rho, rho*, s) values
    # in the domain 1/2 <= sigma <= 1, tau >= 0
    lver_region = ad.compute_best_lver(
        hypotheses,
        Region.from_polytope(Polytope.rect((frac(1,2), 1), (0, Constants.TAU_UPPER_LIMIT))),
        zeta=False
    )

    # Project into a feasible region of (sigma, tau, rho) values
    region = ad.lver_to_lv(lver_region)

    # Take \tau = 6/5 (TODO: replace this step with Huxley subdivision once it is implemented)
    region = region.data.region.substitute({1: frac(6,5)})

    # Constrain to the range 7/10 \leq \sigma \leq 8/10 (\rho unconstrained)
    domain = Polytope([
        [-frac(7,10), 1, 0], [frac(8,10), -1, 0]
    ])
    region.child = [Region.from_polytope(r.child.intersect(domain)) for r in region.child]
    region.child = [r for r in region.child if not r.child.is_empty(include_boundary=False)]

    # Simplify the region
    poly = Polytope.try_union([r.child for r in region.child])

    print("Proved feasible region for (σ, τ, ρ):", poly.to_str("στρ"))

def prove_guth_maynard_intermediate_lvt():
    """
    Prove an intermediate large value theorem from in Guth--Maynard
    using other intermediate large values estimates also due to Guth--Maynard
    """
    # Separate optimization for each value of k
    for k in range(1, Constants.LARGE_VALUES_TRUNCATION):
        hs = Hypothesis_Set()
        hs.add_hypothesis(literature.find_hypothesis(name=f"Guth--Maynard large value estimate 2 with k = {k}"))
        hs.add_hypothesis(literature.find_hypothesis(name=f"Jutila large value estimate with k = 3"))
        hs.add_hypothesis(lv.large_value_estimate_L2)
        hs.add_hypothesis(literature.find_hypothesis(name=f"Huxley large value estimate"))

        # Proof of first theorem: in the region σ >= 7/10, 1 <= τ <= 6/5
        rect = Region.from_polytope(
            Polytope.rect(
                (frac(7,10), 1),
                (1, frac(6,5)),
                (0, Constants.LV_DEFAULT_UPPER_BOUND)
            )
        )
        region = Region.intersect(
            [rect] +
            [h.data.region for h in hs.list_hypotheses(hypothesis_type="Large value estimate")]
        )
        region = region.as_disjoint_union()

        # hypothesized large value estimate
        region2 = lv.convert_bounds(
                [
                    [2, -2, 0],
                    [3, -4, frac(1,2)],
                    [frac(4*k, k+1), -frac(6*k, k+1), frac(k, k+1)],
                    [frac(20*k, 4*k+3), -frac(24*k, 4*k+3), frac(2, 4*k+3)]
                ]
            ).region

        # Check that region is contained in region2
        print(f"Checking hypothesis holds for k = {k}:", region.is_subset_of(region2))

def prove_guth_maynard_intermediate_lvt2():
    """
    Prove another intermediate large value theorem from in Guth--Maynard
    using other intermediate large values estimates also due to Guth--Maynard
    (see also prove_guth_maynard_intermediate_lvt())

    Note: one needs to take k sufficiently large for this verification to
    succeed. The current maximum k value of 99 is sufficient for the current
    sample of sigma, however may be insufficient for other choices of sigma.
    If the verification fails for a specific value of sigma, one could try
    increasing the maximum value of k.
    """
    # Second lvt estimate check (using numerical sampling of sigma)
    sigma_interval = (frac(7,10), frac(1))
    sigma_N = 100
    for i in range(sigma_N + 1):
        sigma = sigma_interval[0] + (sigma_interval[1] - sigma_interval[0]) / sigma_N * i

        regions = []
        for k in range(1, 100): # Need to take k sufficiently high
            region = lv.convert_bounds(
                [
                    [2, -2, 0],
                    [3, -4, frac(1,2)],
                    [frac(4*k, k+1), -frac(6*k, k+1), frac(k, k+1)],
                    [frac(20*k, 4*k+3), -frac(24*k, 4*k+3), frac(2, 4*k+3)]
                ]
            ).region.substitute({0: sigma})
            regions.append(region)

        # L2 mean-value theorem
        L2 = lv.large_value_estimate_L2.data.region.substitute({0: sigma})

        # Huxley large value estimate
        huxley = literature.find_hypothesis(name="Huxley large value estimate")\
                        .data.region.substitute({0: sigma})

        # tau-rho domain
        domain = Region.from_polytope(
            Polytope.rect(
                (1, frac(6,5)), # 1 <= τ <= 6/5
                (0, Constants.LV_DEFAULT_UPPER_BOUND) # 0 <= ρ
            )
        )

        # Compute their intersection
        region1 = Region.intersect([domain, L2, huxley] + regions).as_disjoint_union()

        # Hypothesis that we want to test (a 2-d polytope in (τ, ρ))
        region2 = Region.from_union_of_halfplanes(
            # Halfplanes
            [
                [2 - 2 * sigma, 0, -1], # ρ <= 2 - 2σ
                [3 - 4 * sigma, frac(1,2), -1], # ρ <= 3 - 4σ + τ/2
                [(46 - 60 * sigma) / 5, (30 * sigma - 21) / 5, -1] # ρ <= (46 - 60σ)/5 + (30σ - 21)τ/5
            ],
            # Bounding box constraints
            [
                (-1, 1, 0), # τ >= 1
                (frac(6,5), -1, 0), # τ <= 6/5
                (0, 0, 1), # ρ >= 0
                (Constants.LV_DEFAULT_UPPER_BOUND, 0, -1), # ρ <= some large number
            ]
        )
        print(f"Verifying Guth--Maynard large value estimate for σ = {sigma} ~ {float(sigma)}:", region1.is_subset_of(region2))

def prove_guth_maynard_lvt_from_intermediate_lemmas():
    """
    Prove Theorem 1.1 of Guth--Maynard using intermediate large value estimates found in
    Guth--Maynard
    """
    region = Region.intersect([
        lv.convert_bounds(
                [
                    [2, -2, 0],
                    [3, -4, frac(1,2)],
                    [frac(4*k, k+1), -frac(6*k, k+1), frac(k, k+1)],
                    [frac(20*k, 4*k+3), -frac(24*k, 4*k+3), frac(2, 4*k+3)]
                ]
            ).region
        for k in range(1, Constants.LARGE_VALUES_TRUNCATION)
    ])

    region = region.as_disjoint_union()

    # hypothesized large value estimate
    hypothesis = literature.find_hypothesis(name="Guth--Maynard large value estimate")
    region2 = hypothesis.data.region

    print(f"Checking that hypothesis holds:", region.is_subset_of(region2))

def prove_all_large_value_estimates():
    prove_bourgain_large_values_theorem()
    prove_guth_maynard_large_values_theorem()
    prove_guth_maynard_intermediate_lvt()
    prove_guth_maynard_lvt_from_intermediate_lemmas()

######################################################################################
# Derivations of zero-density estimates for the Riemann zeta-function

def prove_zero_density(
        additional_hypotheses : list,
        verbose : bool,
        sigma_interval : Interval,
        name : str,
        tau0 = frac(3),
        plot : bool = False,
        method : int = 1
    ) -> list[Hypothesis]:

    """
    Prove a zero density estimate for zeta given a set of hypotheses, a range of
    values of sigma, and a choice of tau0.

    Paramters
    ---------
    additional_hypotheses : list of Hypothesis
        The list of hypothesis from the literature to assume (other than classical
        results).
    verbose : bool, optional
        If True, results will be logged to console.
    sigma_interval : Interval
        The range of sigma values to consider.
    name : str
        The name of this density estimate, for plotting/logging purposes.
    tau0 : Number or Affine, optional
        The tau_0 value to use. If method = 1 (i.e. using Corollary 11.7), then
        tau_0 can be any sufficiently large number (default is 3). If method = 2
        (i.e. using Corollary 11.8), tau_0 must be an affine function of sigma,
        represented as an Affine object.
    plot : bool, optional
        If True, a graph of the zero-density estimate will be plotted (default
        is False).
    method : int, optional
        The method of converting large value estimates to zero density theorems.
        If method = 1, then Corollary 11.7 (zd.lv_zlv_to_zd(...)) is used, and if
        method = 2, then Corollary 11.8 (zd.lv_zlv_to_zd2(...)) is used. Note
        this parameter affects the expected type of tau0 parameter.

    Returns
    -------
    list of Hypothesis
        A list of zero density estimates, each valid for a range of sigma.
    """

    # Given additional_hypotheses, a list of new hypothesis (other than classical results),
    # find the best density estimate as a piecewise function, then if 'verbose' is true
    # displays the proof of the piece containing the midpoint of sigma_interval.
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(lv.large_value_estimate_L2)

    # Only add raise to power hypothesis if proof method = 1
    if method == 1:
        for k in range(2, 6):
            hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))
    hypotheses.add_hypotheses(additional_hypotheses)

    if method == 1:
        zdes = zd.lv_zlv_to_zd(hypotheses, sigma_interval, tau0)
    elif method == 2:
        zdes = zd.lv_zlv_to_zd2(hypotheses, sigma_interval, tau0)
    else:
        raise NotImplementedError(f"Unknown proof method: {method}")

    if verbose and len(zdes) > 0:
        sigma = sigma_interval.midpoint()
        hyp = next((h for h in zdes if h.data.interval.contains(sigma)), None)
        if hyp is not None:
            print()
            print(f'Found proof of {name}\'s zero-density estimate')
            hyp.recursively_list_proofs()

    if plot and len(zdes) > 0:
        xs = np.linspace(0.5, 1, 100)
        ys = []
        for x in xs:
            zs = [z for z in zdes if z.data.interval.contains(x)]
            if len(zs) > 0:
                ys.append(min(z.data.at(x) / (1 - x) for z in zs))
            else:
                ys.append(None)
        plt.plot(xs, ys)
        plt.title(name + " zero density estimate")
        plt.xlabel("sigma")
        plt.ylabel("A(sigma)")

    return zdes

def prove_zero_density_ingham_1940(verbose=True):
    """
    Prove Ingham's zero density estimate A(σ) <= 3/(2 - σ) using method 1
    """
    return prove_zero_density([], verbose, Interval(frac(1,2), frac(3,4)), 'Ingham')

def prove_zero_density_ingham_1940_v2(verbose=True):
    """
    Prove Ingham's zero density estimate A(σ) < 3/(2 - σ) using method 2
    """
    sigma = Interval(frac(1,2), frac(3,4))
    tau0 = Affine(-1, 2, sigma)
    return prove_zero_density([], verbose, sigma, 'Ingham', tau0=tau0, method=2)

def prove_zero_density_huxley_1972(verbose=True):
    """
    Prove Huxley's zero density estimate A(σ) <= 3/(3σ - 1) for 3/4 <= σ <= 1
    using method 1
    """
    new_hyps = [
        literature.find_hypothesis(keywords='Huxley large value estimate')
    ]

    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(
        literature.find_hypothesis(data=Bound_mu(frac(1,2), frac(1,6)))
    )
    new_hyps.extend(zlv.mu_to_zlv(hypotheses))

    return prove_zero_density(new_hyps, verbose, Interval(frac(3,4), 1), 'Huxley')

def prove_zero_density_huxley_1972_v2(verbose=True):
    """
    Prove Huxley's zero density estimate A(σ) <= 3/(3σ - 1) for 3/4 <= σ <= 1
    using method 1
    """
    new_hyps = [
        literature.find_hypothesis(keywords='Huxley large value estimate')
    ]

    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(
        literature.find_hypothesis(data=Bound_mu(frac(1,2), frac(1,6)))
    )
    new_hyps.extend(zlv.mu_to_zlv(hypotheses))

    sigma = Interval(frac(3,4), 1)
    tau0 = Affine(3, -1, sigma)
    return prove_zero_density(new_hyps, verbose, sigma, 'Huxley', tau0=tau0, method=2)

def prove_zero_density_jutila_1977(verbose=True):
    """
    Prove Jutila's zero density estimate A(σ) <= 2 for 11/14 <= σ <= 1
    using method 1
    """
    new_hyps = [
        literature.find_hypothesis(keywords='Jutila large value estimate, k = 3')
    ]
    return prove_zero_density(new_hyps, verbose, Interval(frac(11,14), 1), 'Jutila')

def prove_zero_density_jutila_1977_v2(verbose=True):
    """
    Prove Jutila's zero density estimate A(σ) <= 2 for 11/14 <= σ <= 1
    using method 2
    """
    new_hyps = [
        literature.find_hypothesis(keywords='Jutila large value estimate, k = 3')
    ]
    sigma = Interval(frac(11,14), 1)
    tau0 = Affine(0, frac(3,2), sigma)
    return prove_zero_density(new_hyps, verbose, sigma, 'Jutila', tau0=tau0, method=2)

def prove_zero_density_heathbrown_1979a(verbose=True):
    """
    Prove Heath-Browns's zero density estimate A(σ) < 9/(7σ - 1) for 11/14 <= σ <= 1
    using method 1
    """
    new_hyps = [
        literature.find_hypothesis(keywords="Jutila large value estimate, k = 3"),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate",
            keywords="Heath-Brown"
            )
        ]
    return prove_zero_density(new_hyps, verbose, Interval(frac(11,14), 1), 'Heath-Brown')

def prove_zero_density_heathbrown_1979a_v2(verbose=True):
    """
    Prove Heath-Browns's zero density estimate A(σ) < 9/(7σ - 1) for 11/14 <= σ <= 1
    using method 2
    """
    new_hyps = [
        literature.find_hypothesis(keywords="Jutila large value estimate, k = 3"),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate",
            keywords="Heath-Brown"
            )
        ]
    sigma = Interval(frac(11,14), 1)
    tau0 = Affine(frac(7,3), -frac(1,3), sigma)
    return prove_zero_density(new_hyps, verbose, sigma, 'Heath-Brown', tau0=tau0, method=2)

def prove_zero_density_heathbrown_1979b(verbose=True):
    """
    Prove Heath-Browns's zero density estimate A(σ) < max(3/(10σ - 7), 4/(4σ - 1))
    for 20/23 <= σ <= 1 using method 1
    """
    new_hyps = [
        literature.find_hypothesis(keywords="Jutila large value estimate, k = 3"),
        literature.find_hypothesis(keywords="Heath-Brown large value estimate"),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate",
            keywords="Heath-Brown"
            )
        ]
    zdts = []
    zdts.append(prove_zero_density(new_hyps, verbose, Interval(frac(20,23), frac(25,28)), 'part 1/2 of the second Heath-Brown'))
    zdts.append(prove_zero_density(new_hyps, verbose, Interval(frac(25,28), 1), 'part 2/2 of the second Heath-Brown'))
    return zdts

def prove_zero_density_heathbrown_1979b_v2(verbose=True):
    """
    Prove Heath-Browns's zero density estimate A(σ) < max(3/(10σ - 7), 4/(4σ - 1))
    for 20/23 <= σ <= 1 using method 2
    """
    new_hyps = [
        literature.find_hypothesis(keywords="Jutila large value estimate, k = 3"),
        literature.find_hypothesis(keywords="Heath-Brown large value estimate"),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate",
            keywords="Heath-Brown"
            )
        ]
    zdts = []
    sigma = Interval(frac(20,23), frac(25,28))
    tau0 = Affine(10, -7, sigma)
    zdts.append(prove_zero_density(new_hyps, verbose, sigma, 'part 1/2 of the second Heath-Brown', tau0=tau0, method=2))

    sigma = Interval(frac(25,28), 1)
    tau0 = Affine(3, -frac(3,4), sigma)
    zdts.append(prove_zero_density(new_hyps, verbose, sigma, 'part 2/2 of the second Heath-Brown', tau0=tau0, method=2))
    return zdts

def prove_zero_density_ivic_1984():
    """
    Prove Ivic's zero-density estimates e.g.
    A(σ) < 3/(2σ)       3831/4791 <= σ <= 1 (actually, we could do slightly better with better
    choice of exponent pair)
    A(σ) < 9/(7σ -1),   41/53 <= σ <= 1
    A(σ) < 6/(5σ - 1),  13/17 <= σ <= 1
    """
    hs = Hypothesis_Set()
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Exponent pair"))
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Exponent pair transform"))
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Upper bound on beta"))

    hs.add_hypotheses(
        ep.compute_exp_pairs(hs, search_depth=5, prune=True)
    )
    hs.add_hypotheses(ep.exponent_pairs_to_beta_bounds(hs))
    hs.add_hypotheses(ep.compute_best_beta_bounds(hs))
    ephs = ep.beta_bounds_to_exponent_pairs(hs)

    for k in range(2, 20):
        h = zd.ivic_ep_to_zd(ephs, k)
        print(h.data, h.proof)

def prove_zero_density_guth_maynard_2024(verbose=True):
    """
    Prove Guth-Maynards's zero density estimate A(σ) < 15/(5σ + 3)
    """
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate",
            keywords="Guth, Maynard"
            )
        ]
    return prove_zero_density(new_hyps, verbose, Interval(frac(7,10), frac(9,10)), "Guth--Maynard")

def prove_zero_density_guth_maynard_2024_v2(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate",
            keywords="Guth, Maynard"
            )
        ]
    sigma = Interval(frac(7,10), frac(9,10))
    tau0 = Affine(1, frac(3,5), sigma)
    return prove_zero_density(new_hyps, verbose, sigma, "Guth--Maynard", tau0=tau0, method=2)

def prove_zero_density_heathbrown_extended(verbose=True):
    """
    Prove the extended version of Heath-Browns zero density estimate A(σ) < 3/(10σ - 7)
    using method 1
    """
    new_hyps = [
        literature.find_hypothesis(keywords="Jutila large value estimate, k = 3"),
        literature.find_hypothesis(keywords="Heath-Brown large value estimate")
        ]

    # Create a hypothesis representing the (3/40, 31/40) exponent pair
    hs = Hypothesis_Set()
    hs.add_hypothesis(
        ep.derived_exp_pair(
            frac(3,40), frac(31,40),
            'See best_proof_of_exponent_pair(frac(3, 40), frac(31, 40))',
            {})
        )

    # Convert the exponent pair to beta bounds, add the other ZLV assumptions,
    # which will be used to calculate the best zeta large value estimate
    new_hyps.extend(bbeta.exponent_pairs_to_beta_bounds(hs))
    return prove_zero_density(new_hyps, verbose, Interval(frac(20,23), 1), 'extended Heath-Brown', tau0=frac(4), plot=True)

def prove_zero_density_heathbrown_extended_v2(verbose=True):
    """
    Prove the extended version of Heath-Browns zero density estimate A(σ) < 3/(10σ - 7)
    using method 2
    """
    new_hyps = [
        literature.find_hypothesis(keywords="Jutila large value estimate, k = 3"),
        literature.find_hypothesis(keywords="Heath-Brown large value estimate")
        ]

    # Create a hypothesis representing the (3/40, 31/40) exponent pair
    hs = Hypothesis_Set()
    hs.add_hypothesis(
        ep.derived_exp_pair(
            frac(3,40), frac(31,40),
            'See best_proof_of_exponent_pair(frac(3, 40), frac(31, 40))',
            {})
        )

    # Convert the exponent pair to beta bounds, add the other ZLV assumptions,
    # which will be used to calculate the best zeta large value estimate
    new_hyps.extend(bbeta.exponent_pairs_to_beta_bounds(hs))

    sigma = Interval(frac(20,23), 1)
    tau0 = Affine(10, -7, sigma)
    return prove_zero_density(new_hyps, verbose, sigma, 'extended Heath-Brown', tau0=tau0, method=2)

def prove_zero_density_bourgain_improved(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", keywords="Bourgain"
        ),
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", keywords="Jutila, k = 3"
        ),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate", keywords="Heath-Brown"
        )
    ]
    return [
        prove_zero_density(new_hyps, verbose, Interval(frac(3,4), frac(4,5)), "part 1/2 of optimized Bourgain", tau0=frac(3)),
        prove_zero_density(new_hyps, verbose, Interval(frac(4,5), 1), "part 2/2 of optimized Bourgain", tau0=frac(3))
    ]

def prove_zero_density_guth_maynard_improved(verbose=True):
    """
    Prove a improved zero-density estimate using Guth--Maynard
    Proposition 12.1
    """
    new_hyps = [
        literature.find_hypothesis(name="Guth--Maynard large value estimate"),
        literature.find_hypothesis(name="Huxley large value estimate"),
        lv.large_value_estimate_L2
    ]
    for k in range(1, 10):
        new_hyps.append(literature.find_hypothesis(keywords=f"Guth--Maynard large value estimate 2 with k = {k}"))

    zdes = prove_zero_density(new_hyps, verbose, Interval(frac(71,100), frac(8,10)), "Guth--Maynard")
    for zde in zdes:
        print(zde.data)
    return zdes

def compute_best_zero_density():

    """
    Compute the best zero-density estimates from the literature
    """

    hs = Hypothesis_Set()
    hs.add_hypotheses(literature)

    # Add the new zero-density estimates (not part of the literature yet!)
    zd.add_zero_density(hs, "2/(9*x - 6)", Interval("[17/22, 38/49]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "9/(8*(2*x - 1))", Interval("[38/49, 4/5]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "3/(10 * x - 7)", Interval("[701/1000, 1]"), Reference.make("Tao--Trudgian--Yang", 2024))
    hs.add_hypotheses(zd.bourgain_ep_to_zd())

    zd.best_zero_density_estimate(hs, verbose=True)

def prove_all_zero_density_estimates():
    print("Proofs using Corollary 11.8 -------------------------------------------------------")
    prove_zero_density_ingham_1940_v2()
    prove_zero_density_huxley_1972_v2()
    prove_zero_density_jutila_1977_v2()
    prove_zero_density_heathbrown_1979a_v2()
    prove_zero_density_heathbrown_1979b_v2()
    prove_zero_density_guth_maynard_2024_v2()
    prove_zero_density_heathbrown_extended_v2()

    print()
    print("Proofs using Corollary 11.7 -------------------------------------------------------")
    prove_zero_density_ingham_1940()
    prove_zero_density_huxley_1972()
    prove_zero_density_jutila_1977()
    prove_zero_density_heathbrown_1979a()
    prove_zero_density_heathbrown_1979b()
    prove_zero_density_ivic_1984()
    prove_zero_density_guth_maynard_2024()
    prove_zero_density_heathbrown_extended()
    prove_zero_density_bourgain_improved()
    compute_best_zero_density()

#################################################################################################
# Derivations for zero-density energy estimates for the Riemann zeta-function

def prove_heath_brown_energy_estimate():

    # Part 1: \sigma \in [1/2, 3/4]
    hypotheses = Hypothesis_Set()
    for k in range(2, 5):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    # Add classical and literature Large value estimates
    hypotheses.add_hypothesis(lv.large_value_estimate_L2)
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2b"))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))

    # tau_0 as an affine function
    tau0 = Affine(0, 2, Interval(frac(1,2), frac(3,4)))
    LVER_domain = Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-tau0.c, -tau0.m, 1],       # tau >= tau0 = m sigma + c
                [2 * tau0.c, 2 * tau0.m, -1] # tau <= 2 tau0 = 2 m sigma + 2 c
            ])
        )

    # Compute the feasible region for LV*(s, t) as a 3-dimensional polytope
    LV_star_hyp = ad.compute_LV_star(hypotheses, LVER_domain, zeta=False, debug=False)

    # domain representing 2 <= tau <= tau0
    LVER_zeta_domain = Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-2, 0, 1],                     # tau0 >= 2
                [tau0.c, tau0.m, -1],           # tau <= tau0 = m sigma + c
            ])
        )
    # Compute the feasible region for LV_{\zeta}*(s, t) as a 3-dimensional polytope
    LVZ_star_hyp = ad.compute_LV_star(hypotheses, LVER_zeta_domain, zeta=True, debug=False)
    ze.lver_to_energy_bound(LV_star_hyp, LVZ_star_hyp, tau0.domain)


    # Part 2: \sigma \in [3/4, 25/28] ----------------------------------------------------------

    hypotheses = Hypothesis_Set()
    for k in range(2, 5):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    # Add classical and literature Large value estimates
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Huxley large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(hypothesis_type="Zeta large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    # Convert all large value estimates -> large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    # Convert all zeta large value estimates -> zeta large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    # tau_0 as an affine function
    tau0 = Affine(4, -1, Interval(frac(3,4), frac(25,28)))

    LVER_domain = Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-tau0.c, -tau0.m, 1],       # tau >= tau0 = m sigma + c
                [2 * tau0.c, 2 * tau0.m, -1] # tau <= 2 tau0 = 2 m sigma + 2 c
            ])
        )

    # Compute the feasible region for LV*(s, t) as a 3-dimensional
    # polytope for a range of sigma
    LV_star_hyp = ad.compute_LV_star(hypotheses, LVER_domain, zeta=False, debug=False)

    # domain representing 2 <= tau <= tau0
    LVER_zeta_domain = Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-2, 0, 1],                     # tau0 >= 2
                [tau0.c, tau0.m, -1],           # tau <= tau0 = m sigma + c
            ])
        )
    # Compute the feasible region for LV_{\zeta}*(s, t) as a 3-dimensional polytope
    LVZ_star_hyp = ad.compute_LV_star(hypotheses, LVER_zeta_domain, zeta=True, debug=False)
    ze.lver_to_energy_bound(LV_star_hyp, LVZ_star_hyp, tau0.domain)

def prove_improved_heath_brown_energy_estimate():

    # tau_0 as a piecewise affine function
    tau0s = [
        Affine(8, -4, Interval(frac(3,4), frac(4,5)))
    ]

    hypotheses = Hypothesis_Set()

    for k in range(2, 5):
       hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    # Add classical and literature Large value estimates
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Huxley large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    # Convert all large value estimates -> large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))

    # Convert tau_0 into a Region of (sigma, tau)
    # domain representing tau0 <= tau <= 2tau0
    LVER_domain = Region.disjoint_union([
        Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-tau0.c, -tau0.m, 1],       # tau >= tau0 = m sigma + c
                [2 * tau0.c, 2 * tau0.m, -1] # tau <= 2 tau0 = 2 m sigma + 2 c
            ])
        )
        for tau0 in tau0s
    ])
    # Compute the feasible region for LV*(s, t) as a 3-dimensional
    # polytope for a range of sigma
    LV_star_hyp = ad.compute_LV_star(hypotheses, LVER_domain, zeta=False)

    # New set of hypothesis for the zeta LVER computation
    hypotheses = Hypothesis_Set()

    for k in range(2, 3):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Huxley large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(hypothesis_type="Zeta large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    # Convert all large value estimates -> large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    # Convert all zeta large value estimates -> zeta large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    # domain representing 2 <= tau <= tau0
    LVER_zeta_domain = Region.disjoint_union([
        Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-2, 0, 1],                     # tau0 >= 2
                [tau0.c, tau0.m, -1],           # tau <= tau0 = m sigma + c
            ])
        )
        for tau0 in tau0s
    ])
    # Compute the feasible region for LV_{\zeta}*(s, t) as a 3-dimensional polytope
    LVZ_star_hyp = ad.compute_LV_star(hypotheses, LVER_zeta_domain, zeta=True)
    bounds = ze.lver_to_energy_bound(LV_star_hyp, LVZ_star_hyp, Interval(frac(1,2), 1))

def prove_zero_density_energy_2():
    hypotheses = Hypothesis_Set()

    for k in range(2, 5):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    # Add classical and literature Large value estimates
    hypotheses.add_hypothesis(lv.large_value_estimate_L2)
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    # Convert all large value estimates -> large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))

    tau0 = Affine(0, 2, Interval(frac(7,10), frac(3,4)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs:
        print(h.data)
        h.recursively_list_proofs()
    return hs

def prove_zero_density_energy_3():
    hypotheses = Hypothesis_Set()

    for k in range(2, 5):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    # Add classical and literature Large value estimates
    #hypotheses.add_hypothesis(lv.large_value_estimate_L2)
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Huxley large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(hypothesis_type="Zeta large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    # Convert all large value estimates -> large value energy region
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    # tau_0 = 8\sigma - 4
    tau0 = Affine(8, -4, Interval(frac(3,4), frac(4,5)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_4():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 10"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(664,877), frac(31,40)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=True)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_5():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 6"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(42,55), frac(79,103)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_6():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 5"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(79,103), frac(84,109)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_7():
    hypotheses = Hypothesis_Set()

    for k in range(2, 6):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Bourgain optimized large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 5"))
    hypotheses.add_hypothesis(literature.find_hypothesis(hypothesis_type="Zeta large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(8, -4, Interval(frac(84,109), frac(5,6)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_8():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 13"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(37,49), frac(443,586)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=True)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_9():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 12"))
    hypotheses.add_hypothesis(literature.find_hypothesis(hypothesis_type="Zeta large value estimate"))

    # Add Heath-Brown estimates
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(443,586), frac(373,493)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_10():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 11"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(373,493), frac(103,136)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_11():
    hypotheses = Hypothesis_Set()

    for k in range(2, 4):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 11"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 12"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 13"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 14"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 15"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 16"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 17"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 18"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Jutila large value estimate with k = 19"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(3,4), frac(103,136)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=False)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_12():
    hypotheses = Hypothesis_Set()

    for k in range(2, 5):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    sigma_interval = Interval(frac(7,10), frac(42,55))
    tau_0 = 2

    # To speed things up, we will combine LV regions separately first instead of
    # converting them into LVER
    lv_hyps = Hypothesis_Set()
    lv_hyps.add_hypothesis(lv.large_value_estimate_L2)
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 1"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 2"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 3"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 4"))
    #lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 5"))

    # Add the combined large value estimates as a single Hypothesis
    combined_lv = lv.combine_large_value_estimates(
        lv_hyps,
        domain = Region.from_polytope(
            Polytope.rect(
                (sigma_interval.x0, sigma_interval.x1),
                (0, 2 * tau_0), # even with raise-to-power hypotheses, we will only ever need estimates up to 2tau_0
                (0, Constants.LV_DEFAULT_UPPER_BOUND)
            )
        ),
        simplify = True,
        verbose = True
    )
    print(combined_lv.data.region)

    hypotheses.add_hypothesis(combined_lv)
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, tau_0, sigma_interval)
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=True)
    for h in hs: print(h.data)
    return hs

def prove_zero_density_energy_13():
    hypotheses = Hypothesis_Set()

    for k in range(2, 5):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    # To speed things up, we will combine LV regions separately first instead of
    # converting them into LVER
    lv_hyps = Hypothesis_Set()
    lv_hyps.add_hypothesis(lv.large_value_estimate_L2)
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 1"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 2"))
    lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 3"))
    #lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 4"))
    #lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 5"))
    #lv_hyps.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value estimate 2 with k = 6"))

    # Add the combined large value estimates as a single Hypothesis
    combined_lv = lv.combine_large_value_estimates(lv_hyps, simplify=True, verbose=True)
    print(combined_lv.data)
    hypotheses.add_hypothesis(combined_lv)
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))

    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=False))
    hypotheses.add_hypotheses(ad.lv_to_lver(hypotheses, zeta=True))

    tau0 = Affine(0, 2, Interval(frac(7,10), frac(3,4)))
    hs = ze.lver_to_energy_bound(hypotheses, tau0, debug=True)
    for h in hs: print(h.data)
    return hs

def prove_all_zero_density_energy_estimates():
    prove_heath_brown_energy_estimate()
    prove_improved_heath_brown_energy_estimate()
    prove_zero_density_energy_2()
    prove_zero_density_energy_3()
    prove_zero_density_energy_4()
    prove_zero_density_energy_5()
    prove_zero_density_energy_6()
    prove_zero_density_energy_7()
    prove_zero_density_energy_8()
    prove_zero_density_energy_9()
    prove_zero_density_energy_10()

    prove_zero_density_energy_12()

#################################################################################################
# Derivations for prime gap theorems

def prove_prime_gap2():
    hs = Hypothesis_Set()

    # Add zero-density estimates and energy theorems from the literature
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zero density estimate"))
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zero density energy estimate"))

    # New zero-density estimates
    ref = Reference.make("Tao--Trudgian--Yang", 2024)
    zd.add_zero_density(hs, "2/(9*x - 6)", Interval("[17/22, 38/49]"), ref)
    zd.add_zero_density(hs, "9/(8*(2*x - 1))", Interval("[38/49, 4/5]"), ref)
    zd.add_zero_density(hs, "3/(10 * x - 7)", Interval("[701/1000, 1]"), ref)
    hs.add_hypotheses(zd.bourgain_ep_to_zd())
    zd.add_zero_density(hs, "3/(40 * x - 35)", Interval("[39/40, 40/41)"), ref)
    zd.add_zero_density(hs, "2/(13 * x - 10)", Interval("[40/41, 41/42)"), ref)

    # Set of new additive energy estimates
    hs.add_hypotheses([
        ze.literature_zero_density_energy_estimate("5 * (18 - 19 * x) / ((2 * (5 * x + 3)) * (1 - x))", Interval(frac(7,10), 0.7255782330963900973348270455), ref),
        ze.literature_zero_density_energy_estimate("2 * (45 - 44 * x) / ((2 * x + 15) * (1 - x))", Interval(0.7255782330963900973348270455, frac(3,4)), ref),
        ze.literature_zero_density_energy_estimate("(197 - 220 * x) / (8 * (5 * x - 1) * (1 - x))", Interval(frac(3,4), frac(289,380)), ref),
        ze.literature_zero_density_energy_estimate("3 * (29 - 30 * x) / (5 * (5 * x - 1) * (1 - x))", Interval(frac(289,380), 0.7929182893891673924914902646), ref),
        ze.literature_zero_density_energy_estimate("(40 - 36 * x) / ((20 * x - 5) * (1 - x))", Interval(0.7929182893891673924914902646, frac(5,6)), ref)
    ])

    # Compute \theta_{gap, 2}
    pg.compute_gap2(hs, debug=False)


def compute_prime_excep():
    hs = Hypothesis_Set()

    # Add zero-density estimates and energy theorems from the literature
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zero density estimate"))
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zero density energy estimate"))

    # New zero-density estimates
    ref = Reference.make("Tao--Trudgian--Yang", 2024)
    zd.add_zero_density(hs, "2/(9*x - 6)", Interval("[17/22, 38/49]"), ref)
    zd.add_zero_density(hs, "9/(8*(2*x - 1))", Interval("[38/49, 4/5]"), ref)
    zd.add_zero_density(hs, "3/(10 * x - 7)", Interval("[701/1000, 1]"), ref)
    hs.add_hypotheses(zd.bourgain_ep_to_zd())
    zd.add_zero_density(hs, "3/(40 * x - 35)", Interval("[39/40, 40/41)"), ref)
    zd.add_zero_density(hs, "2/(13 * x - 10)", Interval("[40/41, 41/42)"), ref)

    # Set of new additive energy estimates
    hs.add_hypotheses([
        ze.literature_zero_density_energy_estimate("5 * (18 - 19 * x) / ((2 * (5 * x + 3)) * (1 - x))", Interval(frac(7,10), 0.7255782330963900973348270455), ref),
        ze.literature_zero_density_energy_estimate("2 * (45 - 44 * x) / ((2 * x + 15) * (1 - x))", Interval(0.7255782330963900973348270455, frac(3,4)), ref),
        ze.literature_zero_density_energy_estimate("(197 - 220 * x) / (8 * (5 * x - 1) * (1 - x))", Interval(frac(3,4), frac(289,380)), ref),
        ze.literature_zero_density_energy_estimate("3 * (29 - 30 * x) / (5 * (5 * x - 1) * (1 - x))", Interval(frac(289,380), 0.7929182893891673924914902646), ref),
        ze.literature_zero_density_energy_estimate("(40 - 36 * x) / ((20 * x - 5) * (1 - x))", Interval(0.7929182893891673924914902646, frac(5,6)), ref)
    ])

    # Compute prime gap exception at theta
    return pg.prime_excep(hs)


#################################################################################################

def prove_all():
    # van_der_corput_pair(10)
    # prove_hardy_littlewood_mu_bound()
    prove_exponent_pairs()
    # prove_all_large_value_estimates()
    # prove_all_zero_density_estimates()
    # prove_all_zero_density_energy_estimates()
    # prove_prime_gap2()
