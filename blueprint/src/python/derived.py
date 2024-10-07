# This file contains all bounds that can be derived from other bounds, including those from literature.py

from literature import *
import exponent_pair as ep
import zero_density_estimate as zd

import time


# Establish the classical van der Corput exponent pair (\frac{1}{2^k-2}, 1 - \frac{k-1}{2^k-2})
def van_der_corput_pair(k):
    if k < 2:
        raise ValueError("k must be at least 2.")
    A_transform = literature.find_hypothesis(keywords="van der Corput A transform")
    B_transform = literature.find_hypothesis(keywords="van der Corput B transform")
    exp_pair = B_transform.data.transform(trivial_exp_pair)
    for _ in range(k - 2):
        exp_pair = A_transform.data.transform(exp_pair)
    print(f"The van der Corput pair for k = {k} is {exp_pair.desc()}")
    return exp_pair


# Prove the Hardy-Littlewood bound mu(1/2) \leq 1/6 using the van der Corput pair (1/6, 2/3).
def prove_hardy_littlewood_mu_bound():
    HL_bound = literature.find_hypothesis(data=Bound_mu(frac(1, 2), frac(1, 6)))
    A_transform = literature.find_hypothesis(keywords="van der Corput A transform")
    B_transform = literature.find_hypothesis(keywords="van der Corput B transform")
    print(f"We will reprove {HL_bound.desc()}.")
    B_exp_pair = B_transform.data.transform(trivial_exp_pair)
    print(f"We have {B_exp_pair.desc_with_proof()}")
    AB_exp_pair = A_transform.data.transform(B_exp_pair)
    print(f"This implies {AB_exp_pair.desc_with_proof()}")
    mu_bound = obtain_mu_bound_from_exponent_pair(AB_exp_pair)
    print(f"This implies {mu_bound.desc_with_proof()}")
    return mu_bound


def prove_exponent_pair(k, l):

    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on beta")
    )

    # Include all literature exponent pairs and expand the hull using exponent
    # pair transforms
    hypotheses.add_hypothesis(trivial_exp_pair)
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair transform")
    )
    hypotheses.add_hypotheses(compute_exp_pairs(hypotheses, search_depth=1))

    # Perform 1 iteration of exponent pair -> beta bounds -> exponent pair
    hypotheses.add_hypotheses(exponent_pairs_to_beta_bounds(hypotheses))
    hypotheses.add_hypotheses(compute_best_beta_bounds(hypotheses))
    new_exp_pairs = beta_bounds_to_exponent_pairs(hypotheses)

    eph = next((h for h in new_exp_pairs if h.data.k == k and h.data.l == l), None)
    if eph is not None:
        print()
        print(f"Proof of the exponent pair ({k}, {l}) exponent pair:")
        eph.recursively_list_proofs()
    else:
        print('Failed to prove the exponent pair ({k}, {l}).')

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

# Find the shortest proof of the exponent pair (k, l)
def best_proof_of_exponent_pair(k, l, proof_method=Proof_Optimization_Method.DATE, verbose=True):
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


######################################################################################

def prove_bourgain_large_values_theorem():
    
    start_time = time.time()
    lv_hyps = lv.optimize_bourgain_large_value_estimate()
    print("Computed in", time.time() - start_time, "s")
    pieces = []
    for lvh in lv_hyps:
        pieces.extend(lvh.data.bound.pieces)
    bound = Piecewise(pieces)
    bound.plot_domain((1/2, 1), (1, 3))
    
    for lvh in lv_hyps:
        print(lvh.data, lvh.proof)

# Prove the Guth--Maynard large values theorem (Theorem 10.26) using the Large Value Energy 
# Region estimates due to Guth--Maynard
def prove_guth_maynard_large_values_theorem():
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value energy region 1 with k = 2"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value energy region 2"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Guth--Maynard large value energy region 3"))
    
    # Compute the feasible region of (\sigma, \tau, \rho) values 
    region = ad.lver_to_lv(hypotheses)
    
    # Take \tau = 6/5 (TODO: replace this step with Huxley subdivision once it is implemented)
    region = region.substitute({1: frac(6,5)})

    # Constrain to the range 7/10 \leq \sigma \leq 8/10 (\rho unconstrained)
    domain = Polytope([
        [-frac(7,10), 1, 0], [frac(8,10), -1, 0]
    ])
    region.child = [Region.from_polytope(r.child.intersect(domain)) for r in region.child]
    region.child = [r for r in region.child if not r.child.is_empty(include_boundary=False)]

    # Simplify the region
    poly = Polytope.try_union([r.child for r in region.child])

    print("Proved feasible region for (σ, τ, ρ):", poly.to_str("στρ"))

######################################################################################
# Derivations of zero-density estimates for the Riemann zeta-function

# Given additional_hypotheses, a list of new hypothesis (other than classical results),
# find the best density estimate as a piecewise function, then if 'verbose' is true
# displays the proof of the piece containing 'sigma'. 
def prove_zero_density(additional_hypotheses, verbose, sigma, name, tau0=frac(3), plot=False):
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(lv.large_value_estimate_L2)
    for k in range(2, 10):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))
    hypotheses.add_hypotheses(additional_hypotheses)
    
    zdes = zd.lv_zlv_to_zd(hypotheses, Interval(frac(1,2), 1), tau0)
    
    if verbose and len(zdes) > 0:
        # for h in zdes:
        #     print(h.data)
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
    
# Prove Ingham's zero density estimate A(s) < 3/(1-s)
def prove_zero_density_ingham_1940(verbose=True):
    return prove_zero_density([], verbose, frac(1,2), 'Ingham')

# Prove Huxley's zero density estimate A(s) < 3/(3s - 1)
def prove_zero_density_huxley_1972(verbose=True):
    
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type='Large value estimate',
            keywords='Huxley'
            )
        ]
    
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(
        literature.find_hypothesis(data=Bound_mu(frac(1,2), frac(1,6)))
    )
    new_hyps.extend(zlv.mu_to_zlv(hypotheses))
    
    return prove_zero_density(new_hyps, verbose, frac(7,8), 'Huxley')

# Prove Jutila's proof of the density hypothesis for s > 11/14.
def prove_zero_density_jutila_1977(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type='Large value estimate', 
            keywords='Jutila, k = 3'
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(9,10), 'Jutila')

# Prove Heath-Browns's zero density estimate A(s) < 9/(7s - 1)
def prove_zero_density_heathbrown_1979(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Jutila, k = 3"
            ),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate", 
            keywords="Heath-Brown"
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(9,10), 'Heath-Brown')

# Prove Heath-Browns's second zero density estimate A(s) < max(3/(10s - 7), 4/(4s - 1))
def prove_zero_density_heathbrown_1979_2(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Jutila, k = 3"
            ),
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Heath-Brown"
            ),
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate", 
            keywords="Heath-Brown"
            )
        ]
    zdts = []
    zdts.append(prove_zero_density(new_hyps, verbose, frac(20,23), 'part 1/2 of the second Heath-Brown'))
    zdts.append(prove_zero_density(new_hyps, verbose, frac(22,23), 'part 2/2 of the second Heath-Brown'))
    return zdts

# Prove Ivi\'{c}'s zero-density estimates 
# A(s) < 3/(2s)  3831/4791 <= s <= 1 (actually, we could do slightly better with better 
# choice of exponent pair)
# A(s) < 9/(7s -1), 41/53 <= s <= 1
# A(s) < 6/(5s - 1),  13/17 <= s <= 1
def prove_zero_density_ivic_1984():

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
        # h.recursively_list_proofs()

# Prove Guth-Maynards's zero density estimate A(s) < 15/(5 + 3s)
def prove_zero_density_guth_maynard_2024(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Guth, Maynard"
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(3,4), "Guth--Maynard")

# Prove the extended version of Heath-Browns zero density estimate A(s) < 3/(10s - 7)
def prove_zero_density_heathbrown_extended(verbose=True):
    
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Jutila, k = 3"
            ),
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Heath-Brown"
            )
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
    
    # TODO: this proof currently has a missing range \sigma \in [37/40, 1], 
    # which requires Huxley subdivision to prove. We could also prove this by 
    # incorporating more bounds
    return prove_zero_density(new_hyps, verbose, frac(9,10), 'Heath-Brown', tau0=frac(5), plot=True)

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
        #prove_zero_density(new_hyps, verbose, frac(22,29), "part 1/2 of optimized Bourgain", tau0=frac(3)),
        prove_zero_density(new_hyps, verbose, frac(25,32), "part 2/2 of optimized Bourgain", tau0=frac(3))
    ]

# Compute the best zero-density estimates from the literature
def compute_best_zero_density():
    hs = Hypothesis_Set()
    hs.add_hypotheses(literature)

    # Add the new zero-density estimates (not part of the literature yet!)
    zd.add_zero_density(hs, "2/(9*x - 6)", Interval("[17/22, 38/49]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "9/(8*(2*x - 1))", Interval("[38/49, 4/5]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "3/(10 * x - 7)", Interval("[701/1000, 1]"), Reference.make("Tao--Trudgian--Yang", 2024))
    hs.add_hypotheses(zd.bourgain_ep_to_zd())
    # New Pintz-type estimates 
    zd.add_zero_density(hs, "3/(40 * x - 35)", Interval("[39/40, 40/41)"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "2/(13 * x - 10)", Interval("[40/41, 41/42)"), Reference.make("Tao--Trudgian--Yang", 2024))
    
    zd.best_zero_density_estimate(hs, verbose=True)

#################################################################################################
# Derivations for zero-density energy estimates for the Riemann zeta-function

def prove_improved_heath_brown_energy_estimate():
    
    # tau_0 as a piecewise affine function 
    tau0s = [
        # For tracing the bound (18 - 19s) / (6s - 2)
        # Choose tau_0 = 6s - 2 for now
        #Affine(8, -4, Interval(frac(3,4), frac(65,86)))
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
    LV_star_hyp.desc_with_proof()
    
    # New set of hypothesis for the zeta LVER computation
    hypotheses = Hypothesis_Set()

    for k in range(2, 3):
        hypotheses.add_hypothesis(ad.get_raise_to_power_hypothesis(k))

    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Huxley large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(hypothesis_type="Zeta large value estimate"))
    hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2a"))
    #hypotheses.add_hypothesis(literature.find_hypothesis(keywords="Heath-Brown large value energy region 2b"))
    
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
    LVZ_star_hyp.desc_with_proof()

    bounds = ze.compute_best_energy_bound(LV_star_hyp, LVZ_star_hyp, Interval(frac(1,2), 1))

# Prove a zero density energy estimates using Guth--Maynard's LV theorems
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

    for h in hypotheses:
        print(h)
        
    # tau_0 as a piecewise affine function 
    tau0s = [
        Affine(0, 2, Interval(frac(7,10), frac(3,4)))
    ]

    # For each interval of tau_0
    for tau0 in tau0s:
        sigma_interval = tau0.domain

        # domain representing tau0 <= tau <= 2 tau0
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
        LV_star_hyp.desc_with_proof()

        # domain representing 2 <= tau <= tau0
        LVER_zeta_domain = Region.disjoint_union([
            Region.from_polytope(
                Polytope([
                    [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                    [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                    [-2, 0, 1],                  # tau0 >= 2
                    [tau0.c, tau0.m, -1],        # tau <= tau0 = m sigma + c
                ])
            )
            for tau0 in tau0s
        ])
        
        # Compute the feasible region for LV_{\zeta}*(s, t) as a 3-dimensional polytope
        LVZ_star_hyp = ad.compute_LV_star(hypotheses, LVER_zeta_domain, zeta=True)
        if LVZ_star_hyp is not None:
            LVZ_star_hyp.desc_with_proof()

        bounds = ze.compute_best_energy_bound(LV_star_hyp, LVZ_star_hyp, sigma_interval)

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
    sigma_interval = tau0.domain

    # domain representing tau0 <= tau <= 2 tau0
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
    LV_star_hyp = ad.compute_LV_star(hypotheses, LVER_domain, zeta=False)
    LV_star_hyp.desc_with_proof()

    # domain representing 2 <= tau <= tau0
    LVER_zeta_domain = Region.from_polytope(
            Polytope([
                [-tau0.domain.x0, 1, 0],     # sigma >= sigma_interval.x0
                [tau0.domain.x1, -1, 0],     # sigma <= sigma_interval.x1
                [-2, 0, 1],                  # tau0 >= 2
                [tau0.c, tau0.m, -1],        # tau <= tau0 = m sigma + c
            ])
        )
        
    # Compute the feasible region for LV_{\zeta}*(s, t) as a 3-dimensional polytope
    LVZ_star_hyp = ad.compute_LV_star(hypotheses, LVER_zeta_domain, zeta=True)
    if LVZ_star_hyp is not None:
        LVZ_star_hyp.desc_with_proof()

    bounds = ze.compute_best_energy_bound(LV_star_hyp, LVZ_star_hyp, sigma_interval)

#################################################################################################

def prove_exponent_pairs():
    # prove_heathbrown_exponent_pairs()
    prove_exponent_pair(frac(1101653,15854002), frac(12327829,15854002))
    prove_exponent_pair(frac(1959,47230), frac(3975,4723))
    prove_exponent_pair(frac(1175779,38456886), frac(16690288,19228443))
    prove_exponent_pair(frac(89,3478), frac(15327,17390))

    best_proof_of_exponent_pair(frac(1, 6), frac(2, 3))
    best_proof_of_exponent_pair(frac(13, 31), frac(16, 31))
    best_proof_of_exponent_pair(frac(4, 11), frac(6, 11))
    best_proof_of_exponent_pair(frac(2, 7), frac(4, 7))
    best_proof_of_exponent_pair(frac(5, 24), frac(15, 24))
    best_proof_of_exponent_pair(frac(4, 18), frac(11, 18))
    best_proof_of_exponent_pair(frac(3, 40), frac(31, 40), Proof_Optimization_Method.DATE)
    #best_proof_of_exponent_pair(frac(3, 40), frac(31, 40), Proof_Optimization_Method.COMPLEXITY)

def prove_zero_density_estimates():
    # prove_zero_density_ingham_1940()
    # prove_zero_density_huxley_1972()
    # prove_zero_density_jutila_1977()
    # prove_zero_density_heathbrown_1979()
    # prove_zero_density_heathbrown_1979_2()
    # prove_zero_density_ivic_1984()
    # prove_zero_density_guth_maynard_2024()
    # prove_zero_density_heathbrown_extended()
    prove_zero_density_bourgain_improved()
    # compute_best_zero_density()

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
    ze.compute_best_energy_bound(LV_star_hyp, LVZ_star_hyp, tau0.domain)


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
    ze.compute_best_energy_bound(LV_star_hyp, LVZ_star_hyp, tau0.domain)

def prove_all():
    # van_der_corput_pair(10)
    # prove_hardy_littlewood_mu_bound()
    # prove_exponent_pairs()
    # prove_bourgain_large_values_theorem()
    # prove_zero_density_estimates()
    prove_heath_brown_energy_estimate()
    # prove_improved_heath_brown_energy_estimate()
    # prove_guth_maynard_large_values_theorem()
    # prove_zero_density_energy_2()
    # prove_zero_density_energy_3()

prove_all()
