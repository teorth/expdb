# This file contains all bounds that can be derived from other bounds, including those from literature.py

from literature import *
import exponent_pair as ep
import zero_density_estimate as zd


# Establish the classical van der Corput exponent pair (\frac{1}{2^k-2}, 1 - \frac{k-1}{2^k-2})
def van_der_corput_pair(k):
    if k < 2:
        raise ValueError("k must be at least 2.")
    exp_pair = B_transform(trivial_exp_pair)
    for _ in range(k - 2):
        exp_pair = A_transform(exp_pair)
    print(f"The van der Corput pair for k = {k} is {exp_pair.desc()}")
    return exp_pair


# Prove the Hardy-Littlewood bound mu(1/2) \leq 1/6 using the van der Corput pair (1/6, 2/3).
def prove_hardy_littlewood_mu_bound():
    HL_bound = literature.find_hypothesis(data=Bound_mu(frac(1, 2), frac(1, 6)))
    print(f"We will reprove {HL_bound.desc()}.")
    B_exp_pair = B_transform(trivial_exp_pair)
    print(f"We have {B_exp_pair.desc_with_proof()}")
    AB_exp_pair = A_transform(B_exp_pair)
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
    lv.optimize_bourgain_large_value_estimate()



######################################################################################


# Given additional_hypotheses, a list of new hypothesis (other than classical results),
# find the best density estimate as a piecewise function, then if 'verbose' is true
# displays the proof of the piece containing 'sigma'. 
def prove_zero_density(additional_hypotheses, verbose, sigma, name, tau0=frac(3)):
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(lv.large_value_estimate_L2)
    for k in range(2, 10):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))
    hypotheses.add_hypotheses(additional_hypotheses)
    
    zdes = zd.lv_zlv_to_zd(hypotheses, Interval(frac(1,2), 1))
    
    if verbose and len(zdes) > 0:
        # for h in zdes:
        #     print(h.data)
        hyp = next((h for h in zdes if h.data.interval.contains(sigma)), None)
        if hyp is not None:
            print()
            print(f'Found proof of {name}\'s zero-density estimate')
            hyp.recursively_list_proofs()
    return zdes
    
# Prove Ingham's zero density estimate A(s) < 3/(1-s)
def prove_ingham_zero_density(verbose=True):
    return prove_zero_density([], verbose, frac(1,2), 'Ingham')

# Prove Huxley's zero density estimate A(s) < 3/(3s - 1)
def prove_huxley_zero_density(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type='Large value estimate',
            keywords='Huxley'
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(7,8), 'Huxley')

# Prove Jutila's proof of the density hypothesis for s > 11/14.
def prove_jutila_zero_density(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type='Large value estimate', 
            keywords='Jutila, k = 3'
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(9,10), 'Jutila')

# Prove Heath-Browns's zero density estimate A(s) < 9/(7s - 1)
def prove_heathbrown_zero_density(verbose=True):
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
def prove_heathbrown_zero_density2(verbose=True):
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

# Prove Guth-Maynards's zero density estimate A(s) < 15/(5 + 3s)
def prove_guth_maynard_zero_density(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", 
            keywords="Guth, Maynard"
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(3,4), "Guth--Maynard")

# Prove the extended version of Heath-Browns zero density estimate A(s) < 3/(10s - 7)
def prove_extended_heathbrown_zero_density(verbose=True):
    
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
    
    for h in new_hyps:
        h.recursively_list_proofs()
    return prove_zero_density(new_hyps, verbose, frac(9,10), 'Heath-Brown', tau0=frac(5))

# Prove Ivi\'{c}'s zero-density estimates 
# A(s) < 3/(2s)  3831/4791 <= s <= 1 (actually, we could do slightly better with better 
# choice of exponent pair)
# A(s) < 9/(7s -1), 41/53 <= s <= 1
# A(s) < 6/(5s - 1),  13/17 <= s <= 1
def prove_ivic_zero_density():

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

# Compute the best zero-density estimates from the literature
def compute_best_zero_density():
    hs = Hypothesis_Set()
    hs.add_hypotheses(literature)

    # Add the new zero-density estimates (not part of the literature yet!)
    zd.add_zero_density(hs, "2/(9*x - 6)", Interval("[17/22, 38/49]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "9/(8*(2*x - 1))", Interval("[38/49, 4/5]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "3/(10 * x - 7)", Interval("[701/1000, 1]"), Reference.make("Tao--Trudgian--Yang", 2024))

    zd.best_zero_density_estimate(hs, verbose=True)

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
    # prove_ingham_zero_density()
    # prove_huxley_zero_density()
    # prove_jutila_zero_density()
    # prove_heathbrown_zero_density()
    # prove_heathbrown_zero_density2()
    # prove_guth_maynard_zero_density()
    # prove_extended_heathbrown_zero_density()
    # prove_ivic_zero_density()
    compute_best_zero_density()

def prove_all():
    # prove_hardy_littlewood_mu_bound()
    prove_exponent_pairs()
    # prove_zero_density_estimates()
    # prove_bourgain_large_values_theorem()

prove_all()
