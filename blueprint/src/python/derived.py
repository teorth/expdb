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
def best_proof_of_exponent_pair(k, l):
    hyp = ep.find_best_proof(
        frac(3, 40), frac(31, 40), literature, Proof_Optimization_Method.DATE
    )
    if hyp is not None:
        print()
        print(f'Found proof of ({k}, {l}) with complexity = {hyp.proof_complexity()} and date = {hyp.proof_date()}:')
        hyp.recursively_list_proofs()
    else:
        print('Failed to prove the exponent pair ({k}, {l}).')
        
    hyp = ep.find_best_proof(
        frac(3, 40), frac(31, 40), literature, Proof_Optimization_Method.COMPLEXITY
    )
    if hyp is not None:
        print()
        print(f'Found proof of ({k}, {l}) with complexity = {hyp.proof_complexity()} and date = {hyp.proof_date()}:')
        hyp.recursively_list_proofs()
    else:
        print('Failed to prove the exponent pair ({k}, {l}).')


# Given additional_hypotheses, a list of new hypothesis (other than classical results),
# find the best density estimate as a piecewise function, then if 'verbose' is true
# displays the proof of the piece containing 'sigma'. 
def prove_zero_density(additional_hypotheses, verbose, sigma, name):
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(lv.large_value_estimate_L2)
    for k in range(2, 10):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))
    hypotheses.add_hypotheses(additional_hypotheses)
    
    zdes = zd.lv_zlv_to_zd(hypotheses, Interval(frac(1,2), 1))
    
    if verbose and len(zdes) > 0:
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

# Prove Guth-Maynards's zero density estimate A(s) < 15/(5 + 3s)
def prove_guth_maynard_zero_density(verbose=True):
    new_hyps = [
        literature.find_hypothesis(
            hypothesis_type='Large value estimate', 
            keywords='Guth, Maynard'
            )
        ]
    return prove_zero_density(new_hyps, verbose, frac(3,4), 'Guth--Maynard')

def prove_all():
    # prove_hardy_littlewood_mu_bound()
    best_proof_of_exponent_pair(frac(3, 40), frac(31, 40))
    prove_exponent_pair(frac(1101653,15854002), frac(12327829,15854002))
    # prove_heathbrown_exponent_pairs()
    prove_ingham_zero_density()
    prove_huxley_zero_density()
    prove_jutila_zero_density()
    prove_guth_maynard_zero_density()

prove_all()
