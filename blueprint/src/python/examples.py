from literature import *
import large_values as lv
import zero_density_estimate as zd
import zeta_large_values as zlv

# Temporary debugging functionality
import time


# Example code
def mu_bound_examples():
    print("mu bound examples --------------------------------------------------")

    hypotheses = Hypothesis_Set()
    print(f"Assume only the trivial bound.")
    print(best_mu_bound(frac(3, 4), hypotheses).desc_with_proof())

    print()
    HL = literature.find_hypothesis(keywords="Hardy, Littlewood, bound on \\mu(1/2)")
    print(f"In 1913, Hardy and Littlewood proved {HL.desc_with_proof()}")
    bound = Bound_mu(frac(11, 15), frac(1, 15))
    print(
        f"The bound {bound} was obtained in {literature.find_hypothesis(data=bound).desc_with_proof()}"
    )

    print()
    print(f"Now assume the result of Hardy and Littlewood.")
    hypotheses.add_hypothesis(HL)
    print(best_mu_bound(frac(3, 4), hypotheses).desc_with_proof())

    print()
    print(f"Now assume all unconditional bounds known up to 1985")
    hypotheses.add_hypotheses(literature.list_hypotheses("Upper bound on mu", 1985))
    print(best_mu_bound(frac(3, 4), hypotheses).desc_with_proof())

    print()
    print(f"Now assume all unconditional bounds.")
    hypotheses.add_hypotheses(literature.list_hypotheses("Upper bound on mu"))
    start_time = time.time()
    prove_mu_bound(frac(3, 4), frac(13, 168), hypotheses)
    print(f"Computed in {time.time() - start_time} sec")

    print()
    print(f"The computed convex hull is cached, so this second call should be faster.")
    print(f"Once again assume all unconditional bounds.")
    start_time = time.time()
    prove_mu_bound(frac(4, 7), frac(8, 63), hypotheses)
    print(f"Computed in {time.time() - start_time} sec")

    print()
    print(
        f"Now assume all unconditional exponent pairs, and all unconditional exponent pair transforms."
    )
    add_exp_pairs_all(hypotheses)
    hypotheses.add_hypotheses(literature.find_hypothesis(hypothesis_type="Exponent pair transform"))
    ans = prove_mu_bound(frac(3, 4), frac(8, 63), hypotheses)
    print(f"Recursive dependencies used to prove '{ans}':")
    ans.recursive_dependencies().list_proofs()

    print()
    print(
        "Under the same assumptions, the best piecewise-linear bound on mu in [1/3,9/10) is"
    )
    mbs = best_mu_bound_piecewise(frac(1, 3), frac(9, 10), hypotheses)
    for b in mbs:
        print(f"\\mu(x) \\leq {b}")

    print()
    print(f"Now assume the Lindelof hypothesis.")
    hypotheses.add_hypothesis(Lindelof_hypothesis)
    print(best_mu_bound(frac(3, 4), hypotheses).desc_with_proof())

    print(best_mu_bound(0, hypotheses).desc_with_proof())

    print(best_mu_bound(2, hypotheses).desc_with_proof())

    prove_mu_bound(frac(1, 4), frac(1, 4), hypotheses)

    # Currently this is printing a long list of hypothesis, so I have commented it out for now
    # print()
    # print(f'The list of ambient upper bound hypotheses is as follows:')
    # hypotheses.list_proofs()


def exp_pair_examples():
    hypotheses = Hypothesis_Set()
    print(
        f"Assume all known exponent pairs as well as the van der Corput A/B transforms"
    )
    hypotheses.add_hypothesis(ep.trivial_exp_pair)
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair transform")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair to beta bound transform")
    )
    hypotheses.add_hypotheses(
        ep.compute_exp_pairs(hypotheses, search_depth=5, prune=True)
    )
    hyps = ep.compute_convex_hull(hypotheses)

    print(
        f"The computed convex hull has the following {len(hyps)} exponent pairs as vertices."
    )
    for h in hyps:
        h.recursively_list_proofs()

    print()
    print(
        f"Proof that (3/40, 31/40) is an exponent pair, found by searching for the earliest proof"
    )
    eph = ep.find_best_proof(
        frac(3, 40), frac(31, 40), hypotheses, Proof_Optimization_Method.DATE
    )
    eph.recursively_list_proofs()
    print("Proof complexity:", eph.proof_complexity(), "Proof date:", eph.proof_date())

    print()
    print(
        f"Another proof that (3/40, 31/40) is an exponent pair, this time found by searching for the shortest proof"
    )
    eph = ep.find_best_proof(
        frac(3, 40), frac(31, 40), hypotheses, Proof_Optimization_Method.COMPLEXITY
    )
    eph.recursively_list_proofs()
    print("Proof complexity:", eph.proof_complexity(), "Proof date:", eph.proof_date())


def beta_bound_examples():

    hypotheses = Hypothesis_Set()  # Start with an empty hypothesis set

    print("1. Assume all unconditional bounds on \\beta(\\alpha) in the literature:")
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on beta")
    )
    display_beta_bounds(compute_best_beta_bounds(hypotheses))

    print()
    print("2. These beta bounds imply the following exponent pairs:")
    for h in beta_bounds_to_exponent_pairs(hypotheses):
        print(f"\t{h},\t which depends on")
        for h1 in sorted(h.dependencies, key=lambda x: x.data.bound.domain.x0):
            print(f"\t\t-{h1}  {h1.data}\t depends on {h1.dependencies}")

    print()
    print(
        "3. Now assume all exponent pairs in the literature, and the trivial pair (0, 1):"
    )
    hypotheses.add_hypothesis(trivial_exp_pair)
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair transform")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Exponent pair to beta bound transform")
    )
    print(
        "\tTotal exponent pairs:",
        len(hypotheses.list_hypotheses(hypothesis_type="Exponent pair")),
    )
    derived_pairs = compute_exp_pairs(
        hypotheses, search_depth=1
    )  # Compute the hull a bit
    hypotheses.add_hypotheses(derived_pairs)
    print(
        f"\tAlso assume derived exponent pairs with transform depth \\leq 1. Total pairs:",
        len(hypotheses.list_hypotheses(hypothesis_type="Exponent pair")),
    )
    print("\tComputing... (may take a few minutes)")
    start_time = time.time()
    hypotheses.add_hypotheses(exponent_pairs_to_beta_bounds(hypotheses))
    bounds = compute_best_beta_bounds(hypotheses)
    display_beta_bounds(bounds)
    print(f"\tComputed in {time.time() - start_time} sec")

    print()
    print("4. This beta bound implies the following exponent pairs:")
    hypotheses.add_hypotheses(bounds)
    new_exp_pairs = beta_bounds_to_exponent_pairs(hypotheses)
    new_exp_pairs.sort(key=lambda p: p.data.k)
    for h in new_exp_pairs:
        print(f"\t{h}")

    print(
        "5. Recursive list of hypotheses used to derive the (1101653/15854002, 12327829/15854002) exponent pair"
    )
    eph = next(h for h in new_exp_pairs if h.data.k == frac(1101653, 15854002))
    eph.recursively_list_proofs()

    # we could continue iterating, however let's check the implications for mu for now
    print()
    print("6. This implies the following bounds on mu in the range [1/3, 9/10)")
    hypotheses.add_hypotheses(new_exp_pairs)
    mbs = best_mu_bound_piecewise(frac(1, 3), frac(9, 10), hypotheses)
    for b in mbs:
        print(f"\\mu(x) \\leq {b}")


def large_values_examples():

    hypotheses = Hypothesis_Set()  # Start with an empty hypothesis set

    print("1. Assume all unconditional bounds on LV(\\sigma, \\tau) in the literature:")
    hypotheses.add_hypothesis(
        lv.large_value_estimate_L2
    )  # also add classical estimates
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Large value estimate")
    )

    raise_to_second_power = lv.raise_to_power_hypothesis(2)
    transformed_LV_estimate = raise_to_second_power.data.transform(
        lv.large_value_estimate_L2
    )
    print("Before transformation")
    lv.large_value_estimate_L2.recursively_list_proofs()
    print("After transformation")
    transformed_LV_estimate.recursively_list_proofs()

    print("2. Now assume transforms for k = 2, ..., k = 3")
    for k in range(3, 4):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))

    print("\t This implies the following bounds on LV(\\sigma, \\tau):")
    print("\t Computing... (may take a few minutes)")
    print("\t LV(x, y) \\leq {")
    hyps = lv.best_large_value_estimate(hypotheses)
    for h in hyps:
        print("\t\t", h.data.bound.pieces[0])
        h.recursively_list_proofs(3)
    print("\t}")
    print(f"Total {len(hyps)} pieces.")


def zeta_large_values_examples():
    hypotheses = Hypothesis_Set()  # Start with an empty hypothesis set
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on beta")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on mu")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Large value estimate")
    )

    lv_estimates = zlv.beta_to_zlv(hypotheses)
    lv_estimates.extend(zlv.mu_to_zlv(hypotheses))
    lv_estimates.extend(zlv.lv_to_zlv(hypotheses))

    for e in lv_estimates:
        print(e.desc_with_proof())


def plot(zdt, hypotheses, title=None):
    N = 500
    xs = []
    computed_zdt = []
    literature_zdt = []

    for i in range(N):
        sigma = 1 / 2 + 1 / 2 * i / N
        A = next((p.data.at(sigma) for p in zdt if p.data.interval.contains(sigma)), 0)
        xs.append(sigma)
        computed_zdt.append(A / (1 - sigma))
        literature_zdt.append(min(h.data.at(sigma) for h in hypotheses if h.data.interval.contains(sigma)))

    plt.figure(dpi=1200)
    plt.xlabel("σ")
    plt.ylabel("A(σ)")
    plt.plot(xs, computed_zdt, linewidth=0.5, label="Computed zero-density estimate")
    plt.plot(
        xs, literature_zdt, linewidth=0.5, label="Literature zero-density estimate"
    )
    if title is not None:
        plt.title(title)
    plt.legend(loc="lower left")
    plt.show()


def zero_density_estimates_examples():

    # Reset
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypothesis(lv.large_value_estimate_L2)
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Large value estimate")
    )
    for k in range(2, 15):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))

    hypotheses.add_hypotheses(
        zlv.get_trivial_zlv()
    )  # Include the tribial zeta large value estimate
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Zeta large value estimate")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Zero density estimate")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on beta")
    )
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on mu")
    )
    zdt = zd.lv_zlv_to_zd(hypotheses, Interval(frac(1, 2), 1), tau0=frac(5))
    print("Best-known vs computed zero-density estimate")
    print("A(x)(1-x) \\leq")
    for h in zdt:
        h.recursively_list_proofs()
    plot(
        zdt,
        literature.list_hypotheses(hypothesis_type="Zero density estimate"),
        "Best zero density estimate",
    )

    # hypotheses.add_hypothesis(trivial_exp_pair)
    # hypotheses.add_hypotheses(literature.list_hypotheses(hypothesis_type='Exponent pair'))
    # hypotheses.add_hypotheses(literature.list_hypotheses(hypothesis_type='Exponent pair transform'))

    # zdt = zd.ep_to_zd(hypotheses)
    # zdt.recursively_list_proofs()


def zero_density_estimates_examples2():
    hypotheses = Hypothesis_Set()

    # Large value estimates
    hypotheses.add_hypothesis(lv.large_value_estimate_L2)
    hypotheses.add_hypotheses(
        literature.find_hypothesis(
            hypothesis_type="Large value estimate", keywords="Jutila, k = 3"
        )
    )
    for k in range(2, 10):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))

    # Zeta large value estimates
    hypotheses.add_hypothesis(
        literature.find_hypothesis(
            hypothesis_type="Zeta large value estimate", keywords="Heath-Brown"
        )
    )

    zd.approximate_sup_LV_on_tau(hypotheses, lambda s: frac(7, 3) * s - frac(1, 3))
    zdt = zd.compute_best_zero_density_estimate(hypotheses, Interval(frac(1, 2), 1))
    print("Verifying Heath-Brown's zero-density estimate")
    print("A(x)(1-x) \\leq")
    for h in zdt:
        h.recursively_list_proofs()
    plot(zdt, "Heath-Brown zero density estimate")

    hypotheses = Hypothesis_Set()  # Start with an empty hypothesis set
    hypotheses.add_hypothesis(lv.large_value_estimate_L2)
    for k in range(2, 5):
        hypotheses.add_hypothesis(lv.raise_to_power_hypothesis(k))

def more_zero_density_examples():
    hyp = literature.find_hypothesis(hypothesis_type="Large value estimate", keywords="Bourgain")
    print(hyp.data)
    #zd.approx_optimise_bourgain_zero_density_estimate()
    #zd.optimise_bourgain_zero_density_estimate()

def all_examples():
    # mu_bound_examples()
    # exp_pair_examples()
    # beta_bound_examples()
    # large_values_examples()
    # zeta_large_values_examples()
    # zero_density_estimates_examples()
    # zero_density_estimates_examples2()
    more_zero_density_examples()

all_examples()
