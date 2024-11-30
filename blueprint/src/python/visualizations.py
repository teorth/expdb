import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from literature import *
from functions import RationalFunction as RF


def van_der_corput(k, alpha):
    return max(
        alpha + (1 - k * alpha) / (2**k - 2),
        (1 - 2 ** (2 - k)) * alpha - (1 - k * alpha) / (2**k - 2),
    )


def optimized_van_der_corput(alpha):
    bound = 1
    for k in range(2, 30):
        bound = min(bound, van_der_corput(k, alpha))
    return bound


def conjectured_bound(alpha):
    if alpha <= 1:
        return alpha / 2
    else:
        return alpha - 1


def van_der_corput_plot():
    alpha_range = np.linspace(0, 1, 500)
    plt.figure(figsize=(10, 6))

    plt.plot(
        alpha_range,
        [van_der_corput(2, alpha) for alpha in alpha_range],
        label=r"$k=2$ van der Corput",
    )
    plt.plot(
        alpha_range,
        [van_der_corput(3, alpha) for alpha in alpha_range],
        label=r"$k=3$ van der Corput",
    )
    plt.plot(
        alpha_range,
        [van_der_corput(4, alpha) for alpha in alpha_range],
        label=r"$k=4$ van der Corput",
    )
    plt.plot(
        alpha_range,
        [van_der_corput(5, alpha) for alpha in alpha_range],
        label=r"$k=5$ van der Corput",
    )
    plt.plot(
        alpha_range,
        [optimized_van_der_corput(alpha) for alpha in alpha_range],
        label=r"Optimized van der Corput",
    )
    plt.plot(
        alpha_range, [alpha / 2 for alpha in alpha_range], label=r"Conjectured value"
    )
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta(\alpha)$")
    plt.title(r"Classical bounds on $\beta(\alpha)$")
    plt.legend()
    plt.grid(True)
    plt.ylim(0, 1)
    plt.show()


def van_der_corput_plot2():
    alpha_range = np.linspace(0, 1, 500)
    alpha_full_range = np.linspace(0, 2, 1000)

    plt.figure(figsize=(10, 6))

    plt.plot(
        alpha_full_range,
        [optimized_van_der_corput(alpha) for alpha in alpha_full_range],
        label=r"Optimized van der Corput",
    )
    plt.plot(
        alpha_range, [alpha for alpha in alpha_range], label=r"Trivial upper bound"
    )
    plt.plot(
        alpha_full_range,
        [conjectured_bound(alpha) for alpha in alpha_full_range],
        label=r"Conjectured value",
    )
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta(\alpha)$")
    plt.title(r"Classical bounds on $\beta(\alpha)$")
    plt.legend()
    plt.grid(True)
    plt.ylim(0, 1)
    plt.show()


def beta_bound_plot():
    hypotheses = Hypothesis_Set()  # Start with an empty hypothesis set
    hypotheses.add_hypothesis(trivial_beta_bound_1)
    hypotheses.add_hypotheses(
        literature.list_hypotheses(hypothesis_type="Upper bound on beta")
    )
    best_beta_bounds = compute_best_beta_bounds(hypotheses)

    alpha_range = np.linspace(0, 1 / 2, 500)
    plt.figure(figsize=(10, 6))
    plt.plot(
        alpha_range,
        [
            min(
                b.data.bound.at(alpha)
                for b in best_beta_bounds
                if b.data.bound.domain.contains(alpha)
            )
            for alpha in alpha_range
        ],
        label=r"Computed best upper bound on $\beta$",
    )
    plt.plot(
        alpha_range,
        [optimized_van_der_corput(alpha) for alpha in alpha_range],
        label=r"van der Corput bound on $\beta$",
    )
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta(\alpha)$")
    plt.title(r"Best bounds on $\beta(\alpha)$")
    plt.legend()
    plt.grid(True)
    plt.ylim(0, 1 / 2)
    plt.show()

    for b in best_beta_bounds:
        print(b.data)
    print(
        min(
            b.data.bound.at(0)
            for b in best_beta_bounds
            if b.data.bound.domain.contains(0)
        )
    )

def exp_pair_plot():
    exp_pairs = [
        (0, 1),
        (1/980688, 3105401/3105512),
        (1/490512, 9319081/9319728),
        (1/326808, 1552181/1552338),
        (1/163416, 1034765/1034968),
        (1/81720, 155209/155268),
        (1/40872, 775997/776568),
        (1/27230, 258409/258685),
        (1/13614, 258133/258666),
        (1/6806, 64400/64657),
        (1/3402, 7127/7182),
        (1/2890, 1517/1530),
        (1/2432, 5119/5168),
        (1/2025, 2137/2160),
        (1/1666, 3527/3570),
        (1/1352, 359/364),
        (1/1080, 2303/2340),
        (1/847, 907/924),
        (1/650, 1399/1430),
        (1/486, 263/270),
        (1/352, 767/792),
        (1/245, 269/280),
        (1/162, 359/378),
        (1/100, 14/15),
        (89/3478, 15327/17390),
        (10769/351096, 609317/702192),
        (9/217, 1461/1736),
        (2371/43205, 280013/345640),
        (652397/9713986, 7599781/9713986),
        (89/1282, 997/1282),
        (2779/38033, 58699/76066),
        (18/199, 593/796),
        (4742/38463, 35731/51284),
        (13/84, 55/84)
    ]
    n = len(exp_pairs)
    for i in range(n):
        p = exp_pairs[n - i - 1]
        exp_pairs.append((p[1] - 1/2, p[0] + 1/2))
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.add_patch(Polygon(exp_pairs, color=[0.7, 0.7, 0.7]))
    ax.set_xlim([0, 1/2])
    ax.set_ylim([1/2, 1])
    ax.scatter([p[0] for p in exp_pairs], [p[1] for p in exp_pairs], color="k")
    
    plt.ylabel(r"$\ell$")
    plt.xlabel(r"$k$")
    plt.title(r"Convex hull of known exponent pairs $(k, \ell)$")
    plt.show()

# Compare the best literature zero-density estimates against the 
# best-known zero-density estimate
def zero_density_plot():
    # plot zero-density estimate
    N = 500
    sigmas = []
    literature_zd = []
    best_zd = []
    lindelof_zd = []
    
    hs = Hypothesis_Set()
    hs.add_hypotheses(literature)
    
    # Add the new zero-density estimates (not part of the literature yet!)
    zd.add_zero_density(hs, "2/(9*x - 6)", Interval("[17/22, 38/49]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "9/(8*(2*x - 1))", Interval("[38/49, 4/5]"), Reference.make("Tao--Trudgian--Yang", 2024))
    zd.add_zero_density(hs, "3/(10 * x - 7)", Interval("[701/1000, 1]"), Reference.make("Tao--Trudgian--Yang", 2024))
    hs.add_hypotheses(zd.bourgain_ep_to_zd())
    
    def best_zd_at(hypotheses, sigma):
        min_ = 1000000
        for h in hypotheses:
            if h.hypothesis_type == "Zero density estimate":
                if h.data.interval.contains(sigma):
                    q = h.data.at(sigma)
                    if q.is_real and math.isfinite(q) and min_ > q:
                        min_ = q
        return min_
    
    for i in range(N):
        sigma = 1 / 2 + 1 / 2 * i / N
        sigmas.append(sigma)
        literature_zd.append(best_zd_at(literature, sigma))
        best_zd.append(best_zd_at(hs, sigma))
        lindelof_zd.append(2 if sigma <= 3/4 else 0)
    
    for i in range(len(sigmas)):
        print(sigmas[i], literature_zd[i], best_zd[i])
        
    plt.figure(figsize=(10, 6))
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$A(\sigma)$")
    plt.plot(sigmas, literature_zd, label="Literature zero-density estimate")
    plt.plot(sigmas, lindelof_zd, label="Zero-density estimate under the Lindelof hypothesis")
    plt.plot(sigmas, best_zd, label="Best zero-density estimate")
    plt.title("Zero density estimates")
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.show()

# Plot the best zero-density energy estimates known trivially, in the literature, and under the 
# Lindelof hypothesis
def zero_density_energy_plot():
    hypotheses = Hypothesis_Set()
    hypotheses.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zero density estimate"))
    hypotheses.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zero density energy estimate"))

    # add trivial bounds - this uses literature zero-density estimates
    ze.add_trivial_zero_density_energy_estimates(hypotheses)
    
    # List of new derived estimates so far. TODO: replace with actual derivations
    energy_estimates = [
            (RF.parse("1000000"), Interval(frac(1,2), 1)), # default
            (RF.parse("5 * (18 - 19 * x) / (2 * (5 * x + 3) * (1 - x))"), Interval(frac(7,10), 0.7255)),
            (RF.parse("2 * (45 - 44 * x) / ((2 * x + 15) * (1 - x))"), Interval(0.7255, frac(3,4))),
            (RF.parse("(197 - 220 * x) / (8 * (5 * x - 1) * (1 - x))"), Interval(frac(3,4), frac(25,33))),
            (RF.parse("(110 * x - 53) / (16 * (85 * x - 63) * (1 - x))"), Interval(frac(25,33), frac(3039,4010))),
            (RF.parse("(570 * x - 413) / (10 * (85 * x - 63) * (1 - x))"), Interval(frac(3039,4010), frac(254,335))),
            (RF.parse("(66 - 83 * x) / ((69 * x - 50) * (1 - x))"), Interval(frac(254,335), 0.758233)),
            (RF.parse("5 * (18 - 19 * x) / (2 * (5 * x + 3) * (1 - x))"), Interval(0.758233, frac(42,55))),
            (RF.parse("(18 - 19 * x) / (6 * (15 * x - 11) * (1 - x))"), Interval(frac(42,55), frac(97,127))),
            (RF.parse("3 * (18 - 19 * x) / (4 * (4 * x - 1) * (1 - x))"), Interval(frac(97,127), frac(79,103))),
            (RF.parse("(18 - 19 * x) / (2 * (37 * x - 27) * (1 - x))"), Interval(frac(79,103), frac(33,43))),
            (RF.parse("5 * (18 - 19 * x) / (2 * (13 * x - 3) * (1 - x))"), Interval(frac(33,43), frac(84,109))),
            (RF.parse("(18 - 19 * x) / (9 * (3 * x - 2) * (1 - x))"), Interval(frac(33,43), 0.7721)),
            (RF.parse("4 * (10 - 9 * x) / (5 * (4 * x - 1) * (1 - x))"), Interval(0.7721, frac(5, 6))),
        ]
    # as an example, plot the trivial bound and the literature bound 
    sigmas = np.linspace(1/2, 0.999, 1000)
    trivial_bound = []
    lit_bound = []
    best_bound = []
    lindelof_bound = []
    for sigma in sigmas:
        trivial = min(
                    h.data.at(sigma) 
                    for h in hypotheses 
                    if h.name == "Trivial zero density energy estimate" and h.data.interval.contains(sigma)
                )
        trivial_bound.append(trivial)
        
        lit = min(
                h.data.at(sigma)
                for h in literature
                if h.hypothesis_type == "Zero density energy estimate" and h.data.interval.contains(sigma)
            )
        lit_bound.append(min(lit, trivial))
        
        best = min(
                bound.at(sigma)
                for (bound, interval) in energy_estimates
                if interval.contains(sigma)
            )
        best_bound.append(min(best, lit, trivial))
        
        if sigma > 3/4:
            lindelof_bound.append(0)
        else:
            lindelof_bound.append(8 - 4 * sigma)
    
    plt.figure(figsize=(10, 6))
    plt.xlabel(r"$\sigma$")
    plt.ylabel(r"$A^*(\sigma)$")
    plt.plot(sigmas, lit_bound, linewidth=1, label="Literature additive energy estimate")
    plt.plot(sigmas, trivial_bound, linewidth=1, label="\"Trivial\" additive energy estimate")
    plt.plot(sigmas, best_bound, linewidth=1, label="Best additive energy estimate")
    plt.plot(sigmas, lindelof_bound, linewidth=1, label="Additive energy estimate on the Lindelof hypothesis")
    plt.title("Additive energy estimates")
    plt.legend()
    plt.grid(True)
    plt.show()


# van_der_corput_plot2()
# beta_bound_plot()
exp_pair_plot()
# zero_density_plot()
#zero_density_energy_plot()

