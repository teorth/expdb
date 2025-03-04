# code for representing a bound on the exponential sum
# \sum_{n \asymp N}e(TF(n/N)) \ll N^{\beta(\alpha) + o(1)}
# where $\alpha = \log N/\log T.

from constants import *
from hypotheses import *

# from mpmath import fraction as frac
from fractions import Fraction as frac
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
        "Classical bound on \\beta",
        "Upper bound on beta",
        Bound_beta(bound),
        "Classical",
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


def apply_van_der_corput_process_for_beta(bounds):
    # Lemma 4.6 gives a bound of the form \beta(\alpha) <= max( f1(alpha), f2(alpha), f3(alpha) )

    # f3 contains a sup over a term h'. If bounds is piecewise linear with each intercept C >= -1 and continuous then we can take h' = h
    bd0 = bounds[0]
    p0 = bd0.data.bound.deep_copy()
    C = p0.c
    M = p0.m
    if (C < -1) or ( C + M - frac(1,2) < -1):
        return []
    for i in range(len(bounds)-1):
        bd1 = bounds[i]
        bd2 = bounds[i+1]
        p1 = bd1.data.bound.deep_copy()
        p2 = bd2.data.bound.deep_copy()

        v1 = p1.domain.x1
        v2 = p2.domain.x0
        C = p2.c
        M = p2.m
        if p1.at(v1, True) != p2.at(v2, True):
            return []
        if (C < -1) or ( C + M - frac(1,2) < -1):
            return []

    # extend beta bounds to [0.5, 1]
    boundData1 = [[ h.data.bound.domain.x0, h.data.bound.domain.x1, h.data.bound.m, h.data.bound.c ]  for h in bounds ]
    boundData2 = [[ 1 - h.data.bound.domain.x1, 1 - h.data.bound.domain.x0, 1 - h.data.bound.m, h.data.bound.c  + h.data.bound.m - frac(1,2) ]  for h in bounds ]
    boundData2 = boundData2[::-1]
    boundData = boundData1 + boundData2

    newBounds = []
    for bd in bounds:
        p = bd.data.bound.deep_copy()
        v0 = p.domain.x0
        v1 = p.domain.x1
        domain1 = p.domain.deep_copy()
        m1 = p.m
        c1 = p.c

        # We set h = ((1+c2+m2)alpha)/(2+2*c2) + c2/(2+2*c2)
        # In this case f1(alpha) = f3(alpha)
        for i in range(1,len(boundData)):
            dat = boundData[i]
            m2 = dat[2]
            c2 = dat[3]
            X0 = dat[0]
            X1 = dat[1]

            # Skip over what would generate current best bounds
            if frac(1+c2+m2,2+2*c2) == m1 and frac( c2,2+2*c2) == c1:
                continue

            # To apply this bound we require
            # X0 <= alpha/(1 + h - alpha) <= X1
            domain2 = domain1.intersect( Interval( frac( X0, 1 + c2 + m2*X0), frac( X1, 1 + c2 + m2*X1),   True, True ) )
            if domain2.length() == 0:
                continue

            # check that f2(alpha) >= f1(alpha) and is an improvement over the current beta bound
            u0 = domain2.x0
            u1 = domain2.x1
            if  1 + c2 + m2 - 2*m1*(1 + c2) == 0:
                if c2 > 2*c1*(1+c2):
                    continue
            elif  1 + c2 + m2 - 2*m1*(1 + c2) > 0:
                u1 = min( [u1, frac(2*c1*(1+c2) - c2, 1 + c2 + m2 - 2*m1*(1 + c2)) ] )
            else:
                u0 = max( [u0, frac(2*c1*(1+c2) - c2, 1 + c2 + m2 - 2*m1*(1 + c2)) ] )
            if  3*m2 - c2 - 1 < 0:
                u1 = min( [u1, frac( -3*c2, 3*m2 - c2 - 1) ] )
            if u1 > u0:
                newBounds.append( Affine( frac(1+c2+m2,2+2*c2), frac( c2,2+2*c2), Interval( u0, u1, True, True) ) )
        #TODO check if picking h so that f1(alpha) = f2(alpha) achieves any new bounds

    if len(newBounds) > 0:      
        # Merge the same bound that appears over multiple intervals
        newBounds1 = [nb.deep_copy() for nb in newBounds]
        newBounds2 = []
        while len(newBounds1) > 0:
            f1 = newBounds1[0].deep_copy()
            fm = f1.m
            fc = f1.c 
            sameFs0 = [ nb for nb in newBounds1 if (nb.m == fm) and (nb.c == fc) ]
            sameFs0.sort(  key=lambda f: f.domain.x0 )
            sameFs = []
            for SF in sameFs0:
                if not(SF in sameFs):
                    sameFs.append(SF.deep_copy())
            curF = sameFs[0].deep_copy()
            for i in range(1, len(sameFs)):
                newF = sameFs[i].deep_copy()
                u0 = curF.domain.x0
                u1 = curF.domain.x1
                v0 = newF.domain.x0
                v1 = newF.domain.x1
                if v0 <= u1:
                    curF = Affine( curF.m, curF.c, Interval(curF.domain.x0, max( [curF.domain.x1, newF.domain.x1] ), True, True) )
                else:
                    newBounds2.append(curF.deep_copy())
                    curF = sameFs[i].deep_copy()
            newBounds2.append(curF.deep_copy())
            newBounds1 = [ nb for nb in newBounds1 if (nb.m != fm) or (nb.c != fc) ]

        return [nb for nb in newBounds2]
    return []


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


def display_two_sets_of_beta_bounds(hypotheses, newhypotheses):

    for h in hypotheses:
        p = h.data.bound
        plt.plot(
            (p.domain.x0, p.domain.x1),
            (p.at(p.domain.x0, True), p.at(p.domain.x1, True)),
            color="black",
        )

    for h in newhypotheses:
        p = h.data.bound
        pDom = p.domain
        p0 = pDom.x0
        p1 = pDom.x1
        u0 = p1
        u1 = p0
        for m in hypotheses:
            q = m.data.bound
            if (p.m == q.m) and (p.c == q.c):
                qDom = q.domain
                q0 = qDom.x0
                q1 = qDom.x1
                if q0 <= u0:
                    u0 = q0
                if q1 >= u1:
                    u1 = q1
        if u0 > p0:
            plt.plot(
                (p0, u0),
                (p.at(p0, True), p.at(u0, True)),
                color="red",
            )
        if u1 < p1:
            plt.plot(
                (u1, p1),
                (p.at(u1, True), p.at(p1, True)),
                color="red",
            )

    plt.xlim((-0.05, 0.55))
    plt.ylim((-0.05, 0.5))
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\beta(\alpha)$")
    plt.title(r"Best bound on $\beta(\alpha)$")
    plt.show()