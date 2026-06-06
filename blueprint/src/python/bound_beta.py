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

def apply_reflection_beta(bounds: list[Hypothesis]) -> list[Hypothesis]:
    """
    Implements Lemma 4.10 (van der Corput B-process for beta):
    
        beta(1 - alpha) = 1/2 - alpha + beta(alpha),  for 0 < alpha < 1
    
    Mathematical transformation:
        If beta(alpha) <= m*alpha + c  for alpha in [x0, x1]
        Then beta(alpha) <= (1-m)*alpha + (c + m - 1/2)  for alpha in [1-x1, 1-x0]
    
    Parameters:
        bounds: list of Hypothesis objects of type "Upper bound on beta"
    
    Returns:
        list of new Hypothesis objects representing the reflected bounds
    """
    if not bounds:
        return []

    reflected = []

    for bd in bounds:
        # Step 1: Extract the affine bound data
        p = bd.data.bound.deep_copy()
        m  = p.m
        c  = p.c
        x0 = p.domain.x0
        x1 = p.domain.x1

        # Step 2: Apply the reflection transformation (Lemma 4.10)
        # beta(alpha') <= (1-m)*alpha' + (c + m - 1/2)
        # for alpha' in [1 - x1, 1 - x0]
        new_m  = 1 - m
        new_c  = c + m - frac(1, 2)
        new_x0 = 1 - x1
        new_x1 = 1 - x0

        # Step 3: Build the new Affine bound with reflected interval
        new_affine = Affine(new_m, new_c, Interval(new_x0, new_x1, True, True))

        # Step 4: Wrap in a derived Hypothesis, preserving the dependency chain
        new_hyp = derived_bound_beta(
            new_affine,
            f'Follows from Lemma 4.10 (reflection) applied to "{bd.name}"',
            {bd},
        )
        reflected.append(new_hyp)

    # Step 5: Reverse the list because reflecting [x0, x1] -> [1-x1, 1-x0]
    # reverses the order of intervals on the number line
    reflected.reverse()

    return reflected

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


def apply_van_der_corput_process_for_beta(bounds: list[Hypothesis]) -> list[Hypothesis]:
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
def apply_van_der_corput_f1_eq_f2(bounds: list[Hypothesis]) -> list[Affine]:
    """
    Implements the optimized van der Corput process from Lemma 4.6 
    by picking h(alpha) such that f1(alpha) = f2(alpha).
    """
    if not bounds:
        return []

    # 1. Extend beta bounds to [0, 1] using reflection lemma
    boundData1 = [
        [h.data.bound.domain.x0, h.data.bound.domain.x1, h.data.bound.m, h.data.bound.c, h]
        for h in bounds
    ]
    boundData2 = [
        [1 - h.data.bound.domain.x1, 1 - h.data.bound.domain.x0,
         1 - h.data.bound.m, h.data.bound.c + h.data.bound.m - frac(1, 2), h]
        for h in bounds
    ]
    boundData2 = boundData2[::-1]
    boundData_all = boundData1 + boundData2

    newBounds = []

    # 2. Match each current bound with all auxiliary bounds
    for bd in bounds:
        p  = bd.data.bound.deep_copy()
        v0 = p.domain.x0
        v1 = p.domain.x1
        m1 = p.m
        c1 = p.c

        for dat in boundData_all:
            X0, X1, m2, c2 = dat[0], dat[1], dat[2], dat[3]

            if m2 + c2 + 1 == 0:
                continue

            # Helper to solve the quadratic equation exactly in Q
            def solve_t(alpha):
                L1_val = m1 * alpha + c1
                if L1_val <= 0:
                    return None
                A_val = (1 - alpha) * m2
                B_val = (1 - alpha) * c2 + L1_val * (1 - 2 * alpha)
                C_val = alpha * L1_val

                if A_val == 0:
                    if B_val == 0:
                        return None
                    return frac(C_val, B_val)
                else:
                    disc = B_val * B_val + 4 * A_val * C_val
                    if disc < 0:
                        return None
                    
                    import math
                    num, den = disc.numerator, disc.denominator
                    sqrt_num, sqrt_den = math.isqrt(num), math.isqrt(den)
                    
                    # Ensure perfect square to maintain exact fraction arithmetic
                    if sqrt_num * sqrt_num != num or sqrt_den * sqrt_den != den:
                        return None  
                    
                    sqrt_disc = frac(sqrt_num, sqrt_den)
                    t = frac(-B_val + sqrt_disc, 2 * A_val)
                    if t < 0:
                        t = frac(-B_val - sqrt_disc, 2 * A_val)
                    return t if t >= 0 else None

            def B_val(alpha, t):
                return (1 - alpha) * (m2 * t + c2)

            def h_val(alpha, t):
                return frac(alpha * (1 + t), t) - 1 if alpha != 0 else None

            # Evaluate at endpoints to build the conservative affine interpolant
            results = []
            for alpha_test in [v0, v1]:
                t_sol = solve_t(alpha_test)
                if t_sol is None or t_sol < X0 or t_sol > X1:
                    results.append(None)
                    continue
                h_sol = h_val(alpha_test, t_sol)
                if h_sol is None or h_sol <= 0:
                    results.append(None)
                    continue
                results.append((alpha_test, B_val(alpha_test, t_sol)))

            if results[0] is None or results[1] is None:
                continue

            (a0, b0), (a1, b1) = results
            if a0 == a1:
                continue

            # Linear interpolation between endpoints
            new_m = frac(b1 - b0, a1 - a0)
            new_c = b0 - new_m * a0

            if new_m == m1 and new_c >= c1:
                continue  

            # Find sub-interval where this new bound actually improves the current one
            u0, u1 = v0, v1
            dm = new_m - m1
            dc = c1 - new_c
            if dm > 0:
                u1 = min(u1, frac(dc, dm))
            elif dm < 0:
                u0 = max(u0, frac(dc, dm))

            if u1 > u0:
                newBounds.append(Affine(new_m, new_c, Interval(u0, u1, True, True)))

    if not newBounds:
        return []

    # 3. Clean up and merge identical segments
    newBounds.sort(key=lambda f: (f.m, f.c, f.domain.x0))
    merged = []
    cur = newBounds[0].deep_copy()
    for nb in newBounds[1:]:
        if nb.m == cur.m and nb.c == cur.c and nb.domain.x0 <= cur.domain.x1:
            cur = Affine(cur.m, cur.c, Interval(cur.domain.x0, max(cur.domain.x1, nb.domain.x1), True, True))
        else:
            merged.append(cur)
            cur = nb.deep_copy()
    merged.append(cur)

    return merged

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
    def combine_beta_bounds(hypotheses: list) -> list:
    """Compute the pointwise minimum and return a non-redundant list of Hypothesis."""
    if not hypotheses:
        return []

    # 1. Collect all domain endpoints and intersection points
    breakpoints = set()
    bounds_data = []
    
    for hyp in hypotheses:
        aff = hyp.data.bound
        x0, x1, m, c = aff.domain.x0, aff.domain.x1, aff.m, aff.c
        bounds_data.append((m, c, x0, x1, hyp))
        breakpoints.add(x0)
        breakpoints.add(x1)

    n = len(bounds_data)
    for i in range(n):
        mi, ci, x0i, x1i, _ = bounds_data[i]
        for j in range(i + 1, n):
            mj, cj, x0j, x1j, _ = bounds_data[j]
            if mi != mj:
                alpha_cross = frac(cj - ci, mi - mj)
                if x0i <= alpha_cross <= x1i and x0j <= alpha_cross <= x1j:
                    breakpoints.add(alpha_cross)

    breakpoints = sorted(breakpoints)

    # 2. Find which hypothesis achieves the minimum on each sub-interval
    result_pieces = []
    for k in range(len(breakpoints) - 1):
        left, right = breakpoints[k], breakpoints[k + 1]
        if left == right:
            continue

        mid = frac(left + right, 2)
        best_val, best_hyp, best_m, best_c = None, None, None, None

        for (m, c, x0, x1, hyp) in bounds_data:
            if x0 <= mid <= x1:
                val = m * mid + c
                if best_val is None or val < best_val:
                    best_val, best_hyp, best_m, best_c = val, hyp, m, c

        if best_hyp is not None:
            result_pieces.append((best_m, best_c, left, right, best_hyp))

    if not result_pieces:
        return []

    # 3. Merge adjacent intervals sharing identical affine functions
    merged = []
    cur_m, cur_c, cur_left, cur_right, cur_hyp = result_pieces[0]

    for i in range(1, len(result_pieces)):
        m, c, left, right, hyp = result_pieces[i]
        if m == cur_m and c == cur_c and hyp is cur_hyp and left == cur_right:
            cur_right = right
        else:
            merged.append((cur_m, cur_c, cur_left, cur_right, cur_hyp))
            cur_m, cur_c, cur_left, cur_right, cur_hyp = m, c, left, right, hyp

    merged.append((cur_m, cur_c, cur_left, cur_right, cur_hyp))

    # 4. Build the final clean Hypothesis objects list
    result = []
    for (m, c, x0, x1, source_hyp) in merged:
        new_affine = Affine(m, c, Interval(x0, x1, True, True))
        if new_affine == source_hyp.data.bound:
            result.append(source_hyp)
        else:
            result.append(
                derived_bound_beta(
                    new_affine,
                    f"Pointwise minimum of beta bounds; tightest bound on [{x0}, {x1}] is from: {source_hyp.name}",
                    {source_hyp},
                )
            )

    return result
