from fractions import Fraction as frac


class Constants:

    # The maximum number of beta bounds in a series we will consider
    BETA_TRUNCATION = 20

    # The maximum number of exponent pairs in a series we will consider
    EXP_PAIR_TRUNCATION = 20

    # The maximum number of large value theorems in a series we will consider
    LARGE_VALUES_TRUNCATION = 20

    # The maximum value of the alpha parameter (reciprocal of tau)
    ALPHA_UPPER_LIMIT = frac(100)

    # Impose an upper limit on the \tau variable (appears in large value estimates)
    # to ensure that all polytopes are finite subsets of R^n (used to compute centroids).
    TAU_UPPER_LIMIT = frac(1000000)

    # The default upper bound on LV(sigma, tau), representing no upper bound
    LV_DEFAULT_UPPER_BOUND = frac(1000000)

    # In the calculation of zero-density theorems, only those estimates holding on an 
    # an interval [A, B] with A <= ZERO_DENSITY_SIGMA_LIMIT are considered. 
    ZERO_DENSITY_SIGMA_LIMIT = frac("0.999")


# Often many proofs exist for a single problem - this class allows us to select
# the objective function used to choose the 'best' proof.
class Proof_Optimization_Method:

    # No preference given to any proof
    NONE = 0

    # Find the proof that minimise the date of the most recent result in
    # the (recursive) dependency set of the proof
    DATE = 1

    # Find the proof with the lowest number of dependencies
    COMPLEXITY = 2
