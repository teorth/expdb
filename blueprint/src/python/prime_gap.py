# This file contains calculations related to prime gaps (Ch. 13)

from fractions import Fraction as frac
from functions import Interval, RationalFunction as RF
from hypotheses import *
import zero_density_estimate as zd
import zero_density_energy_estimate as ze
import sympy

#: Anything below this in absolute value counts as zero imaginary part.  A real
#: root of a cubic is often returned in the casus irreducibilis form, i.e. as a
#: sum of complex radicals, so an exactly-zero imaginary part cannot be relied on.
_IMAG_TOL = 1e-9


def _real_value(expr, x, p):
    """The value of `expr` at `x = p` as a float, discarding a numerical zero
    imaginary part.  float() alone raises on a real root written with complex
    radicals."""
    v = complex(sympy.N(expr.subs(x, p)))
    if abs(v.imag) > _IMAG_TOL:
        raise ValueError(f"expected a real value at x = {p}, got {v}")
    return v.real


def _critical_points(expr, x, interval):
    """The stationary points of `expr` lying in `interval`, together with its
    endpoints.  Roots whose realness sympy cannot decide symbolically are kept
    when they evaluate to a real number."""
    pts = set()
    for p in sympy.solve(sympy.diff(expr, x), x):
        v = complex(sympy.N(p))
        if abs(v.imag) > _IMAG_TOL:
            continue
        if interval.contains(v.real):
            pts.add(p)
    pts.update([interval.x0, interval.x1])
    return pts


# Compute the best estimate of \theta_{gap, 2} using Proposition 15.9
def compute_gap2(hypotheses, debug=False):
    if not isinstance(hypotheses, Hypothesis_Set):
        raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")

    best_zd = zd.best_zero_density_estimate(hypotheses, verbose=False)
    for zdh in best_zd:
        print(zdh.data)

    best_ze = ze.compute_best_energy_bound(hypotheses)
    for zeh in best_ze:
        print(zeh.data)

    # Set of endpoints
    crits = {frac(1,2), frac(1)}
    crits.update(zh.data.interval.x0 for zh in best_zd)
    crits.update(zh.data.interval.x1 for zh in best_zd)
    crits.update(zh.data.interval.x0 for zh in best_ze)
    crits.update(zh.data.interval.x1 for zh in best_ze)
    crits = list(crits)
    crits.sort()

    x = RF.x

    for i in range(1, len(crits)):
        interval = Interval(crits[i - 1], crits[i])
        mid = interval.midpoint()

        zdh = next(h for h in best_zd if h.data.interval.contains(mid))
        zeh = next(h for h in best_ze if h.data.interval.contains(mid))

        # For simplicity - just work with sympy objects for now (in the future
        # these will be replaced with RationalFunction instances)
        A = zdh.data.bound.num / zdh.data.bound.den
        B = zeh.data.bound.num / zeh.data.bound.den

        alpha = 4 * x - 2 + 2 * (B * (1 - x) - 1) / (B - A)
        beta = 4 * x - 2 + (B * (1 - x) - 1) / A

        # Find all critical points and evaluate alpha at each
        statpts = _critical_points(alpha, x, interval)
        sup_alpha = max(_real_value(alpha, x, p) for p in statpts)
        if debug:
            print("alpha -------------------------------------------------")
            for p in statpts:
                print(p, _real_value(x, x, p), alpha.subs(x, p), _real_value(alpha, x, p))

        # Do the same for beta
        statpts = _critical_points(beta, x, interval)
        sup_beta = max(_real_value(beta, x, p) for p in statpts)
        if debug:
            print("beta --------------------------------------------------")
            for p in statpts:
                print(p, _real_value(x, x, p), beta.subs(x, p), _real_value(beta, x, p))

        print(interval, max(sup_alpha, sup_beta), alpha.simplify(), beta.simplify())


def prime_excep(hypotheses, DISCRETIZATION=100):
  """
  Bound the exponent mu_PNT(theta) for the number of exceptions to the prime number theorem
  in short intervals [x,x+x^theta], using the inequality provided by Gafni-Tao, and plot the
  bound over theta.  This computation is currently numerical with a discretization error; a
  future project would be to perform this computation symbolically.
  """
  if not isinstance(hypotheses, Hypothesis_Set):
      raise ValueError("Parameter hypotheses must be of type Hypothesis_Set")

  best_zd = zd.best_zero_density_estimate(hypotheses, verbose=False)
  for zdh in best_zd:
    print(zdh.data)

  best_ze = ze.compute_best_energy_bound(hypotheses)
  for zeh in best_ze:
    print(zeh.data)

  thetas = []
  bounds = []
  for j in range(DISCRETIZATION):
    theta = j / DISCRETIZATION
    bound = 0

    for i in range(0, DISCRETIZATION):
      sigma = i / DISCRETIZATION
      if sigma <= 1/2:
        A = 1 / (1 - sigma)
        Ae = 3 / (1 - sigma)
      else:
        A = next((b.data.at(sigma) for b in best_zd if b.data.interval.contains(sigma)), 0)
        Ae = next((b.data.at(sigma) for b in best_ze if b.data.interval.contains(sigma)), 0)
      if A < 1/(1-theta) - 1/DISCRETIZATION:
        continue
      mu2 = (1-theta) * (1-sigma) * A + 2*sigma - 1
      mu4 = (1-theta) * (1-sigma) * Ae + 4*sigma - 3
      bound = max(bound, min(mu2, mu4))

    bound = min(bound,1)
    thetas.append(theta)
    bounds.append(bound)

  import matplotlib.pyplot as plt

  plt.plot(thetas, bounds)
  plt.xlabel(r"$\theta$")
  plt.ylabel(r"$\mu_{PNT}(\theta)$")
  plt.title(r"Bound on $\mu_{PNT}(\theta)$")
  plt.grid()
  plt.show()
