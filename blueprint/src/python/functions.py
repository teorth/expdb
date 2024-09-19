# Code for some commonly encountered functions
import cdd
import copy
from fractions import Fraction as frac
import itertools
import numbers
import matplotlib.pyplot as plt
from polytope import *
import sympy


# Represents an interval on the real line. Occasionally we need to distinguish
# between closed and open intervals, hence the need for a custom class instead
# of implementing the in-built python interval class.
class Interval:

    # Creates an object representing the real interval [x0, x1). There are two
    # initialisation schemes:
    # Scheme 1: (x0, x1, include_lower, include_upper), with parameters
    #   - x0: (numerical type) the lower limit of the interval
    #   - x1: (numerical type) the upper limit of the interval
    #   - include_lower: (boolean type) If include_lower is True (defaults to
    #                   True), then the lower limit x0 is included in the interval.
    #   - include_upper: (boolean type) If include_upper is True (defaults to
    #                   False), then the upper limit x1 is included in the interval.
    #
    # Scheme 2: (expr), with parameters
    #   - expr (string type) a string expression like (a, b) or [a, b]
    def __init__(self, x0, x1=None, include_lower=True, include_upper=False):

        if isinstance(x0, str):
            i = Interval.parse(x0)
            self.x0 = i.x0
            self.x1 = i.x1
            self.include_lower = i.include_lower
            self.include_upper = i.include_upper
        else:
            # if not isinstance(x0, numbers.Number):
            #     raise ValueError('x0 must be of type number', x0)
            # if not isinstance(x1, numbers.Number):
            #     raise ValueError('x1 must be of type number', x1)

            if x0 - x1 > 0:
                raise ValueError("x0 must be <= x1")
            self.x0 = x0
            self.x1 = x1
            self.include_lower = include_lower
            self.include_upper = include_upper

    def __repr__(self):
        lower_bracket = "[" if self.include_lower else "("
        upper_bracket = "]" if self.include_upper else ")"
        return f"{lower_bracket}{self.x0},{self.x1}{upper_bracket}"

    def __eq__(self, other):
        # All empty intervals are considered equal
        if self.is_empty() and other.is_empty():
            return True
        return (
            self.x0 == other.x0
            and self.x1 == other.x1
            and self.include_lower == other.include_lower
            and self.include_upper == other.include_upper
        )

    # -------------------------------------------------------------------------
    # Public static functions
    def parse(s):
        inner = s[1:-1].split(",")
        return Interval(
            frac(inner[0].strip()), frac(inner[1].strip()), s[0] == "[", s[-1] == "]"
        )

    # -------------------------------------------------------------------------
    def deep_copy(self):
        return Interval(self.x0, self.x1, self.include_lower, self.include_upper)

    def length(self):
        return self.x1 - self.x0

    # Returns true if the interval contains any elements
    def is_empty(self):
        if self.length() == 0:
            return not (self.include_lower and self.include_upper)
        return False

    # Returns the midpoint of the interval
    def midpoint(self):
        return (self.x0 + self.x1) / 2

    # Returns true if the interval contains the point x
    def contains(self, x):
        # previous implementation required 6 comparisons, new implementation only 
        # requires 4 comparisons (but more cases)
        if self.include_lower:
            if self.include_upper:
                return x - self.x0 >= 0 and self.x1 - x >= 0
            else:
                return x - self.x0 >= 0 and self.x1 - x > 0

        if self.include_upper:
            return x - self.x0 > 0 and self.x1 - x >= 0
        else:
            return x - self.x0 > 0 and self.x1 - x > 0


    # Returns an Interval object representing the intersection of this interval
    # with another.
    def intersect(self, other):
        lower = max(self.x0, other.x0)
        upper = min(self.x1, other.x1)

        # Degenerate interval
        if lower > upper:
            return Interval(0, 0, False, False)  # Default

        if self.x0 > other.x0:
            include_lower = self.include_lower
        elif self.x0 < other.x0:
            include_lower = other.include_lower
        else:
            include_lower = self.include_lower and other.include_lower

        if self.x1 < other.x1:
            include_upper = self.include_upper
        elif self.x1 > other.x1:
            include_upper = other.include_upper
        else:
            include_upper = self.include_upper and other.include_upper

        return Interval(lower, upper, include_lower, include_upper)


# Represents an affine function
class Affine:

    # Creates an object representing the function f(x) = mx + c where m, c are
    # constants, in the interval x \in domain.
    # Parameters:
    #   - m: (numerical type) the slope of the affine function
    #   - c: (numerical type) the intercept of the affine function
    def __init__(self, m, c, domain, label=None):
        self.m = m
        self.c = c
        self.domain = domain
        self.label = label

    def __repr__(self):
        s1 = "x" if self.m == 1 else f"{self.m}x"
        s2 = ""
        if self.c > 0:
            s2 = f" + {self.c}"
        elif self.c < 0:
            s2 = f" - {abs(self.c)}"
        s3 = f"   ({self.label})" if self.label is not None else ""
        return f"{s1}{s2}  for  x \\in {self.domain}{s3}"

    def __eq__(self, other):
        if isinstance(other, Affine):
            # Special handling of singleton sets
            if (
                self.domain.length() == 0
                and other.domain.length() == 0
                and self.domain.x0 == other.domain.x0
            ):
                return self.at(self.domain.x0) == other.at(self.domain.x0)

            return (self.m, self.c, self.domain) == (other.m, other.c, other.domain)
        return NotImplemented

    def deep_copy(self):
        return Affine(self.m, self.c, self.domain.deep_copy())

    # Evaluates the function at x. If extend_domain is True, then the domain of
    # the function will be ignored.
    def at(self, x, extend_domain=False):
        if extend_domain or self.domain.contains(x):
            return self.m * x + self.c
        return None

    # Returns true if this function is the same as another (possibly on different)
    # domains
    def function_equals(self, other):
        return self.m == other.m and self.c == other.c

    # Returns the point of intersection with another Affine object if they
    # intersect, and None if not
    def intersection(self, other):
        if self.m == other.m:
            return None
        inter = (other.c - self.c) / (self.m - other.m)
        if self.domain.contains(inter) and other.domain.contains(inter):
            return inter
        return None

    # comparator: if comparator(self, f2) > 0 then self is more favourable than f2
    # where f2 is a list of Affine objects
    # returns a list of Affine objects representing the piecewise minimum
    #
    # TODO: we could make this function more efficient if we can assume
    # that f2 is already a piecewise-defined function, i.e. no two Affine objects
    # in f2 overlap in domain. Maybe a future optimisation to do?
    def compare_with(self, f2, domain, comparator):

        def compare_f(p, q, domain):
            if p is None and q is None:
                return None
            if q is None:
                f = p
            elif p is None:
                f = q
            else:
                x = domain.midpoint()
                f = p if comparator(p.at(x), q.at(x)) else q
            return Affine(f.m, f.c, f.domain.intersect(domain), f.label)

        # Construct a unique sorted list of critical points
        crits = set({self.domain.x0, self.domain.x1})
        crits.update(p.domain.x0 for p in f2)
        crits.update(p.domain.x1 for p in f2)

        # Add intersections
        for p in f2:
            inter = self.intersection(p)
            if inter is not None:
                crits.add(inter)
        for i in range(len(f2)):
            for j in range(i):
                inter = f2[i].intersection(f2[j])
                if inter is not None:
                    crits.add(inter)

        # crop to domain
        if domain is not None:
            crits = set(b for b in crits if domain.contains(b))
            crits.update([domain.x0, domain.x1])

        crits = list(crits)
        crits.sort()

        # Given boundary points [b0, b1, ... bn], iterate through sets of
        # the form {b0} (b0, b1), {b1}, (b1, b2), ... , {bn} and compute the
        # maximum on each set, as an Affine, storing them in the pieces list
        pieces = []
        i = 0
        processed_boundary = False
        while (i < len(crits) - 1) or not processed_boundary:
            if not processed_boundary:
                # Consider the singleton interval [x, x]
                x = crits[i]
                interval = Interval(x, x, True, True)
                if domain is None or domain.contains(x):
                    p = self if self.domain.contains(x) else None
                    for f in f2:
                        if f.domain.contains(x):
                            p = compare_f(p, f, interval)
                    if p is not None:
                        pieces.append(p)
                processed_boundary = True
                continue

            # Consider the interval (crits[i], crits[i + 1])
            interval = Interval(crits[i], crits[i + 1], False, False)
            test_x = interval.midpoint()
            p = self if self.domain.contains(test_x) else None
            for f in f2:
                if f.domain.contains(test_x):
                    p = compare_f(p, f, interval)
            if p is not None:
                pieces.append(p)
            processed_boundary = False
            i += 1

        # Simplify the pieces, which are arranged in increasing order of domain
        j = 0
        while j < len(pieces) - 1:
            # Determine whether piece j and piece j + 1 can be combined into a
            # single piece
            p1 = pieces[j]
            p2 = pieces[j + 1]
            match = (
                p1.domain.x1 == p2.domain.x0
                and (p1.domain.include_upper or p2.domain.include_lower)
                and pieces[j].label == pieces[j + 1].label
                and (
                    (
                        p1.domain.length() == 0
                        and p1.at(p1.domain.x1) == p2.at(p2.domain.x0)
                    )
                    or (
                        p2.domain.length() == 0
                        and p1.at(p1.domain.x1) == p2.at(p2.domain.x0)
                    )
                    or (p1.function_equals(p2))
                )
            )

            # Extend the interval of p1 to include p2 (no need to update argmin[j])
            if match:
                pieces.pop(j + 1)
                p1.domain.x1 = p2.domain.x1
                p1.domain.include_upper = p2.domain.include_upper
            else:
                j += 1

        return pieces

    # Returns a list of Affine objects representing the minimum between this
    # function and a list of other Affine functions, on a given domain. If a
    # function is not defined on any subset of the domain, it is treated as +\infty.
    # If domain is not specified, it is taken to be the union of the domains of the two
    # functions.
    def min_with(self, other, domain=None):
        return self.compare_with(other, domain, lambda x, y: x < y)

    # Returns a list of Affine objects representing the maximum between this
    # function and a list of other Affine functions, on a given domain. If a
    # function is not defined on any subset of the domain, it is treated as -\infty.
    # If domain is not specified, it is taken to be the union of the domains of the two
    # functions.
    def max_with(self, other, domain=None):
        return self.compare_with(other, domain, lambda x, y: x > y)


# Represents a d-dimensional affine functions
class Affine2:

    # Represents the function f(x) = a^Tx, where a, x are (d + 1)-dimensional
    # vectors. The last d entries of a are the d coefficients of the function, and
    # the first entry of a is the constant term.
    # domain is a d-dimensional polytope representing the domain on which the
    # function is defined.
    def __init__(self, a, domain, label=None):
        if not isinstance(domain, Polytope):
            raise ValueError("domain must be of type Polytope")
        self.a = a
        self.domain = domain
        self.label = label

    # Shallow copy of self
    def __copy__(self):
        return Affine2(copy.copy(self.a), copy.copy(self.domain), copy.copy(self.label))

    # By default if the dimension of the domain is <= 4, then the variable names
    # are (in order): x, y, z, w. Higher dimensional domains take the variable
    # names x_1, ..., x_n.
    def __repr__(self):
        # Get symbols
        ch = list("xyzw")
        if len(self.a) > 5:
            ch = [f"x_{i}" for i in range(1, len(self.a))]
        return Affine2.to_string(self.a, ch) + f" on {str(self.domain)}"

    # Returns an expression given a list of coefficients and a string of symbols
    # (in order of appearance)
    def to_string(coeffs, symbols):
        if symbols is None:
            raise ValueError("symbols")
        if len(symbols) < len(coeffs) - 1:
            raise ValueError("Not enough symbols")

        f = f"{coeffs[0]}"
        for i in range(1, len(coeffs)):
            if coeffs[i] == 0:
                continue
            if coeffs[i] == 1:
                f += " + " + symbols[i - 1]
            elif coeffs[i] == -1:
                f += " - " + symbols[i - 1]
            elif coeffs[i] < 0:
                f += " - " + str(abs(coeffs[i])) + symbols[i - 1]
            else:
                f += " + " + str(coeffs[i]) + symbols[i - 1]

        return f"{f}"

    # Evaluates the function at a particular point in R^d. x must be a list or tuple
    def at(self, x):
        if len(x) != len(self.a) - 1:
            raise ValueError
        return self.a[0] + sum(self.a[i + 1] * x[i] for i in range(len(x)))

    # Computes the intersection between this affine object and another, as a hyperplane
    def intersection(self, other):
        if len(self.a) != len(other.a):
            raise ValueError("Dimensions do not match")
        return Hyperplane([self.a[i] - other.a[i] for i in range(len(self.a))])

    # Returns new Affine2 function with dimension i scaled by a factor
    def scale(self, i, factor, additional_constraints):
        return Affine2(
            self.a.copy(), self.domain.scale(i, factor, additional_constraints)
        )

    # Computes the minimum between this function and a list of 2-d affine functions
    # on the intersection of the domains of the two functions. Everything outside
    # the intersection of domains is ignored (even if one of the functions is
    # defined on it)
    def min_with(self, other):
        pieces = []
        for f in other:

            # Compute the intersection of the functions' domains
            domain = self.domain.intersect(f.domain)

            if domain.is_empty(include_boundary=False):
                continue

            # Check if the functions are equal. If they are equal, then no domain
            # divisions are necessary.
            # If two functions are identical, then we add the one with a smaller
            # label. This is to establish a replicable order of functions.
            if all(self.a[i] == f.a[i] for i in range(len(self.a))):
                fn = self if self.label < f.label else f
                pieces.append(Affine2(fn.a, domain, fn.label))  # Shallow copying
                continue

            # If the functions are not identically equal, check if/where they
            # intersect. Compute the line of intersection as a hyperplane, then
            # check if the hyperplane intersects the domain.
            plane = Hyperplane([self.a[i] - f.a[i] for i in range(len(self.a))])

            if domain.intersects(plane):
                # Representing the constraint (a - a')x + (b - b')y + (c - c') > 0
                d1 = domain.intersect(
                    Polytope([[self.a[i] - f.a[i] for i in range(len(self.a))]])
                )
                d2 = domain.intersect(
                    Polytope([[f.a[i] - self.a[i] for i in range(len(self.a))]])
                )

                # If (a - a')x + (b - b')y + (c - c') > 0 then this > other
                if not d1.is_empty(include_boundary=False):
                    pieces.append(Affine2(list(f.a), d1, f.label))

                if not d2.is_empty(include_boundary=False):
                    pieces.append(Affine2(list(self.a), d2, self.label))
            else:
                # Otherwise - one of the functions is \leq the other function on
                # the entire domain - so just add that function (we have already
                # checked that the domain is nonempty)
                test_x = domain.get_centroid()
                fn = self if self.at(test_x) < f.at(test_x) else f
                pieces.append(Affine2(list(fn.a), domain, fn.label))
        return pieces


# Represents a piecewise-affine function in d dimensions
class Piecewise:

    def __init__(self, pieces):
        self.pieces = pieces

    def __repr__(self):
        if len(self.pieces) < 5:
            return str(self.pieces)
        return f"[{len(self.pieces)} pieces]"

    def __copy__(self):
        return Piecewise([copy.copy(p) for p in self.pieces])

    # --------------------------------------------------------------------------
    # Check that the piecewise function is not one-to-many, and is defined everywhere
    # on the domain
    def check(self, xlim, ylim):
        for i in range(len(self.pieces)):
            for j in range(0, i):
                p = self.pieces[i]
                q = self.pieces[j]
                if not p.domain.intersect(q.domain).is_empty(include_boundary=False):
                    print("Error: domain overlap")
                    print(p.domain)
                    print(q.domain)
                    return False

        # Check that the function is defined everywhere - by sampling
        N = 100
        for i in range(N):
            for j in range(N):
                x = xlim[0] + (xlim[1] - xlim[0]) * i / N
                y = ylim[0] + (ylim[1] - ylim[0]) * j / N
                if self.at([x, y]) is None:
                    print("Error: undefined at", x, y)
                    return False
        return True

    # For 2-D functions only, plot the domain of the function
    def plot_domain(self, xlim, ylim, resolution=300, title=None, xlabel=None, ylabel=None, variables="xy"):
        series = {}
        pieces = {}
        for p in self.pieces:
            pieces[tuple(p.a)] = p
            series[tuple(p.a)] = [[], []]

        N = resolution
        for i in range(N):
            for j in range(N):
                x = xlim[0] + (xlim[1] - xlim[0]) * i / N
                y = ylim[0] + (ylim[1] - ylim[0]) * j / N

                for p in self.pieces:
                    if p.domain.contains([x, y]):
                        s = series[tuple(p.a)]
                        s[0].append(x)
                        s[1].append(y)

        for key in pieces:
            ser = series[key]
            if len(ser[0]) == 0:
                continue
            plt.plot(series[key][0], series[key][1], label=Affine2.to_string(key, variables))
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3)
        #plt.legend(loc="upper right", ncol=3)
        if title is not None:
            plt.title(title)
        if xlabel is not None:
            plt.xlabel(xlabel)
        if ylabel is not None:
            plt.ylabel(ylabel)
        plt.show()

    # --------------------------------------------------------------------------

    # Crops the domain (in place) of this function to the specified domain
    def crop(self, domain, simplify=True):
        if not isinstance(domain, Polytope):
            raise ValueError("domain must be of type Polytope")

        for p in self.pieces:
            p.domain = p.domain.intersect(domain)

        if simplify:
            for i in range(len(self.pieces) - 1, -1, -1):
                p = self.pieces[i]
                if p.domain.is_empty(include_boundary=False):
                    self.pieces.pop(i)

    # Given the two functions
    # f_1 : D_1 \subseteq R^d \mapsto R
    # f_2 : D_2 \subseteq R^d \mapsto R
    # returns min(f_1, f_2) : D_1 \cap D_2 \mapsto R.
    #
    # The returned function is only defined on the intersection of the domain of
    # this function and that of the other function.
    def min_with(self, other):
        pieces = []
        for p in self.pieces:
            pieces.extend(p.min_with(other.pieces))
        return Piecewise(pieces)

    def at(self, x):
        for p in self.pieces:
            if p.domain.contains(x):
                return p.at(x)
        return None

    # Returns a new Piecewise object with each piece having dimension i scaled by factor
    def scale(self, i, factor, additional_constraints):
        return Piecewise(
            [p.scale(i, factor, additional_constraints) for p in self.pieces]
        )

    # Simplify the function by computing unions of polytopes
    def simplify(self, max_grouping=3):
        
        # a single simplification iteration, which tries to represent multiple
        # Affine2 objects as a single object. The parameter n represents the
        # number of objects we try to combine at a time .
        def iteration(matchable, n):
            
            # Flag to indicate whether any changes have been made
            changed = False
            
            # Prevent modifying set when iterating
            keys = [k for k in matchable]
            for key in keys:
                # a list of compatible Affine2 objects (same function, same label)
                group = matchable[key]
                for c in itertools.combinations(range(len(group)), n):
                    union = Polytope.try_union([group[fi].domain for fi in c])
                    if union is not None:
                        # Add new element at the end of list 
                        f = group[c[0]]
                        new_f = Affine2(f.a, union, f.label)
                        group.append(new_f)
                        
                        # Remove indices of c from group
                        matchable[key] = [group[gi] for gi in range(len(group)) if gi not in c]
                        
                        # Since the group has now been altered, we can no longer
                        # continue iterating through this group. Instead move
                        # onto the next group.
                        changed = True
                        
                        break
            return changed
        
        # First - compute the matchable sets (in terms of the function and label)
        matchable = {}
        for i in range(len(self.pieces)):
            p = self.pieces[i]
            # Create a tuple representing the function and label, to act as a hash
            key = tuple([p.label] + p.a)
            if key in matchable:
                matchable[key].append(p)
            else:
                matchable[key] = [p]
        
        # Often it is possible to aggregate the entire set - try that first
        # Unfortunately, this method is simply too slow
        # keys = [k for k in matchable]
        # for key in keys:
        #     group = matchable[key]
        #     if len(group) <= 1: continue
        #     union = Polytope.try_union([f.domain for f in group])
        #     if union is not None:
        #         f = group[0]
        #         matchable[key] = [Affine2(f.a, union, f.label)]
            
        # For now - only match 2 - 3 at a time
        for n in range(2, max_grouping + 1):
            while iteration(matchable, n):
                pass

        # Finally: re-group the matchable and set self.pieces
        pieces = []
        for key in matchable:
            pieces.extend(matchable[key])
        
        self.pieces = pieces


# Auxiliary functions for sympy
class SympyHelper:

    # Handles conversions between Sympy object types and Fraction, if possible
    def to_frac(expr):
        if isinstance(expr, sympy.core.numbers.Zero):
            return frac(0)
        if isinstance(expr, sympy.core.numbers.Half):
            return frac(1, 2)
        if isinstance(expr, sympy.core.numbers.One):
            return frac(1)
        if isinstance(expr, sympy.core.numbers.Rational):
            n, d = sympy.fraction(expr)
            if isinstance(n, sympy.core.numbers.Integer) and isinstance(
                n, sympy.core.numbers.Integer
            ):
                return frac(int(n), int(d))
        if isinstance(expr, sympy.core.numbers.Float):
            return float(expr)
        return expr

# Alternative implementation of the univariate function P(x)/Q(x) where P, Q are polynomials
# This class is currently experimental.
# Differences to RationalFunction class:
# 1) Interval is used to represent domain of definition
# 2) Less reliance on symbolic algebra - interval storage uses tuple of numbers instead of 
#   symbolic object
# 3) Object is now immutable 
class RationalFunction2:
    
    # coeffs: either: 
    # - list of numbers, representing the coefficients of P(x)
    # - tuple of list of numbers, representing the coefficients of P(x), Q(x)
    # domain: Interval
    def __init__(self, coeffs, domain):
        if not isinstance(coeffs, list) and not isinstance(coeffs, tuple):
            raise ValueError("Parameter coeffs must be of type list of tuple or tuple")
        if not isinstance(domain, Interval):
            raise ValueError("Parameter domain must be of type Interval")
        
        if len(coeffs) == 0:
            raise ValueError("Parameter coeffs must have at least 1 element")

        if len(coeffs) == 1:
            self._num_coeffs = tuple(coeffs)
            self._den_coeffs = tuple([1])
        else:
            self._num_coeffs = tuple(coeffs[0])
            self._den_coeffs = tuple(coeffs[1])
        self.domain = domain
    
    def __repr__(self):

        # handle special cases first 
        if self.is_zero(): 
            return f"0 on {self.domain}"
        if self.is_constant(): 
            return f"{self._num_coeffs[0] / self._den_coeffs[0]} on {self.domain}"

        num = RationalFunction2._poly_str(self._num_coeffs, "x")
        den = RationalFunction2._poly_str(self._den_coeffs, "x")

        if den == "1":
            return f"{num} on {self.domain}"
 
        # Add brackets if required 
        if sum(1 for c in self._num_coeffs if c != 0) > 1: num = "(" + num + ")"
        if sum(1 for c in self._den_coeffs if c != 0) > 1: den = "(" + den + ")"

        return f"{num}/{den} on {self.domain}"

    # Static functions ------------------------------------------------------------
    # Pretty-print a polynomial given the coefficients and the variable name 
    def _poly_str(coeffs, var_name):
        deg = len(coeffs) - 1
        s = []
        front = True
        for c in coeffs:
            if front: 
                if c < 0: s.append("-")
            else:
                s.append("-" if c < 0 else "+")
            s.append(RationalFunction2._mon_str(c, var_name, deg))
            deg -= 1
            front = False

        if len(s) == 0: return "0"
        return " ".join(s)

    # Pretty-print a monomial term |c|x^{deg}
    def _mon_str(c, var, deg):
        if c == 0: return ""
        if deg == 0: return f"{abs(c)}"
        stem = f"{var}" if deg == 1 else f"{var}^{deg}"
        if abs(c) == 1: return stem
        return f"{abs(c)}" + stem
    
    # Instance methods -----------------------------------------------------------
    # Returns whether this polynomial is zero
    def is_zero(self):
        return all(c == 0 for c in self._num_coeffs)

    def is_constant(self):
        return all(c == 0 for c in self._num_coeffs[1:]) and \
                all(c == 0 for c in self._den_coeffs[1:])

    def is_affine(self):
        return all(c == 0 for c in self._num_coeffs[2:]) and \
                all(c == 0 for c in self._den_coeffs[1:])

    
# Represents the univariate function P(x)/Q(x) where P, Q are polynomials. Currently
# this is just a wrapper around the sympy symbolic algebra library.
class RationalFunction:

    # Static variable to be re-used
    x = sympy.symbols("x")

    # Create a rational function by providing the coefficients (arranged in order
    # from highest degree to lowest degree) for both the numerator and denominator
    def __init__(self, num_coeffs, den_coeffs=None):

        if den_coeffs == None:
            den_coeffs = [1]

        num = 0
        xn = 1
        for coeff in num_coeffs[::-1]:
            num += coeff * xn
            xn *= RationalFunction.x

        den = 0
        xn = 1
        for coeff in den_coeffs[::-1]:
            den += coeff * xn
            xn *= RationalFunction.x
        self.num = num  # the numerator
        self.den = den  # the denominator

    def __repr__(self):
        return f"{sympy.latex(sympy.simplify(self.num / self.den))}"

    def __eq__(self, other):
        return (self.num * other.den - self.den * other.num).equals(0)


    # Static functions -------------------------------------------------

    def parse(expr):
        s = sympy.parsing.sympy_parser.parse_expr(expr)
        (num, den) = sympy.fraction(s)
        r = RationalFunction([])
        r.num = num
        r.den = den
        return r

    # ------------------------------------------------------------------

    def at(self, x):
        f = self.num / self.den
        if isinstance(f, numbers.Number):
            return f
        return f.subs(RationalFunction.x, x)

    def add(self, other):
        if isinstance(other, numbers.Number):
            r = RationalFunction([])
            r.num = self.num + self.den * other
            r.den = self.den
            return r

        # Hacky initialisation
        r = RationalFunction([])
        r.num = self.num * other.den + self.den * other.num
        r.den = self.den * other.den
        return r

    def mul(self, other):
        if isinstance(other, numbers.Number):
            # Hacky initialisation
            r = RationalFunction([])
            r.num = self.num * other
            r.den = self.den
            return r
        raise NotImplementedError

    # Computes this function divided by another
    def div(self, other):
        # Hacky initialisation
        r = RationalFunction([1])
        r.num = self.num * other.den
        r.den = self.den * other.num
        return r

    # find all critical points
    # TODO: cache the critical points
    def _compute_critical_pts(self, interval):
        # Find all critical points
        crit = [interval.x0, interval.x1]

        # Stationary points of P(x) / Q(x) occur when P(x)Q'(x) - P'(x)Q(x) = 0
        dP = sympy.diff(self.num, RationalFunction.x)
        dQ = sympy.diff(self.den, RationalFunction.x)
        df = self.num * dQ - dP * self.den

        if not df.is_constant():
            crit.extend(
                SympyHelper.to_frac(r)
                for r in sympy.real_roots(df)
                if interval.contains(r)
            )

        return crit

    # Computes all the intersection points between this RationalFunction and another
    # on an interval
    def intersections(self, other, interval):
        # Intersections of P/Q = R/S when PS - RQ = 0
        f = self.num * other.den - self.den * other.num
        if not f.is_constant():
            return [
                SympyHelper.to_frac(r)
                for r in sympy.real_roots(f)
                if interval.contains(r)
            ]
        return []

    def roots(self, interval):
        f = self.num / self.den
        if f.is_constant():
            return []
        return [r for r in sympy.real_roots(f) if interval.contains(r)]

    # Compute the maximum over an interval of type Interval
    def maximise(self, interval):
        crit = self._compute_critical_pts(interval)
        i = max(enumerate(crit), key=lambda u: self.at(u[1]))[0]
        return (crit[i], self.at(crit[i]))

    # Compute the minimum over an interval of type Interval
    def minimise(self, interval):
        crit = self._compute_critical_pts(interval)
        i = min(enumerate(crit), key=lambda u: self.at(u[1]))[0]
        return (crit[i], self.at(crit[i]))
