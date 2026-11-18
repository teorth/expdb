from functions import RationalFunction2

# Basic univariate dense polynomial implementation with rational coefficients
class Polynomial:

    # Initialize with coefficients as a list, ordered from 
    # lowest degree to highest degree. If len(coefficients) = 0
    # then the polynomial is initialized to 0. 
    def __init__(self, coefficients):
        if coefficients is None or not isinstance(coefficients, tuple):
            raise ValueError("Parameter coefficients must be of type tuple")

        # leading coefficient must be non-zero
        if len(coefficients) > 0 and coefficients[-1] == 0:
            raise ValueError("Leading coefficient must be non-zero")
        
        self.coefficients = coefficients

    def __repr__(self):
        return self.to_str("x")

    def __copy__(self):
        return Polynomial(self.coefficients)

    def __eq__(self, other):
        return self.coefficients == other.coefficients

    def __add__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.subtract(other)

    def __mul__(self, other):
        return self.multiply(other)

    # Returns the string representation of this polynomial, with a specified 
    # variable name. 
    def to_str(self, var):
        deg = 0
        s = []
        front = True
        for c in self.coefficients:
            if front: 
                if c < 0: s.append("-")
            else:
                s.append("-" if c < 0 else "+")

            mon = Polynomial._monomial_str(c, var, deg)
            if len(mon) > 0:
                front = False
                s.append(mon)
            deg += 1

        if len(s) == 0: return "0"
        return " ".join(s)

    # Returns a string representation of a monomial 
    def _monomial_str(c, var, deg):
        if c == 0: return ""
        if deg == 0: return f"{abs(c)}"
        stem = f"{var}" if deg == 1 else f"{var}^{deg}"
        if abs(c) == 1: return stem
        return f"{abs(c)}" + stem

    # Removes the last elements of a tuple, if they are zero
    def _trim_trailing_zeroes(c):
        n = len(c)
        while n > 0:
            if c[n - 1] != 0:
                return c if n == len(c) else c[:n]
            n -= 1
        return tuple() # empty tuple - all elements are zero

    def degree(self):
        n = len(self.coefficients)
        if n == 0: return float('-inf')
        return n - 1

    def is_zero(self):
        return len(self.coefficients) == 0

    # Returns a new polynomial equal to the sum of this polynomial and another
    def add(self, other):
        if not isinstance(other, Polynomial):
            raise ValueError("Parameter other must be of type Polynomial.")
        p, q = self.coefficients, other.coefficients
        m, n = len(p), len(q)
        if m == n:
            c = tuple(p[i] + q[i] for i in range(m))
        elif m > n:
            c = tuple(p[i] + q[i] for i in range(n)) + p[n:m]
        else:
            c = tuple(p[i] + q[i] for i in range(m)) + q[m:n]
        return Polynomial(Polynomial._trim_trailing_zeroes(c))
    
    # Returns a new polynomial equal to this polynomial subtract another
    def subtract(self, other):
        if not isinstance(other, Polynomial):
            raise ValueError("Parameter other must be of type Polynomial.")
        p, q = self.coefficients, other.coefficients
        m, n = len(p), len(q)
        if m == n:
            c = tuple(p[i] - q[i] for i in range(m))
        elif m > n:
            c = tuple(p[i] - q[i] for i in range(n)) + p[n:m]
        else:
            c = tuple(p[i] - q[i] for i in range(m)) + tuple(-a for a in q[m:n])
        return Polynomial(Polynomial._trim_trailing_zeroes(c))
    
    # Returns a new polynomial equal to this polynomial multiplied by another
    def multiply(self, other):
        if not isinstance(other, Polynomial):
            raise ValueError("Parameter other must be of type Polynomial.")
        if self.is_zero() or other.is_zero():
            return Polynomial(())

        p, q = self.coefficients, other.coefficients
        m, n = len(p), len(q)
        c = [0] * (m + n - 1)
        for i in range(m):
            for j in range(n):
                c[i + j] += p[i] * q[j]
        return Polynomial(Polynomial._trim_trailing_zeroes(tuple(c)))

    # Returns this polynomial divided by another 
    def divide(self, other):
        if not isinstance(other, Polynomial):
            raise ValueError("Parameter other must be of type Polynomial.")
        if other.is_zero():
            raise ValueError("Division by zero")

        # Returns 
        return RationalFunction2(self, other)