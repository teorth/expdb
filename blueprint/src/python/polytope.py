import cdd
import copy
from fractions import Fraction as frac
import itertools
import matplotlib.pyplot as plt
import numpy as np



###############################################################################
# Represents a d-dimensional linear constraint of the form (a^T)x > 0
# Deprecated
class Constraint:

    # inequality types
    GREATER_EQUALS = 0
    GREATER = 1
    EQUALS = 2

    # Creates an object representing the linear constraint c + a^Tx \geq 0
    #
    # Parameters:
    #   - coefficients: (Number list) a d + 1 dimensional vector containing [c a]
    #   - inequality_type:  0: \geq, 1: >, 2: =
    def __init__(self, coefficients, inequality_type=GREATER_EQUALS):
        self.coefficients = coefficients
        self.inequality_type = inequality_type

    def __repr__(self):

        # Get symbols
        ch_set = list("xyzwu")
        if len(self.coefficients) > 5:
            ch_set = [f"x_{i}" for i in range(1, len(self.coefficients) - 1)]

        f = ""
        if self.coefficients[0] != 0:
            f += str(self.coefficients[0])
        for i in range(1, len(self.coefficients)):
            v = ch_set[i - 1]
            c = self.coefficients[i]
            if c == 0:
                continue
            if len(f) != 0:
                f += " + " if c > 0 else " - "
            if c == 1 or c == -1:
                f += v
            else:
                f += str(abs(c)) + v

        if self.inequality_type == Constraint.GREATER_EQUALS:
            return f"{f} >= 0"
        if self.inequality_type == Constraint.GREATER:
            return f"{f} > 0"
        return f"{f} = 0"

    def __eq__(self, other):
        if (
            len(self.coefficients) != len(other.coefficients)
            or self.inequality_type != other.inequality_type
        ):
            return False
        return all(
            self.coefficients[i] == other.coefficients[i]
            for i in range(len(self.coefficients))
        )

    def __hash__(self):
        # Borrow the hashing function of tuples
        return hash(tuple(self.coefficients) + (self.inequality_type,))

    # Returns whether the point x (dimension d vector) satisfies the constraint
    def satisfies(self, x):
        d = len(self.coefficients) - 1
        if len(x) != d:
            raise ValueError(f"x must have length {d}")

        q = (
            sum(self.coefficients[i + 1] * x[i] for i in range(0, d))
            + self.coefficients[0]
        )

        if self.inequality_type == Constraint.GREATER_EQUALS:
            return q >= 0
        elif self.inequality_type == Constraint.GREATER:
            return q > 0
        elif self.inequality_type == Constraint.EQUALS:
            return q == 0


###############################################################################
# Represents a d-dimensional hyperplane residing in (d + 1)-dimensional Euclidean space
class Hyperplane:
    # Define a hyperplane as a[0] + a[1:] \cdot x = 0 where a is a (d + 1)-dimensional
    # vector (represented as a list)
    def __init__(self, a):
        self.constraint = Constraint(a, inequality_type=Constraint.EQUALS)

    def __repr__(self):
        return str(self.constraint)

    # Determine whether a point x lies on this plane
    def contains(self, x):
        return self.constraint.satisfies(x)


###############################################################################
# Represents a d-dimensional convex polytope
class Polytope:

    # Build Polytope object in terms of its H representation
    # Parameters:
    #   - constraints: a list of lists, each of which has the form [a_0 a_1 ... a_d]
    #           representing the inequality a_0 + a_1 x_1 + ... + a_d x_d \geq 0
    #   - linear: (boolean) if True, then the constraints are equalities, otherwise 
    #           the constraints are inequalities
    #   - canonicalize: (boolean) if True, then the list of constraints will be 
    #           canonicalized. In general this operation is expensive, however certain
    #           operations that depending on a matrix's lin_set require it.   
    def __init__(self, constraints, linear=False, canonicalize=False):
        
        # There is a strange bug when constraint is a list of zeroes - where 
        # the resulting matrix adds to the linset. Remove such redundant constraints
        # before constructing the matrix 
        constraints = [c for c in constraints if not all(ci == 0 for ci in c)]
        
        self.mat = cdd.Matrix(constraints, linear=linear, number_type="fraction")
        self.mat.rep_type = cdd.RepType.INEQUALITY
        if canonicalize:
            self.mat.canonicalize()

        # Cache whether this polytope is guaranteed to be canonical
        self.is_canonical = canonicalize
        # Cache whether this polytope is empty
        self._is_empty_incl_boundary = None 
        self._is_empty_excl_boundary = None 
        
        self.polyhedron = cdd.Polyhedron(self.mat)

        # Workspaces for V representation
        self.vertices = None
        self.rays = None

    def __copy__(self):
        p = Polytope._from_mat(self.mat.copy())
        p.is_canonical = self.is_canonical
        return p

    def __repr__(self):
        if self.dimension() <= 5:
            return self.to_str("xyzwu") # Default variable names
        return self.to_str([f"x_{i}" for i in range(1, self.dimension() + 1)])
        
    def to_str(self, variables):

        if len(variables) < self.dimension():
            raise ValueError(f"Parameter variables must be a list of length >= {self.dimension()}.")
        
        def _to_str(row, is_equality):
            f = ""
            if row[0] != 0:
                f += str(row[0])
            for i in range(1, len(row)):
                v = variables[i - 1]
                c = row[i]
                if c == 0:
                    continue
                elif c > 0:
                    if len(f) != 0: f += " + "
                else:
                    f += " - " if len(f) != 0 else "-"
                if c == 1 or c == -1:
                    f += v
                else:
                    f += str(abs(c)) + v
            return f"{f} = 0" if is_equality else f"{f} >= 0"

        return str(
            [_to_str(self.mat[i], i in self.mat.lin_set) for i in range(len(self.mat))]
        )

    def __eq__(self, other):
        if self.vertices is None:
            self.compute_V_rep()
        if other.vertices is None:
            other.compute_V_rep()
        # Construct set of vertices from V representation (as tuples which are hashable)
        return set(tuple(v) for v in self.vertices) == set(
            tuple(v) for v in other.vertices
        )

    # -------------------------------------------------------------------------
    # internal/private functions

    # Check if the point x satisfies the inequality or equality given in row
    # row_index of the given matrix
    def _satisfies(matrix, row_index, x):
        r = matrix[row_index]
        s = r[0] + sum(r[i + 1] * x[i] for i in range(len(x)))
        return (s == 0) if row_index in matrix.lin_set else (s >= 0)
    
    # Construct a Polytope object directly from a matrix instead of from constraints
    def _from_mat(matrix):
        p = Polytope([])  # Construct an empty polytope temporarily
        p.mat = matrix
        p.polyhedron = cdd.Polyhedron(matrix)
        p.vertices = None
        p.rays = None
        return p
    
    # Returns the constraints defining this polytope as a list of tuples
    # If lin_set is True, then only return equalities. If lin_set is False, then
    # only return inequalities 
    def _matrix_as_list(matrix, lin_set):
        return [
            matrix[r]
            for r in range(matrix.row_size)
            if (r in matrix.lin_set) == lin_set
        ]
    
    # -------------------------------------------------------------------------
    # Static public functions

    # Compute the intersection of a list of convex polytopes, returning the 
    # result as a convex polytope. 
    def intersection(polys, canonicalize=True):
        if not isinstance(polys, list) or not all(isinstance(p, Polytope) for p in polys):
            raise ValueError("Parameter polys must be a list of Polytope.")
        
        if len(polys) == 0:
            raise ValueError("Must have at least one polytope.")

        if len(polys) == 1:
            return polys[0]
        
        # Copy the constraints of the first polytope 
        mat = polys[0].mat.copy()
        
        # Loop through the remaining polytopes and add their constraints
        ineqs = []  # the inequality constraints from the other polytopes
        eqs = []    # the equality constraints from the other polytopes
        for i in range(1, len(polys)):
            ineqs.extend(Polytope._matrix_as_list(polys[i].mat, False))
            eqs.extend(Polytope._matrix_as_list(polys[i].mat, True))
        
        if len(ineqs) > 0: mat.extend(ineqs, linear=False)
        if len(eqs) > 0: mat.extend(eqs, linear=True)

        if canonicalize:
            mat.canonicalize()
        
        p = Polytope._from_mat(mat)
        p.is_canonical = canonicalize
        return p

    # Compute the union of a list of polytopes, if the resulting union is
    # a polytope (otherwise, this function returns None).
    #
    # This implementation is based on the union algorithm in
    # A. Bemporad et al. "Convexity recognition of the union of polyhedra" (2001)
    # Algorithm 4.1
    #
    # This implementation only works with finite polytopes - for infinite polytopes 
    # the wrong result is returned. See e.g. the test cases in test_polytope. 
    # TODO: implement check and error for infinite polytopes. 
    def try_union(polys):
        if not isinstance(polys, list) or len(polys) == 0:
            return None

        if len(polys) == 1:
            return polys[0]

        # Performance optimisation: see if intersection is nonempty first. This 
        # is relative cheap compared to the full union operation, since it only
        # requires solving 1 LP (without canonicalisation, since boundaries are
        # included) and is expected to occur quite frequently in practice, 
        # especially for high dimensional polytopes. However, this only works 
        # when taking the union of 2 polytopes 
        if len(polys) == 2 and \
            Polytope.intersection(polys, canonicalize=False).is_empty(include_boundary=True):
            return None
        
        env = set()  # Stores the constraints used to compute the envelope
        neg = []  # Stores the negation of the removed constraints (1 list per polytope)

        # Classifies the constraint by checking if the constraint in row 'row_index'
        # of 'matrix' is satisfied by all polytopes except for polytope i
        def satisfies(matrix, row_index, i):
            for j in range(len(polys)):
                if i != j:
                    if not all(
                        Polytope._satisfies(matrix, row_index, v)
                        for v in polys[j].get_vertices()
                    ):
                        return False
            return True

        # Iterate through the polytopes and decide whether each constraint is
        # of env type or neg type
        for i in range(len(polys)):
            p = polys[i]
            neg_p = []
            for r in range(p.mat.row_size):
                if satisfies(p.mat, r, i):
                    # Use hashability of tuples to ensure uniqueness in set
                    env.add(tuple(p.mat[r])) 
                    # If this a linear constraint, we also need to add its negation
                    if r in p.mat.lin_set:
                        env.add(tuple(-x for x in p.mat[r]))
                else:
                    # Add the negation of the constraint
                    neg_p.append([-x for x in p.mat[r]])
                    if r in p.mat.lin_set:
                        neg_p.append([x for x in p.mat[r]])
                
            neg.append(neg_p)

        # Choose 1 constraint from each negated constraint set, and
        # compute the intersection between them and env
        env = list(env)
        for combination in itertools.product(*neg):
            # Check polytope formed (equivalent to the feasibility LP). If
            # the polytope formed is non-empty, then A U B < Env, so the
            # resulting union is not a polytope. Returns None in this case
            p = Polytope(env + [c for c in combination])
            if not p.is_empty(include_boundary=False):
                return None
        
        # Otherwise, envelop is the union
        return Polytope(env, canonicalize=True)

    # Constructs a d-dimensional box as a polytope
    # lims is a ordered list of limits for each variable
    def rect(*lims):
        bounds = []
        dim = len(lims)
        for i in range(dim):
            lim = lims[i]
            
            b = [-lim[0]] + ([0] * dim)
            b[i + 1] = 1
            bounds.append(b)   # x >= lim[0]
            
            b = [lim[1]] + ([0] * dim)
            b[i + 1] = -1
            bounds.append(b)   # x <= lim[1]

        return Polytope(bounds)

    # Create a polytope from its vertices
    # TODO: add support for extreme rays
    def from_V_rep(verts):

        # cdd does not support empty polytopes, so we return None
        if len(verts) == 0:
            return None

        V = cdd.Matrix([[1] + v for v in verts], number_type="fraction")
        V.rep_type = cdd.RepType.GENERATOR
        poly = cdd.Polyhedron(V)
        H = poly.get_inequalities()

        p = Polytope([])  # Construct an empty polytope temporarily
        p.mat = H
        p.polyhedron = poly
        p.vertices = copy.copy(verts)  # Shallow copy is fine for now
        p.rays = []
        return p

    # -------------------------------------------------------------------------
    # public instance functions

    # Set the internal matrix to canonical form
    def canonicalize(self):
        self.mat.canonicalize()
        self.is_canonical = True

    # The dimension of the space in which this polytope resides 
    # (this is not the rank of the H-representation matrix) 
    def dimension(self):
        return self.mat.col_size - 1

    # Returns true if this polytope cannot be embedded into a lower-dimension
    def is_full_dim(self):
        if not self.is_canonical:
            self.mat.canonicalize()
            self.is_canonical = True
        return len(self.mat.lin_set) == 0
        
    def num_constraints(self):
        return self.mat.row_size

    # Computes the V representation of this polytope
    def compute_V_rep(self):
        v = self.polyhedron.get_generators()
        self.vertices = [r[1:] for r in v if r[0] == 1]  # unpack matrix output
        self.rays = [r[1:] for r in v if r[0] == 0]

    # Returns the vertices of this polytope
    def get_vertices(self):
        if self.vertices is None:
            self.compute_V_rep()
        return self.vertices

    # Returns the centeriod of this polytope
    def get_centroid(self):

        # Check if V-representation is cached
        if self.vertices is None:
            self.compute_V_rep()

        cent = [0] * self.dimension()
        for v in self.vertices:
            for i in range(len(v)):
                cent[i] += v[i]

        den = frac(len(self.vertices), 1)
        for i in range(self.dimension()):
            cent[i] /= den
        return cent

    # Returns this polytope as a set of Constraint objects
    def get_constraints(self):
          # The inequality constraints
        return [Constraint(list(r), Constraint.GREATER_EQUALS) for r in Polytope._matrix_as_list(self.mat, False)] \
                + [Constraint(list(r), Constraint.EQUALS) for r in Polytope._matrix_as_list(self.mat, True)]

    # Returns the edges of this polytope, returned as 2-tuples of vertices
    def get_edges(self):
        if not self.is_canonical:
            self.canonicalize() 
            # Need to reset polyhedron and recompute vertices
            self.polyhedron = cdd.Polyhedron(self.mat)
            self.compute_V_rep()
        
        # Check if V rep is already computed 
        if self.vertices is None:
            self.compute_V_rep()
        
        # Get adjacency list
        adj = self.polyhedron.get_adjacency()

        edges = []
        for i in range(len(adj)):
            v1 = self.vertices[i]
            for j in adj[i]:
                v2 = self.vertices[j]
                edges.append((v1, v2))
        return edges

    # Returns a new polytope with the ith dimension scaled by a factor
    # Note that i is 1-indexed
    def scale(self, i, factor, additional_constraints):

        ineq_rows = [list(r) for r in Polytope._matrix_as_list(self.mat, False)]
        for c in ineq_rows:
            c[i] /= factor
        ineq_rows.extend(additional_constraints)

        eq_rows = [list(r) for r in Polytope._matrix_as_list(self.mat, True)]
        for c in eq_rows:
            c[i] /= factor

        # The matrix will not have dimensions initialised if the number of inequality
        # rows = 0, hence this strange duplicated code
        if len(ineq_rows) > 0:
            p = Polytope(ineq_rows, linear=False, canonicalize=False)
            if len(eq_rows) > 0:
                p.mat.extend(eq_rows, linear=True)

        elif len(eq_rows) > 0:
            p = Polytope(eq_rows, linear=True, canonicalize=False)

        p.mat.canonicalize()

        # reset polytope object
        p.polyhedron = cdd.Polyhedron(p.mat)
        return p
    
    # Given a list of length self.dimension(), scales the i-th dimension by factors[i]
    def scale_all(self, factors):
        if not isinstance(factors, list):
            raise ValueError("Parameter factors must be of type list")
        if len(factors) != self.dimension():
            raise ValueError(f"Parameter factors must have length {self.dimension()}")
        
        D = len(factors)
        ineq_rows = Polytope._matrix_as_list(self.mat, False)
        eq_rows = Polytope._matrix_as_list(self.mat, True)

        mat = None
        if len(ineq_rows) > 0:
            ineq_rows = [[r[0]] + [r[i + 1] / factors[i] for i in range(D)] for r in ineq_rows]
            mat = cdd.Matrix(ineq_rows, linear=False, number_type="fraction")
            mat.rep_type = cdd.RepType.INEQUALITY
        
        if len(eq_rows) > 0:
            eq_rows = [[r[0]] + [r[i + 1] / factors[i] for i in range(D)] for r in eq_rows]
            if mat is None:
                mat = cdd.Matrix(eq_rows, linear=True, number_type="fraction")
                mat.rep_type = cdd.RepType.INEQUALITY
            else:
                mat.extend(eq_rows, linear=True)
        
        return Polytope._from_mat(mat)

    # Returns whether this region contains the point x
    def contains(self, x):
        if len(x) != self.dimension():
            raise ValueError(f"Dimension of x must be {self.dimension()}")
            
        for ri in range(self.mat.row_size):
            row = self.mat[ri]
            q = row[0] + sum(row[i + 1] * x[i] for i in range(len(x)))
            if ri in self.mat.lin_set:
                if q != 0: return False
            else:
                if q < 0: return False
        return True

    # Returns whether this polytope intersects with a plane. Parameter plane must
    # be of type Hyperplane
    def intersects(self, plane):
        if not isinstance(plane, Hyperplane):
            raise ValueError("Parameter plane must be of type Hyperplane")

        if self.vertices is None:
            self.compute_V_rep()

        # Check if all vertices are on the same side of the hyperplane
        sign = 0
        a = plane.constraint.coefficients
        for vert in self.vertices:
            new_sign = a[0] + sum(a[i + 1] * vert[i] for i in range(len(vert)))
            if sign * new_sign < 0:
                return True  # points are on opposite sides of the plane
            if new_sign != 0:
                sign = new_sign
        return False

    # Computes the intersection of this polytope with another. This polytope is 
    # unchanged 
    def intersect(self, other):
        return Polytope.intersection([self, other], canonicalize=True)

    # Computes A \ B where A is this polytope and B is another polytope. Returns 
    # a list of polytopes that are guaranteed to be disjoint, and whose union 
    # equals A \ B. Note that since this polytope implementation ignores boundaries
    # we ignore any possible ambiguity regarding whether points lying on the 
    # boundary of the resulting polytopes belong to A \ B. 
    def set_minus(self, other, simplify=True):
        if not isinstance(other, Polytope):
            raise ValueError("Parameter other must be of type Polytope")
        
        # Iterate through the constraints of B. For each constraint, take its 
        # complement, then add the original constraint to the set of constraints
        # for the new region. All linear constraints are ignored (see above for 
        # rationale)
        A = self.mat.copy()
        B = other.mat
        polys = []
        for i in range(B.row_size):
            if i in B.lin_set: continue
            c = B[i]
            # Invert the constraint c then compute A \intersect c'
            c_comp = [-x for x in c]
            A_copy = A.copy()
            A_copy.extend([c_comp], linear=False)
            if simplify:
                A_copy.canonicalize()

            part = Polytope._from_mat(A_copy)
            if simplify:
                if not part.is_empty(include_boundary=False):
                    polys.append(part)
            else:
                polys.append(part)

            # For the next part
            if i < B.row_size - 1:
                A.extend([c])

        return polys

    # Returns whether the polytope is empty, i.e. whether the constraints are consistent
    # if include_boundary is True, the boundary of the polytope is considered part of 
    # the polytope
    def is_empty(self, include_boundary=True):

        if include_boundary:
            # Case 1: ------------------------------------------------------------------
            if self._is_empty_incl_boundary is not None:
                return self._is_empty_incl_boundary

            # If the polytope is closed, then one can check for emptiness by solving the 
            # linear program Ax <= b with an arbitrary objective function to check for 
            # feasibility. Here we use the objective function (0, 0, ..., 0)
            # Note: this routine no longer passes unit tests if objective function is set 
            # to e.g. (1, 1, ..., 1) - possibly because of different lp.status values?
            A = self.mat.copy()
            A.obj_type = cdd.LPObjType.MAX
            A.obj_func = tuple([0] * (1 + self.dimension()))
            lp = cdd.LinProg(A)
            lp.solve()

            # If status is "OPTIMAL", then a solution was found for this LP and 
            # the region represented by this polytope is non-empty
            self._is_empty_incl_boundary = (lp.status != cdd.LPStatusType.OPTIMAL)
            return self._is_empty_incl_boundary

        else:
            # Case 2: ------------------------------------------------------------------
            if self._is_empty_excl_boundary is not None:
                return self._is_empty_excl_boundary

            # If we do not consider the boundary, we need to check the feasibility of the 
            # system Ax < b (with strict inequality). To do this we solve the LP 
            # 
            # max epsilon s.t. Ax + epsilon <= b, 0 <= epsilon <= 1
            # 
            # If the solution is epsilon = 0 then the polytope (excluding boundaries) is 
            # empty. 

            # Edge case: if there are any linear constraints, the polytope is empty
            if len(self.mat.lin_set) > 0:
                return True

            # coefficient matrix represented as [M 1] . [x epsilon]' where 1 is a column 
            # of 1's and M is the set of inequality constraints
            A = cdd.Matrix(
                [list(row) + [-1] for row in Polytope._matrix_as_list(self.mat, False)] + 
                [[0] * self.mat.col_size + [1]] +               # epsilon >= 0
                [[1] + [0] * (self.mat.col_size - 1) + [-1]],   # epsilon <= 1
                linear=False
            )
            A.obj_type = cdd.LPObjType.MAX
            A.obj_func = tuple([0] * (self.dimension() + 1) + [1]) # objective is epsilon
            lp = cdd.LinProg(A)
            lp.solve()

            # If no solution for epsilon, the region is empty
            if lp.status != cdd.LPStatusType.OPTIMAL:
                return True

            # If epsilon = 0, then the feasible region = the boundary of the region, 
            # i.e. there is no interior
            if lp.primal_solution[-1] == 0:
                return True
            return False    

    # Returns whether this polytope is a subset of another.
    def is_subset_of(self, other):
        if not isinstance(other, Polytope):
            raise ValueError("Parameter other must be of type Polytope")
        
        # Returns the H rep of a polytope as inequalities only (representing 
        # equality constraints via two-sided constraints)
        def as_ineq(mat):
            ineq = Polytope._matrix_as_list(mat, False)
            ineq.extend(Polytope._matrix_as_list(mat, True))
            ineq.extend([tuple(-x for x in r) 
                for r in Polytope._matrix_as_list(mat, True)])
            return ineq 
        
        Ab = as_ineq(other.mat) # Represent other polytope as Ax <= b

        # Construct matrix
        mat = self.mat.copy()
        mat.obj_type = cdd.LPObjType.MIN

        # Iterate through the rows of Ab
        for ab in Ab:
            # If there is a point x in P1 that is not in P2 then P1 is not a 
            # subset of P2. 
            # If a point is in P2 it must satisfy each constraint of the form 
            # ab >= 0. If we cannot find such a point within P1, then it is not 
            # in P2, so P1 is not a subset of P2.  
            mat.obj_func = ab
            lp = cdd.LinProg(mat)
            lp.solve()
            if lp.obj_value < 0:
                return False
        return True

    # Returns true if this polytope is contained in the union of a list of polytopes
    # This method implements the Algorithm 3 in Mato Baoti\'{c}
    # "Polytopic Computations in Constrained Optimal Control" (2009) 
    def is_covered_by(self, polys):
        if not isinstance(polys, list):
            raise ValueError("Parameter polys must be of type list.")
        
        # Given a list of polytopes Q1, ..., Qn, with q1, ..., qn constraints
        # respectively, return False if there exist n constraints c1, ..., cn
        # where ci belongs to qi, such that the region P intersect neg c1, ..., cn
        # is nonempty (not including boundaries)
        # Strategy is to try and find this region using DFS. 

        # Terminal condition
        if len(polys) == 0:
            return False
        
        # Convert the first polytope Q into a list of inequality constraints 
        # (convert equality constraints to inequality constraints)
        Q = polys[0]
        cons = Polytope._matrix_as_list(Q.mat, False)
        for eq in Polytope._matrix_as_list(Q.mat, True):
            cons.append(eq)
            cons.append(tuple(-x for x in eq))

        # Iterate through the constraints of Q, looking for violations
        for con in cons:
            # Intersect P with the negation of a constraint
            mat = self.mat.copy()
            mat.extend([[-x for x in con]], linear=False)
            poly = Polytope._from_mat(mat)

            # If there are points in P not in Q, check whether 
            # these points are also not in the other Qi's
            if not poly.is_empty(include_boundary=False):
                if not poly.is_covered_by(polys[1:]):
                    return False
        
        return True

    # Given a dictionary "values" of the form {i:v} where i is a non-negative integer and
    # v is a Number, compute a polytope formed by taking i-th variable as v in this Polytope. 
    def substitute(self, values):

        # The rows of the shrunk matrix with the specified dimensions projected out
        new_rows = []
        for r in self.mat:
            new_row = [r[0] + sum(r[i + 1] * values[i] for i in values)]
            new_row.extend(r[i] for i in range(1, len(r)) if i - 1 not in values)
            new_rows.append(new_row)
        
        new_mat = cdd.Matrix(new_rows)
        new_mat.rep_type = cdd.RepType.INEQUALITY
        new_mat.lin_set = self.mat.lin_set # copy over the lin set
        return Polytope._from_mat(new_mat)

    # Given a set of integers representing the indices of dimensions, returns a 
    # new Polytope object projected onto those dimensions.
    #
    # Note: this method uses the V representation of a polytope so only works 
    # for finite polytopes. 
    def project(self, dims):
        if not isinstance(dims, set):
            raise ValueError("Parameter dims must of type set.")
        
        if self.vertices is None:
            self.compute_V_rep()
        
        projected_verts = []
        for v in self.vertices:
            projected_verts.append([v[i] for i in range(len(v)) if i in dims])
        proj = Polytope.from_V_rep(projected_verts)
        return proj

    # Given a list (var) x[0], x[1], ..., x[N - 1], with N > self.dimension(), 
    # returns a polytope formed by lifting this polytope to N dimensions. 
    # Each x[i] is either
    # - an integer, in which case the i-th dimension of the resulting polytope 
    #   is the x[i]-th dimension of this polytope, or
    # - a tuple (a, b) of numbers representing the limits along the (new) i-th 
    #   dimension of the resulting polytope
    def lift(self, var):
        if not isinstance(var, list):
            raise ValueError("Parameter var must be of type list")
        if not all(isinstance(x, int) or isinstance(x, tuple) for x in var):
            raise ValueError("Each entry of parameter var must be either of type int or tuple")

        N = len(var)

        # The inequality constraints of the form ax >= 0
        ineq_constraints = [self.mat[i] for i in range(self.mat.row_size) if i not in self.mat.lin_set] 
        # the equality constraints of the form ax = 0
        eq_constraints = [self.mat[i] for i in range(self.mat.row_size) if i in self.mat.lin_set]

        # Convert existing constraints into one in N dimensions
        def f(constraint):
            # Account for the constant term in the front
            return [constraint[0]] + [(constraint[x + 1] if isinstance(x, int) else 0) for x in var]

        ineq_constraints = [f(c) for c in ineq_constraints]
        eq_constraints = [f(c) for c in eq_constraints]

        # Add constraints due to additional variables 
        for i in range(N):
            if isinstance(var[i], int): 
                continue
            (lower, upper) = var[i]
            # The lower bound
            c = [-lower] + ([0] * N)
            c[i + 1] = 1
            ineq_constraints.append(c)

            # The upper bound 
            c = [upper] + ([0] * N)
            c[i + 1] = -1
            ineq_constraints.append(c)
        
        # Construct matrix 
        mat = cdd.Matrix(ineq_constraints, linear=False, number_type="fraction")
        mat.rep_type = cdd.RepType.INEQUALITY
        if len(eq_constraints) > 0:
            mat.extend(eq_constraints, linear=True)
        return Polytope._from_mat(mat)


    # Plot the polytope (current only 2D polytopes are supported)
    def plot(self, resolution=100):
        if self.dimension() != 2:
            raise NotImplementedError("Currently, only 2-dimensional polytopes are supported")
        
        # Compute the bounds 
        if self.vertices is None:
            self.compute_V_rep()
        
        (xmin, xmax) = min(v[0] for v in self.vertices), max(v[0] for v in self.vertices)
        (ymin, ymax) = min(v[1] for v in self.vertices), max(v[1] for v in self.vertices)

        xvals = np.linspace(xmin, xmax, resolution)
        yvals = np.linspace(ymin, ymax, resolution)

        xs = []
        ys = []
        for x in xvals:
            for y in yvals:
                if self.contains([x, y]):
                    xs.append(x)
                    ys.append(y)
        
        plt.plot(xs, ys)
        plt.show()

        
