import cdd
import copy
from fractions import Fraction as frac
import itertools


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
        ch_set = list('xyzw')
        if len(self.coefficients) > 5:
            ch_set = [f'x_{i}' for i in range(1, len(self.coefficients) - 1)]
            
        f = ''
        if self.coefficients[0] != 0:
            f += str(self.coefficients[0])
        for i in range(1, len(self.coefficients)):
            v = ch_set[i - 1]
            c = self.coefficients[i]
            if c == 0: continue
            if len(f) != 0:
                f += (' + ' if c > 0 else ' - ')
            if c == 1 or c == -1:
                f += v
            else:
                f += str(abs(c)) + v
            
        if self.inequality_type == Constraint.GREATER_EQUALS:
            return f'{f} >= 0'
        if self.inequality_type == Constraint.GREATER:
            return f'{f} > 0'
        return f'{f} = 0'
      
    def __eq__(self, other):
        if len(self.coefficients) != len(other.coefficients) or \
            self.inequality_type != other.inequality_type:
            return False
        return all(self.coefficients[i] == other.coefficients[i] for i in range(len(self.coefficients)))
        
    # Returns whether the point x (dimension d vector) satisfies the constraint
    def satisfies(self, x):
        d = len(self.coefficients) - 1
        if len(x) != d:
            raise ValueError(f'x must have length {d}')
            
        q = sum(self.coefficients[i + 1] * x[i] for i in range(0, d)) + self.coefficients[0]
       
        if self.inequality_type == Constraint.GREATER_EQUALS: return q >= 0
        elif self.inequality_type == Constraint.GREATER: return q > 0
        elif self.inequality_type == Constraint.EQUALS: return q == 0
    
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
    
    # build Polytope object in terms of its H representation
    # Parameters:
    #   - constraints: a list of lists, each of which has the form [a_0 a_1 ... a_d]
    #           representing the inequality a_0 + a_1 x_1 + ... + a_d x_d \geq 0
    def __init__(self, constraints, linear=False, canonicalize=False):
        self.mat = cdd.Matrix(constraints, linear=linear, number_type='fraction')
        self.mat.rep_type = cdd.RepType.INEQUALITY
        if canonicalize:
            self.mat.canonicalize()
        self.polyhedron = cdd.Polyhedron(self.mat)
        
        # Workspaces for V representation
        self.vertices = None
        self.rays = None
    
    def __copy__(self):
        p = Polytope([]) # Construct an empty polytope temporarily
        p.mat = self.mat.copy()
        p.polyhedron = cdd.Polyhedron(p.mat)
        p.vertices = None
        p.rays = None
        return p
    
    def __repr__(self):
        def to_str(row, is_lin):
            ch_set = list('xyzw')
            if len(row) > 5:
                ch_set = [f'x_{i}' for i in range(1, len(row) - 1)]
                
            f = ''
            if row[0] != 0:
                f += str(row[0])
            for i in range(1, len(row)):
                v = ch_set[i - 1]
                c = row[i]
                if c == 0: continue
                if len(f) != 0:
                    f += (' + ' if c > 0 else ' - ')
                if c == 1 or c == -1:
                    f += v
                else:
                    f += str(abs(c)) + v
            if is_lin:
                return f'{f} = 0'
            return f'{f} >= 0'
        
        return str([to_str(self.mat[i], i in self.mat.lin_set) for i in range(len(self.mat))])
    
    def __eq__(self, other):
        if self.vertices is None:
            self.compute_V_rep()
        if other.vertices is None:
            other.compute_V_rep()
        # Construct set of vertices from V representation (as tuples which are hashable)
        return set(tuple(v) for v in self.vertices) == set(tuple(v) for v in other.vertices)
    
    # -------------------------------------------------------------------------
    # internal functions

    # Check if the point x satisfies the inequality or equality given in row 
    # row_index of the given matrix
    def _satisfies(matrix, row_index, x):
        r = matrix[row_index]
        s = r[0] + sum(r[i + 1] * x[i] for i in range(len(x)))
        return (s == 0) if row_index in matrix.lin_set else (s >= 0)


    # -------------------------------------------------------------------------
    # Static public functions 
    
    # Compute the joint union with multiple polytopes, if the resulting union is 
    # a polytope (otherwise, this function returns None). 
    #
    # This implementation is based on the union algorithm in 
    # A. Bemporad et al. "Convexity recognition of the union of polyhedra" (2001) 
    # Algorithm 4.1
    def try_union(polys):
        if not isinstance(polys, list) or len(polys) == 0:
            return None
        
        if len(polys) == 1: return polys[0]
        
        env = [] # Stores the constraints used to compute the envelope
        neg = [] # Stores the negation of the removed constraints (1 list per polytope)
        
        # Classifies the constraint by checking if the constraint in row 'row_index'
        # of 'matrix' is satisfied by all polytopes
        # except for polytope i
        def satisfies(matrix, row_index, i):
            for j in range(len(polys)):
                if i != j: 
                    if not all(Polytope._satisfies(matrix, row_index, v) \
                               for v in polys[j].get_vertices()):
                        return False
            return True
            
        # Iterate through the polynomials and decide whether each constraint is 
        # of env type or neg type
        for i in range(len(polys)):
            p = polys[i]
            neg_p = []
            for r in range(p.mat.row_size):
                if satisfies(p.mat, r, i):
                    env.append(p.mat[r])
                else:
                    # Add the negation of the constraint
                    neg_p.append([-x for x in p.mat[r]])
            neg.append(neg_p)
        
        # Choose 1 constraint from each negated constraint set, and
        # compute the intersection between them and env
        for combination in itertools.product(*neg):
            # Check polytope formed (equivalent to the feasibility LP). If 
            # the polytope formed is non-empty, then A U B < Env, so the 
            # resulting union is not a polytope. Returns None in this case
            p = Polytope(env + [c for c in combination])
            if not p.is_empty(include_boundary=False):
                return None
        
        # Otherwise, envelop is the union
        return Polytope(env, canonicalize=True)
        
    # Constructs a rectangle as a polytope
    def rect(xlim, ylim):
        return Polytope([
                [-xlim[0], 1, 0],     # x >= xlim[0]
                [xlim[1], -1, 0],     # x <= xlim[1]
                [-ylim[0], 0, 1],     # y >= ylim[0]
                [ylim[1], 0, -1]      # y <= ylim[1]
            ])
    
    # Create a polytope from its vertices 
    # TODO: add support for extreme rays
    def from_V_rep(verts):
        V = cdd.Matrix([[1] + v for v in verts], number_type='fraction')
        V.rep_type = cdd.RepType.GENERATOR
        poly = cdd.Polyhedron(V)
        H = poly.get_inequalities()
        
        p = Polytope([]) # Construct an empty polytope temporarily
        p.mat = H
        p.polyhedron = poly
        p.vertices = copy.copy(verts) # Shallow copy is fine for now
        p.rays = None
        return p
        
    # -------------------------------------------------------------------------
    # public instance functions 

    # Computes the V representation of this polytope
    def compute_V_rep(self):
        v = self.polyhedron.get_generators()
        self.vertices = [r[1:] for r in v if r[0] == 1] # unpack matrix output
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
        return [Constraint(list(r), Constraint.GREATER_EQUALS) for r in self.mat]
    
    # Returns a new matrix with the ith dimension scaled by a factor 
    def scale(self, i, factor, additional_constraints):
        
        ineq_rows = [list(self.mat[r]) for r in range(self.mat.row_size) if r not in self.mat.lin_set]
        for c in ineq_rows:
            c[i] /= factor
        ineq_rows.extend(additional_constraints)
                
        eq_rows = [list(self.mat[r]) for r in range(self.mat.row_size) if r in self.mat.lin_set]
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
        
    # Returns whether this region contains the point x
    def contains(self, x):
        if len(x) != self.dimension():
            raise ValueError(f'Dimension of x must be {self.dimension()}')
        for row in self.mat:
            if row[0] + sum(row[i + 1] * x[i] for i in range(len(x))) < 0:
                return False
        return True
    
        
    def dimension(self):
        return self.mat.col_size - 1
    
    def num_constraints(self):
        return self.mat.row_size
    
    # Computes the intersection of this polytope with another
    def intersect(self, other):
        if not isinstance(other, Polytope):
            raise ValueError('Parameter other must be of type Polytope')
        mat = self.mat.copy()
        
        # TODO: check lin_set? Currently passes unit tests
        rows = [r for r in other.mat]
        if len(rows) > 0:
            mat.extend(rows, linear=False)
            
        return Polytope(mat, canonicalize=True)
    
    # Computes A \ B where A is this polytope and B is another polytope
    def set_minus(self, other):
        if not isinstance(other, Polytope):
            raise ValueError('Parameter other must be of type Polytope')
        # A \ B = A intersect (B complement)
        mat = self.mat.copy()
        mat.extend([[-x for x in r] for r in other.mat if r not in other.lin_set], linear=False)
        # ignore linear constraint since we do not consider boundaries
        # mat.extend([r for r in other.mat if r in other.lin_set], linear=True)
        return Polytope(mat, canonicalize=True)
        
    # Returns whether this polytope intersects with a plane. Parameter plane must 
    # be of type Hyperplane
    def intersects(self, plane):
        if not isinstance(plane, Hyperplane):
            raise ValueError('Parameter plane must be of type Hyperplane')
        
        if self.vertices is None:
            self.compute_V_rep()
        
        # Check if all vertices are on the same side of the hyperplane
        sign = 0
        a = plane.constraint.coefficients
        for vert in self.vertices:
            new_sign = a[0] + sum(a[i + 1] * vert[i] for i in range(len(vert)))
            if sign * new_sign < 0:
                return True # points are on opposite sides of the plane
            if new_sign != 0:
                sign = new_sign
        return False
    
    # Returns whether the polytope is empty, i.e. whether the constraints are consistent
    # This is determined by checking the number of vertices in its V representation
    def is_empty(self, include_boundary=True):
        if self.vertices is None:
            self.compute_V_rep()
        
        if include_boundary:
            return len(self.vertices) == 0
        
        # If not including boundary, strict subspace vertex regions are considered empty
        return len(self.vertices) <= self.dimension()
        