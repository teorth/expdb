
import copy
import itertools
import matplotlib.pyplot as plt
import numpy as np
from polytope import Polytope

class Region_Type:
    
    # Supports many regions 
    INTERSECT = 0
    UNION = 1
    DISJOINT_UNION = 2
    
    # Supports one region or polytope
    POLYTOPE = 3
    COMPLEMENT = 4


    
    def to_str(region_type):
        if region_type == Region_Type.INTERSECT:
            return "Intersection"
        elif region_type == Region_Type.UNION:
            return "Union"
        elif region_type == Region_Type.DISJOINT_UNION:
            return "Disjoint union"
        elif region_type == Region_Type.POLYTOPE:
            return "Polytope"
        elif region_type == Region_Type.COMPLEMENT:
            return "Complement"
        raise ValueError("Unknown type", region_type)
    
    
    
# A region represents a boolean combination of polytopes 
class Region:
    def __init__(self, region_type, children):
        if not isinstance(region_type, int):
            print(region_type)
            raise ValueError("Parameter region_type must be of type int.")
        if not isinstance(children, Polytope) and \
            not isinstance(children, list):
            raise ValueError("Parameter polytope must either be of type Polytope or list of Polytope")
        
        self.region_type = region_type
        self.child = children

    def __repr__(self):
        return self.to_str(0)

    def __copy__(self):
        return Region(self.region_type, copy.copy(self.child))
        
    def to_str(region, indentation=0):
        s = ("\t" * indentation) + Region_Type.to_str(region.region_type) + "\n"
        if isinstance(region.child, Polytope):
            s += ("\t" * indentation) + str(region.child) + "\n"
        else:
            for r in region.child:
                s += Region.to_str(r, indentation + 1)
        return s
    
    # Private methods --------------------------------------------------------

    # If the region is a union of polytopes, try to simplify it
    # Returns True if changes were made
    def _try_simplify_union_of_polys(region):

        if region.region_type not in {Region_Type.UNION, Region_Type.DISJOINT_UNION} or \
            not all(r.region_type == Region_Type.POLYTOPE for r in region.child):
            return False
        
        polys = [r.child for r in region.child]
        polys = Region._simplify_union_of_polys(polys)
        region.child = [Region(Region_Type.POLYTOPE, p) for p in polys]  # repack

        return True

    # Simplify a union of polytopes 
    def _simplify_union_of_polys(polys):
        
        # The list of polytopes will be expanded, without removal
        # Instead we keep track of the indexes of the polytopes included in the union
        # in this set 
        index = set(i for i in range(len(polys)))

        # Keep track of the tuples (int, int) of indexes which do not form a pairwise 
        # union
        non_unions = set()
        changed = True
        while changed:

            changed = False
            for c in itertools.combinations(list(index), 2):
                (i, j) = c
                # First check if i, j have not been removed from a previous iteration
                if i not in index or j not in index:
                    continue

                # Skip checking the expensive try_union operation if the combination
                # is known to be non-unionable
                if (i, j) in non_unions or (j, i) in non_unions: 
                    continue

                union = Polytope.try_union([polys[i], polys[j]])
                if union is not None:
                    # Union is successful: remove indexes in c and add new index
                    index.remove(i)
                    index.remove(j)
                    index.add(len(polys))
                    polys.append(union)
                    changed = True
                else:
                    non_unions.add((i, j))

        return [polys[i] for i in index]
    
    # Static methods ---------------------------------------------------------

    # Construct a region from a polytope
    def from_polytope(polytope):
        return Region(Region_Type.POLYTOPE, polytope)

    # Construct a region as the complement of another region 
    def complement(region):
        return Region(Region_Type.COMPLEMENT, region)

    # Compute the union of regions 
    def union(regions):
        return Region(Region_Type.UNION, regions)

    # Compute the union of regions, assuming that they are disjoint 
    def disjoint_union(regions):
        return Region(Region_Type.DISJOINT_UNION, regions)

    # Compute the intersection of regions
    def intersect(regions):
        return Region(Region_Type.INTERSECT, regions)
    
    # Instance methods -------------------------------------------------------

    # Returns whether this region contains the point x
    def contains(self, x):
        if self.region_type == Region_Type.COMPLEMENT:
            return not self.child.contains(x)
        if self.region_type == Region_Type.POLYTOPE:
            return self.child.contains(x)
        if self.region_type == Region_Type.UNION or self.region_type == Region_Type.DISJOINT_UNION:
            return any(c.contains(x) for c in self.child)
        if self.region_type == Region_Type.INTERSECT:
            return all(c.contains(x) for c in self.child)
        raise NotImplementedError(self.region_type)

    # Returns a new region object where the substitute function is called on each polytope
    def substitute(self, values):
        if self.region_type == Region_Type.POLYTOPE:
            return Region(Region_Type.POLYTOPE, self.child.substitute(values))
        # Handle regions with single child
        if self.region_type == Region_Type.COMPLEMENT: 
            return Region(Region_Type.COMPLEMENT, self.child.substitute(values))
        # Handle regions with multiple children
        return Region(self.region_type, [c.substitute(values) for c in self.child])
    
    # Given a list of numbers of length dimension(), returns a new region where each 
    # dimension is scaled by the given factor 
    def scale_all(self, factors):
        if not isinstance(factors, list):
            raise ValueError("Parameter factors must be of type list")
        # Handle regions with single child
        if self.region_type in {Region_Type.POLYTOPE, Region_Type.COMPLEMENT}:
            return Region(self.region_type, self.child.scale_all(factors))
        # Handle regions with multiple children
        return Region(self.region_type, [c.scale_all(factors) for c in self.child])
    
    # Project this region onto a set of dimensions (set of integers)
    # Here we make the assumption that the projection of a union (resp. intersection) 
    # of sets is the union (resp. intersection) of the projection of each set. 
    def project(self, dims):
        if self.region_type == Region_Type.POLYTOPE:
            return Region(Region_Type.POLYTOPE, self.child.project(dims))
        
        # Unfortunately the disjointness of unions are not preserved under projection
        if self.region_type == Region_Type.UNION or self.region_type == Region_Type.DISJOINT_UNION:
            return Region(Region_Type.UNION, [c.project(dims) for c in self.child])
        
        if self.region_type == Region_Type.INTERSECTION:
            return Region(Region_Type.INTERSECTION, [c.project(dims) for c in self.child])
        
        # In particular, complements are not yet supported. 
        raise NotImplementedError()
    
    # Returns a representation of this region as a disjoint union of convex polytopes.
    # Returns the result as a Region object of type DISJOINT_UNION
    def as_disjoint_union(self):
        polys = self._as_disjoint_union_poly()
        return Region.disjoint_union([Region(Region_Type.POLYTOPE, p) for p in polys])

    # Internal method 
    # Same as the as_disjoint_union(self) method except it returns the result as a 
    # list of Polytope objects.
    def _as_disjoint_union_poly(self):
        if self.region_type == Region_Type.COMPLEMENT:
            raise NotImplementedError(self.region_type) # TODO: implement this
        if self.region_type == Region_Type.POLYTOPE:
            return [self.child]
        if self.region_type == Region_Type.INTERSECT:
            # Compute intersection of regions using distributive property of set 
            # algebra, i.e. 
            # A ∩ (B1 u B2 u ... u Bn) = (A ∩ B1) u (A ∩ B2) u ... u (A ∩ Bn)
            SIMPLIFY_EVERY = 10
        
            # Sets of lists of polytopes
            child_sets = [r._as_disjoint_union_poly() for r in self.child]
            disjoint = child_sets[0]
            for i in range(1, len(child_sets)):
                new_disjoint = []
                for p in disjoint:
                    for q in child_sets[i]:
                        inter = p.intersect(q)
                        if not inter.is_empty(include_boundary=False):
                            new_disjoint.append(inter)
                
                # Every few rounds, simplify
                if i % SIMPLIFY_EVERY == 0:
                    prevlen = len(new_disjoint)
                    new_disjoint = Region._simplify_union_of_polys(new_disjoint)
                    print(prevlen, "->", len(new_disjoint))

                print(len(new_disjoint), i + 1, "of", len(child_sets))
                disjoint = new_disjoint

            return disjoint
        if self.region_type == Region_Type.UNION:
            raise NotImplementedError(self.region_type) # TODO: implement this
        if self.region_type == Region_Type.DISJOINT_UNION:
            result = []
            for r in self.child:
                result.extend(r._as_disjoint_union_poly())
            return result
        raise NotImplementedError(self.region_type)
    
    # If this region is 2-dimensional, plot it (see also the plot method in Polytope)
    def plot2d(self, xlim, ylim, resolution=100):
        xvals = np.linspace(xlim[0], xlim[1], resolution)
        yvals = np.linspace(ylim[0], ylim[1], resolution)
        xs = []
        ys = []
        for x in xvals:
            for y in yvals:
                if self.contains((x, y)):
                    xs.append(x)
                    ys.append(y)
        
        plt.plot(xs, ys)
        plt.show()

    # Simplify the Region
    # TODO: add more simplification routines 
    def simplify(self):

        changed = False
        # If this region is a union of polytopes, it may be simplified
        changed = changed or Region._try_simplify_union_of_polys(self)

        return changed
    
