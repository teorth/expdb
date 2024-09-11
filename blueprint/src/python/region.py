
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
    def _try_simplify_union_of_polytopes(region, max_groupsize=3):
        if region.region_type not in {Region_Type.UNION, Region_Type.DISJOINT_UNION} or \
            not all(r.region_type == Region_Type.POLYTOPE for r in region.child):
            return False
        
        polys = [r.child for r in region.child]
        
        # a single simplification iteration, which tries to represent multiple
        # Affine2 objects as a single object. The parameter groupsizen represents 
        # the number of objects we try to combine at a time.
        def iteration(objs, groupsize):
            for c in itertools.combinations(range(len(objs)), groupsize):
                union = Polytope.try_union([objs[i] for i in c])
                if union is not None:
                    # Remove indices of c from group, add new element at the end of list 
                    return [objs[i] for i in range(len(objs)) if i not in c] + [union]
            return None
        
        changed = False

        for n in range(2, max_groupsize + 1):
            while True:
                new_polys = iteration(polys, n)
                if new_polys is None:
                    break
                polys = new_polys
                changed = True

        region.child = [Region(Region_Type.POLYTOPE, p) for p in polys]  # repack
            
    
    
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
            
            # Sets of lists of polytopes
            child_sets = [r._as_disjoint_union_poly() for r in self.child]
            disjoint = child_sets[0]
            for i in range(1, len(child_sets)):
                new_disjoint = []
                for p in disjoint:
                    for q in child_sets[i]:
                        x = (0.5526333971071384, 2.6127523920876934, 1.4688895492136673, 1.2858300913524678, 3.1622832001931718)
                        before = p.contains(x) and q.contains(x)

                        inter = p.intersect(q)

                        after = inter.contains(x)
                        if before != after:
                            print(before, after)
                            print("components")
                            print(p)
                            print(q)
                            print("intersect")
                            print(inter)
                            raise ValueError()
                        if not inter.is_empty(include_boundary=False):
                            new_disjoint.append(inter)
                print(len(new_disjoint))
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
        changed = changed or Region._try_simplify_union_of_polytopes(self, 2)

        return changed
    

