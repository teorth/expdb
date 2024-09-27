
import copy
import matplotlib.pyplot as plt
import numpy as np
from helpers.region_helper import Region_Helper
from polytope import Polytope
import time

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
        if isinstance(region.child, Polytope):
            s = ("\t" * indentation) + str(region.child) + "\n"
        else:
            s = ("\t" * indentation) + Region_Type.to_str(region.region_type) + "\n"
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
        polys = Region_Helper.simplify_union_of_polys(polys, Polytope.try_union)
        region.child = [Region(Region_Type.POLYTOPE, p) for p in polys]  # repack

        return True

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

    # Returns whether this region is a subset of another region
    # TODO: implement other cases of this method
    def is_subset_of(self, region) -> bool:
        if not isinstance(region, Polytope) and not isinstance(region, Region):
            raise ValueError("Parameter region must be of type Polytope or Region")
        
        # Currently: only polytope and union types are supported
        if self.region_type == Region_Type.POLYTOPE:
            if isinstance(region, Polytope):
                return self.child.is_subset_of(region)
            
            if region.region_type == Region_Type.POLYTOPE:
                return self.child.is_subset_of(region.child)
            elif region.region_type == Region_Type.UNION or region.region_type == Region_Type.DISJOINT_UNION:
                return self.child.is_covered_by(region.child)
            elif region.region_type == Region_Type.INTERSECT:
                return all(self.is_subset_of(p) for p in region.child)
        
        elif self.region_type == Region_Type.UNION or self.region_type == Region_Type.DISJOINT_UNION:
            return all(r.is_subset_of(region) for r in self.child)
        
        # Everything else is not implemented yet
        raise NotImplementedError()

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
    def project(self, dims: set[int]) -> 'Region':
        if not isinstance(dims, set):
            raise ValueError("Parameter dims must be of type set[int]")
        
        if self.region_type == Region_Type.POLYTOPE:
            p = self.child.project(dims)
            if p is None: return None
            return Region(Region_Type.POLYTOPE, p)
        
        # Unfortunately the disjointness of unions are not preserved under projection
        if self.region_type == Region_Type.UNION or self.region_type == Region_Type.DISJOINT_UNION:
            ps = [c.project(dims) for c in self.child]
            ps = [p for p in ps if p is not None]
            if len(ps) == 0: return None
            if len(ps) == 1: return ps[0]
            return Region(Region_Type.UNION, ps)
        
        if self.region_type == Region_Type.INTERSECT:
            ps = [c.project(dims) for c in self.child]
            ps = [p for p in ps if p is not None]
            if len(ps) == 0: return None
            if len(ps) == 1: return ps[1]
            return Region(Region_Type.INTERSECT, ps)
        
        # In particular, complements are not yet supported. 
        raise NotImplementedError()
    
    # See Polytope.lift(var)
    def lift(self, var) -> 'Region':
        
        # Need to check carefully that all operations are preserved under lifting
        if self.region_type not in {Region_Type.POLYTOPE, Region_Type.UNION, 
                                    Region_Type.DISJOINT_UNION, Region_Type.INTERSECT, Region_Type.COMPLEMENT}:
            raise NotImplementedError("lifting operation is not implemented for the region_type {self.region_type}.")
        
        if self.region_type in {Region_Type.POLYTOPE, Region_Type.COMPLEMENT}:
            return Region(Region_Type.POLYTOPE, self.child.lift(var))
        
        return Region(self.region_type, [c.lift(var) for c in self.child])

    # Returns a representation of this region as a disjoint union of convex polytopes.
    # Returns the result as a Region object of type DISJOINT_UNION
    def as_disjoint_union(self) -> 'Region':
        polys = self._as_disjoint_union_poly()
        return Region.disjoint_union([Region(Region_Type.POLYTOPE, p) for p in polys])

    # Internal method 
    # Same as the as_disjoint_union(self) method except it returns the result as a 
    # list of Polytope objects.
    def _as_disjoint_union_poly(self) -> list[Polytope]:
        if self.region_type == Region_Type.COMPLEMENT:
            raise NotImplementedError(self.region_type) # TODO: implement this
        if self.region_type == Region_Type.POLYTOPE:
            return [self.child]
        if self.region_type == Region_Type.INTERSECT:

            # Compute intersection of regions using distributive property of set 
            # algebra, i.e. 
            #
            # A ∩ (B1 u B2 u ... u Bn) = (A ∩ B1) u (A ∩ B2) u ... u (A ∩ Bn)
            #
            # This implementation contains some ad-hoc optimizations which improved
            # performance on some test problems. Once every SIMPLIFY_EVERY iterations 
            # the number of polytopes in the union is reduced by combining multiple
            # polytopes into a single polytope (this is an expensive operation but on
            # balance helps with performance). In addition, every time we take the 
            # intersection of a polytope with a union of polytopes, if the number of 
            # nonempty regions is no greater than UNION_THRESHOLD, we try to simplify 
            # the union of the regions into a single polytope. 
            SIMPLIFY_EVERY = 10
            UNION_THRESHOLD = 5

            start_time = time.time()

            # Sets of lists of polytopes
            child_sets = [r._as_disjoint_union_poly() for r in self.child]
            A = child_sets[0]
            i = 0
            for B in child_sets[1:]:
                i += 1

                # # [[ad-hoc performance optimisation]] that did not improve performance 
                # # check if A is a subset of B, so that A intersect B = A and we may
                # # avoid computing the intersection altogether. 
                # if all(p.is_covered_by(B) for p in A):
                #     continue
                
                new_A = []
                for p in A:
                    inters = []
                    for q in B:
                        # [[Ad hoc performance optimisation]] - take care of easy cases
                        # first, which tend to occur frequently in practice
                        if p.is_subset_of(q):
                            inters.append(p)
                        elif q.is_subset_of(p):
                            inters.append(q)
                        else:
                            inter = p.intersect(q)
                            if not inter.is_empty(include_boundary=False):
                                inters.append(inter)

                    # [[Ad hoc performance optimisation]] - frequently in practice p is a 
                    # subset of B, so it is a good candidate for merging
                    if 2 <= len(inters) and len(inters) <= UNION_THRESHOLD:
                        union = Polytope.try_union(inters)
                        if union is not None:
                            new_A.append(union)
                            continue
                    
                    # Otherwise, just add the pieces as-is
                    new_A.extend(inters)

                
                # [[Ad hoc performance optimisation]]
                # Every few rounds, simplify if the number of polytopes is too large
                if i % SIMPLIFY_EVERY == 0 and len(new_A) > 100:
                    prevlen = len(new_A)
                    new_A = Region_Helper.simplify_union_of_polys(new_A, Polytope.try_union)
                    print(prevlen, "->", len(new_A))

                print("Processed", i, "of", len(child_sets) - 1, ":", len(new_A), "pieces, computed in", time.time() - start_time, "sec")
                A = new_A
            return A
        
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
    
