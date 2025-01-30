
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
            not isinstance(children, Region) and \
            not isinstance(children, list):
            raise ValueError("Parameter polytope must either be of type Polytope, Region or list of Polytope/Region")
            
        self.region_type = region_type
        self.child = children

    def __repr__(self):
        return self.to_str(0)

    def __copy__(self):
        return Region(self.region_type, copy.copy(self.child))
        
    def to_str(region, use_indentation=True, indentation=0, variables=None):
        # Returns a long-form expression
        if use_indentation:
            if isinstance(region.child, Polytope):
                p = str(region.child) if variables is None else region.child.to_str(variables)
                return ("\t" * indentation) + p + "\n"
            else:
                s = ("\t" * indentation) + Region_Type.to_str(region.region_type) + "\n"
                for r in region.child:
                    s += Region.to_str(r, use_indentation, indentation + 1, variables)
                return s
            
        # Returns a short form expression
        else:
            if isinstance(region.child, Polytope):
                return str(region.child) if variables is None else region.child.to_str(variables)
            else:
                return Region_Type.to_str(region.region_type) + " of {" + \
                    ", ".join(Region.to_str(r, use_indentation=False, variables=variables) for r in region.child) + "}"
    

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
    
    # Computes a Region object representing the union of a list of polytopes, 
    # where each polytope is defined as the intersection of 
    # - a box, represented as a list of lists (a list of constraints)
    # - a halfplane, represented as a single list (a constraint)
    # The resulting Region object is represented as a DISJOINT_UNION, which 
    # greatly improves performance over a UNION representation
    # 
    # This method is useful for quickly initializing many commonly encountered 
    # regions in the study of large value energy regions, since they correspond to 
    # a region implied by a single max() function. 
    def from_union_of_halfplanes(halfplanes, box):
        # Once a halfplane has been added, include its complement in the list of 
        # neg_ineq, to add as a constraint to all remaining polytopes to be 
        # constructed. 
        neg_ineq = []
        polys = []
        for hp in halfplanes:
            polys.append(
                Region(
                    Region_Type.POLYTOPE, 
                    Polytope(box + [hp] + neg_ineq, canonicalize=True)
                )
            )
            neg_ineq.append([-x for x in hp])
        return Region.disjoint_union(polys)

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

    def set_label(self, label):
        """
        Sets a label for each Polytope in this Region object 
        """
        if self.region_type == Region_Type.POLYTOPE:
            self.child.label = label
            return 
        
        # Recursively set the label attribute
        if isinstance(self.child, Region):
            self.child.set_label(label)
        else:
            for r in self.child:
                r.set_label(label)


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
                # This assumes that the children of region are regions of type Region_Type.POLYTOPE
                if not all(
                        isinstance(r, Region) and r.region_type == Region_Type.POLYTOPE 
                        for r in region.child
                    ):
                    raise NotImplementedError("Parameter region is of a currently unsupported type.")
                return self.child.is_covered_by([r.child for r in region.child])
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

    def as_disjoint_union(self, verbose=False, track_dependencies=False) -> 'Region':
        """
        Computes a representation of this region as a union of convex polytopes. 
        
        Parameters 
        ----------
        verbose : boolean, optional
            If True, debugging information will be logged to the console 
            (default is False)
        track_dependencies : boolean, optional
            If True, each polytope in the output will contain a "dependencies" 
            attribute that contains the "label" attribute of all polytopes used 
            to generate that polytope (default is False).

        Returns
        -------
        Region
            A region object of type "DISJOINT_UNION"
        """
        polys = self._as_disjoint_union_poly(verbose=verbose, track_dependencies=track_dependencies)
        return Region.disjoint_union([Region(Region_Type.POLYTOPE, p) for p in polys])

    def _as_disjoint_union_poly(self, verbose=False, track_dependencies=False) -> list[Polytope]:
        """
        Computes a representation of this region as a union of convex polytopes. 
        
        Parameters 
        ----------
        verbose : boolean, optional
            If True, debugging information will be logged to the console 
            (default is False)
        track_dependencies : boolean, optional
            If True, each polytope in the output will contain a "dependencies" 
            attribute that contains the "label" attribute of all polytopes used 
            to generate that polytope (default is False).

        Raises
        ------
        NotImplementedError
            If the disjoint union operation is not yet implemented for the 
            region type. 
        
        Returns
        -------
        list
            a list of Polytope objects, whose union equals this Region.  
        """

        if self.region_type == Region_Type.COMPLEMENT:
            raise NotImplementedError(self.region_type) # TODO: implement this
        if self.region_type == Region_Type.POLYTOPE:
            if track_dependencies:
                # If tracking labels, create a copy so as to not add the 
                # dependencies attribute to the original polytope object. 
                poly = copy.copy(self.child)
                poly.dependencies = {self.child.label}
                return [poly]
            else:
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
            SIMPLIFY_EVERY = 20
            UNION_THRESHOLD = 5

            if verbose:
                start_time = time.time()

            # Sets of lists of polytopes
            child_sets = [
                r._as_disjoint_union_poly(
                    verbose=verbose, 
                    track_dependencies=track_dependencies
                ) 
                for r in self.child
            ]
            A = child_sets[0]
            i = 0
            for B in child_sets[1:]:
                i += 1

                # # [[Ad-hoc performance optimisation]] that did not improve performance 
                # # check if A is a subset of B, so that A intersect B = A and we may
                # # avoid computing the intersection altogether. 
                # if all(p.is_covered_by(B) for p in A):
                #     continue
                
                new_A = []
                for p in A:
                    inters = []
                    for q in B:
                        # [[Ad-hoc performance optimisation]] - take care of easy cases
                        # first, which tend to occur frequently in practice
                        if p.is_subset_of(q):
                            inters.append(p)
                        elif q.is_subset_of(p):
                            inters.append(q)
                        else:
                            inter = p.intersect(q)
                            if not inter.is_empty(include_boundary=False):
                                if track_dependencies:
                                    inter.dependencies = p.dependencies | q.dependencies
                                inters.append(inter)

                    # [[Ad-hoc performance optimisation]] - frequently in practice p is a 
                    # subset of B, so it is a good candidate for merging
                    if 2 <= len(inters) and len(inters) <= UNION_THRESHOLD:
                        union = Polytope.try_union(inters)
                        if union is not None:
                            if track_dependencies:
                                deps = set()
                                for a in inters: deps = deps | a.dependencies
                                union.dependencies = deps
                            new_A.append(union)
                            continue
                    
                    # Otherwise, just add the pieces as-is. Dependencies (if they are 
                    # being tracked) will propagate automatically
                    new_A.extend(inters)

                
                # [[Ad-hoc performance optimisation]]
                # Every few rounds, simplify if the number of polytopes is too large
                if ((i % SIMPLIFY_EVERY == 0 and len(new_A) >= 100) or len(new_A) >= 500):
                    prevlen = len(new_A)
                    new_A = Region_Helper.simplify_union_of_polys(new_A, Polytope.try_union, track_dependencies)
                    if verbose:
                        print("Simplifying", prevlen, "->", len(new_A))

                if verbose:
                    print("Processed", i, "of", len(child_sets) - 1, ":", len(new_A), "pieces, computed in", time.time() - start_time, "sec")
                A = new_A
            return A
        
        if self.region_type == Region_Type.UNION:
            raise NotImplementedError(self.region_type) # TODO: implement this
        if self.region_type == Region_Type.DISJOINT_UNION:
            result = []
            for r in self.child:
                result.extend(r._as_disjoint_union_poly(
                    verbose=verbose, 
                    track_dependencies=track_dependencies
                ))
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
    
