from polytope import Polytope
import copy

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

    # Returns a representation of this region as a disjoint union of convex polytopes.
    # Returns the result as a list of disjoint Polytope objects
    def to_disjoint_union(self):
        if self.region_type == Region_Type.COMPLEMENT:
            raise NotImplementedError() # TODO: implement this
        if self.region_type == Region_Type.POLYTOPE:
            return [self.child]
        if self.region_type == Region_Type.INTERSECT:
            # Sets of lists of polytopes
            child_sets = [r.to_disjoint_union() for r in self.child]
            disjoint = child_sets[0]
            for i in range(1, len(child_sets)):
                new_disjoint = []
                for p in disjoint:
                    for q in child_sets[i]:
                        inter = p.intersect(q)
                        if not inter.is_empty(include_boundary=False):
                            new_disjoint.append(inter)
                print(len(new_disjoint))
                disjoint = new_disjoint
            return disjoint
        if self.region_type == Region_Type.UNION:
            raise NotImplementedError() # TODO: implement this
        if self.region_type == Region_Type.DISJOINT_UNION:
            result = []
            for r in self.child:
                result.extend(r.to_disjoint_union())
            return result

