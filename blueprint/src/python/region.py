from polytope import Polytope

class Region_Type:
    
    # Supports many regions 
    INTERSECT = 0,
    UNION = 1,
    
    # Supports one region or polytope
    POLYTOPE = 2,
    COMPLEMENT = 3
    
    def to_str(region_type):
        if region_type == Region_Type.INTERSECT:
            return "Intersection"
        elif region_type == Region_Type.UNION:
            return "Union"
        elif region_type == Region_Type.POLYTOPE:
            return "Polytope"
        elif region_type == Region_Type.COMPLEMENT:
            return "Complement"
        raise ValueError("Unknown type", region_type)
    
    
    
# A region represents a boolean combination of polytopes 
class Region:
    
    # Construct a region from a polytope
    def __init__(self, region_type, region):
        if not isinstance(region, Polytope) and \
            not isinstance(region, Region):
            raise ValueError("Parameter polytope must either be of type Polytope or Region")
            
        self.region_type = region_type
        self.child = region

    def __repr__(self):
        return self.to_str(0)

    def to_str(region, indentation=0):
        s = Region_Type.to_str(region.region_type) + r"\n"
        if isinstance(region.child, Region):
            s += ("\t" * indentation) + str(region.child)
        else:
            for r in region.child:
                s += ("\t" * indentation) + r + r"\n"
        return s
        
    # Static methods ---------------------------------------------------------

    # Construct a region from a polytope
    def from_polytope(polytope):
        return Region(Region_Type.POLYTOPE, polytope)

    # Construct a region as the complement of another region 
    def complement(region):
        return Region(Region_Type.COMPLEMENT, region)

    # Compute the union of regions 
    def union(*regions):
        return Region(Region_Type.UNION, list(regions))

    # Compute the intersection of regions
    def intersect(*regions):
        return Region(Region_Type.INTERSECT, list(regions))

    # Instance methods -------------------------------------------------------

    # Returns whether this region contains the point x
    def contains(self, x):
        if self.region_type == Region_Type.COMPLEMENT:
            return not self.child.contains(x)
        if self.region_type == Region_Type.POLYTOPE:
            return self.child.contains(x)
        if self.region_type == Region_Type.UNION:
            return any(c.contains(x) for c in self.child)
        if self.region_type == Region_Type.INTERSECT:
            return all(c.contains(x) for c in self.child)

        