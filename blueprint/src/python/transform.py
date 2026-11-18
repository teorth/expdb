# This class represents a transformation from one Hypothesis to another 
#
# Example use cases: 
# van der Corput A/B transform: 'Exponent pair' -> 'Exponent pair'
# zeta LV theorem transform: 'Large value estimate' -> 'Zeta large value estimate'
# 
# This is only used to simple transformations that takes 1 instance of the 
# input Hypothesis object. For instance, computing bounds on beta via computation
# of a convex hull combines many Hypothesis, and thus will not be implemented 
# as a Transform (instead, it is implemented as a method in the bound_beta
# file)

class Transform:
    
    # Parameters:
    # name: (string) the name of this transformation (e.g. van der Corput A transform)
    # transform: (func Hypothesis -> list of Hypothesis) the transformation function
    def __init__(self, name, transform):
        if not isinstance(name, str):
            raise ValueError('Parameter name must be of type string')
        self.name = name
        self.transform = transform
        
    def __repr__(self):
        return self.name
        