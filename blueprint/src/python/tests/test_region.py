from ..polytope import Polytope
from ..region import Region
import random as rd


def test_union():

    # Generate some random regions
    for i in range(100):
        xlims = (rd.uniform(0, 3), rd.uniform(0, 3))
        ylims = (rd.uniform(0, 3), rd.uniform(0, 3))
        R1 = Region.from_polytope(Polytope.rect((min(xlims), max(xlims)), (min(ylims), max(ylims))))

        xlims = (rd.uniform(0, 3), rd.uniform(0, 3))
        ylims = (rd.uniform(0, 3), rd.uniform(0, 3))
        R2 = Region.from_polytope(Polytope.rect((min(xlims), max(xlims)), (min(ylims), max(ylims))))

        union = Region.union([R1, R2])
        x = (rd.uniform(0, 3), rd.uniform(0, 3))
        assert union.contains(x) == (R1.contains(x) or R2.contains(x))

def test_intersect():
    # Generate some random regions
    for i in range(100):
        xlims = (rd.uniform(0, 3), rd.uniform(0, 3))
        ylims = (rd.uniform(0, 3), rd.uniform(0, 3))
        R1 = Region.from_polytope(Polytope.rect((min(xlims), max(xlims)), (min(ylims), max(ylims))))

        xlims = (rd.uniform(0, 3), rd.uniform(0, 3))
        ylims = (rd.uniform(0, 3), rd.uniform(0, 3))
        R2 = Region.from_polytope(Polytope.rect((min(xlims), max(xlims)), (min(ylims), max(ylims))))

        intersection = Region.intersect([R1, R2])
        x = (rd.uniform(0, 3), rd.uniform(0, 3))
        assert intersection.contains(x) == (R1.contains(x) and R2.contains(x))

test_union()
test_intersect()
