from polytope import Polytope
from region import Region, Region_Type
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

def test_as_disjoint_union():

    # A union of overlapping rectangles, rewritten as a disjoint union: the
    # pieces must cover exactly the same points, and no point may lie in the
    # interior of two of them.
    for i in range(50):
        rects = []
        for _ in range(3):
            xlims = (rd.uniform(0, 3), rd.uniform(0, 3))
            ylims = (rd.uniform(0, 3), rd.uniform(0, 3))
            rects.append(
                Polytope.rect((min(xlims), max(xlims)), (min(ylims), max(ylims)))
            )

        union = Region.union([Region.from_polytope(r) for r in rects])
        pieces = union.as_disjoint_union()

        for _ in range(20):
            x = (rd.uniform(0, 3), rd.uniform(0, 3))
            assert pieces.contains(x) == any(r.contains(x) for r in rects)

        # Pairwise disjoint: any two pieces may share a face, not an interior.
        polys = [piece.child for piece in pieces.child]
        for j in range(len(polys)):
            for k in range(j + 1, len(polys)):
                overlap = polys[j].intersect(polys[k])
                assert overlap.is_empty(include_boundary=False)

test_union()
test_intersect()
test_as_disjoint_union()


def test_empty_boolean_combinations_raise():
    try:
        Region.union([])
        assert False
    except ValueError:
        pass
    try:
        Region.intersect(())
        assert False
    except ValueError:
        pass


def test_region_accepts_a_tuple_of_children():
    r = Region.from_polytope(Polytope.rect((0, 1), (0, 1)))
    u = Region.union((r, r))
    assert u.contains((0.5, 0.5))


test_empty_boolean_combinations_raise()
test_region_accepts_a_tuple_of_children()


def test_lift_preserves_complement():
    box = Polytope.rect((0, 1), (0, 1))
    R = Region.complement(Region.from_polytope(box))
    lifted = R.lift([0, 1, (0, 1)])
    assert lifted.region_type == Region_Type.COMPLEMENT
    assert isinstance(lifted.child, Region)
    # inside the lifted box => outside the complement
    assert not lifted.contains((0.5, 0.5, 0.5))
    assert lifted.contains((2, 2, 0.5))


test_lift_preserves_complement()
