from polytope import Polytope
from region import Region


def test_simplification_preserves_isolated_point():
    point = Polytope.rect((2, 2), (3, 3))
    region = Region.union([Region.from_polytope(point)])
    assert region.contains((2, 3))
    region.simplify()
    assert region.contains((2, 3))
    assert not region.contains((2, 4))


def test_simplification_preserves_segment_outside_full_dimensional_component():
    rectangle = Polytope.rect((0, 1), (0, 1))
    segment = Polytope.rect((2, 2), (0, 1))
    region = Region.union([Region.from_polytope(rectangle), Region.from_polytope(segment)])
    points = [(0, 0), (0.5, 0.5), (2, 0), (2, 0.5), (2, 1), (1.5, 0.5)]
    before = [region.contains(point) for point in points]
    region.simplify()
    assert [region.contains(point) for point in points] == before
