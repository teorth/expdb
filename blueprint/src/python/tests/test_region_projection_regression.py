from polytope import Polytope
from region import Region, Region_Type as T

def box(x, y):
    return Region(T.POLYTOPE, Polytope.rect(x, y))

def run_regressions():
    a = box((0, 1), (0, 1))
    b = box((0, 1), (2, 3))
    assert Region(T.INTERSECT, [a, b]).project({0}) is None
    face = Region(T.INTERSECT, [a, box((0, 1), (1, 2))]).project({0})
    assert face.contains([0]) and face.contains([1])
    assert not face.contains([2])
    union = Region(T.UNION, [a, b])
    projected = Region(T.INTERSECT, [union, box((0, 2), (2, 4))]).project({0})
    assert projected.contains([0]) and not projected.contains([2])
    empty = Region(T.UNION, [])
    assert Region(T.INTERSECT, [a, empty]).project({0}) is None
    assert Region(T.INTERSECT, [Region(T.INTERSECT, [a, b]), a]).project({0}) is None

run_regressions()
