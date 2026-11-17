from polytope import Polytope

def run_regressions():
    box = Polytope.rect((0, 1), (0, 1))
    assert not Polytope.rect((2, 2), (0, 0)).is_covered_by([box])
    line = Polytope.rect((0, 2), (0, 0))
    assert not line.is_covered_by([box])
    assert line.is_covered_by([box, Polytope.rect((1, 2), (-1, 0))])
    assert not line.is_covered_by([Polytope.rect((0, 0), (0, 0)),
                                   Polytope.rect((2, 2), (0, 0))])
    assert Polytope.rect((1, 1), (1, 1)).is_covered_by([box])
    assert not box.is_covered_by([line])
    assert box.is_covered_by([box])
    assert not box.is_covered_by([])
    assert Polytope([[-1, 1], [0, -1]]).is_covered_by([])
    assert not Polytope([[0, 1]]).is_covered_by([Polytope.rect((0, 1))])

run_regressions()
