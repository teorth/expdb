from polytope import Polytope

def run_regressions():
    positive = Polytope([[0, 1]])
    negative = Polytope([[0, -1]])
    assert positive != negative
    assert positive == Polytope([[0, 7]])
    assert Polytope([[1, 0]]) != positive
    assert Polytope([[1, 0]]) == Polytope([[2, 0]])
    assert Polytope([[0, 1, 0]]) != Polytope([[0, 0, 1]])
    assert Polytope([[0, 1, 0]]) == Polytope([[0, 3, 0]])
    assert Polytope.rect((0, 1)) == Polytope([[0, 2], [3, -3]])
    assert Polytope.rect((0, 0)) != Polytope.rect((0, 0), (0, 0))
    assert positive != object()
    empty = Polytope([[-1, 1], [0, -1]])
    assert empty == Polytope([[-2, 1], [0, -1]])
    assert empty != positive

run_regressions()
