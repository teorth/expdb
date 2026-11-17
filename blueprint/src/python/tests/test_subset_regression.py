from polytope import Polytope

def run_regressions():
    positive = Polytope([[0, 1]])
    negative = Polytope([[0, -1]])
    whole = Polytope([[1, 0]])
    assert not positive.is_subset_of(negative)
    assert not whole.is_subset_of(positive)
    assert positive.is_subset_of(whole)
    assert Polytope([[0, 2]]).is_subset_of(positive)
    assert not positive.is_subset_of(Polytope([[0, 1]], linear=True))
    assert Polytope.rect((0, 0)).is_subset_of(positive)
    empty = Polytope([[-1, 1], [0, -1]])
    assert empty.is_subset_of(positive)
    assert empty.is_subset_of(empty)
    assert not positive.is_subset_of(empty)
    assert Polytope.rect((0, 1)).is_subset_of(Polytope.rect((-1, 2)))

run_regressions()
