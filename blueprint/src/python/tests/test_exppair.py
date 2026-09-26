from functions import *
from hypotheses import *
from literature import *
from exponent_pair import *
from reference import *


def run_exp_pair_transform_tests():

    transforms = literature.list_hypotheses(hypothesis_type="Exponent pair transform")
    A_process = next(t for t in transforms if t.name == "van der Corput A transform")
    B_process = next(t for t in transforms if t.name == "van der Corput B transform")
    E = literature_exp_pair(frac(1, 6), frac(2, 3), Reference.make("Test", 2024))

    # B(1/6, 2/3) = (1/6, 2/3)
    BE = B_process.data.transform(E)
    assert BE.data.k == frac(1, 6) and BE.data.l == frac(2, 3)

    # A(1/6, 2/3) = (1/14, 11/14)
    AE = A_process.data.transform(E)
    assert AE.data.k == frac(1, 14) and AE.data.l == frac(11, 14)

run_exp_pair_transform_tests()


def test_prove_convex_combination_of_two_pairs():
    left = trivial_exp_pair
    right = literature_exp_pair(frac(1, 2), frac(1, 2), Reference.make("Test", 2024))
    hypotheses = Hypothesis_Set([left, right])
    for deep_search in (False, True):
        result = construct_proof(frac(1, 4), frac(3, 4), hypotheses,
                                     reduce_dependencies=True, deep_search=deep_search)
        assert result.data == Exp_pair(frac(1, 4), frac(3, 4))
        assert result.dependencies == {left, right}
    assert construct_proof(frac(1, 4), frac(1, 2), hypotheses) is None
    assert construct_proof(0, 1, Hypothesis_Set([left])) is left


test_prove_convex_combination_of_two_pairs()


def test_triangle_containment_is_independent_of_vertex_order():
    import itertools
    for vertices in itertools.permutations([(0, 0), (1, 0), (0, 1)]):
        assert in_triangle(*vertices, (frac(1, 4), frac(1, 4)))
        assert in_triangle(*vertices, (frac(1, 2), frac(1, 2)))
        assert not in_triangle(*vertices, (1, 1))
    assert in_triangle((0, 0), (1, 0), (2, 0), (frac(1, 2), 0))
    assert not in_triangle((0, 0), (1, 0), (2, 0), (3, 0))
    assert not in_triangle((0, 0), (0, 0), (0, 0), (1, 1))


test_triangle_containment_is_independent_of_vertex_order()
