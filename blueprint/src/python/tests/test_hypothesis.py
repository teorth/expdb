import os
import sys
from copy import copy

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_TESTS_DIR, os.path.dirname(_TESTS_DIR)]

from hypotheses import Hypothesis, Hypothesis_Set
from reference import Reference


def _hyp(name="leaf"):
    return Hypothesis(name, "Upper bound on mu", type("D", (), {"__repr__": lambda self: "d"})(), "leaf", Reference.classical())


def test_proof_depth_on_a_leaf():
    # max() of an empty generator used to raise here.
    h = _hyp()
    assert h.proof_depth() == 1
    assert h.proof_complexity() == 1
    assert h.proof_date() == -1


def test_proof_depth_on_a_tree():
    leaf = _hyp("a")
    mid = Hypothesis("mid", "Upper bound on mu", type("D", (), {"__repr__": lambda self: "d"})(), "mid", Reference.derived(1990))
    mid.add_dependency(leaf)
    root = Hypothesis("root", "Upper bound on mu", type("D", (), {"__repr__": lambda self: "d"})(), "root", Reference.derived(2001))
    root.add_dependency(mid)
    assert mid.proof_depth() == 2
    assert root.proof_depth() == 3
    assert root.proof_complexity() == 3
    assert root.proof_date() == 2001


def test_hypothesis_set_copy_does_not_share_cached_data():
    hs = Hypothesis_Set()
    hs.data["hull"] = [1]
    other = copy(hs)
    other.data["hull"] = [2]
    assert hs.data["hull"] == [1]


def test_add_hypotheses_accepts_a_tuple():
    a = _hyp("a")
    b = _hyp("b")
    hs = Hypothesis_Set()
    hs.add_hypotheses((a, b))
    assert len(hs) == 2


def test_find_hypothesis_returns_none_quietly():
    hs = Hypothesis_Set()
    assert hs.find_hypothesis(name="missing") is None


def test_is_match_accepts_year_as_string():
    h = Hypothesis("x", "Upper bound on mu", type("D", (), {"__repr__": lambda self: "d"})(), "p", Reference.make("X", 1999))
    assert h.is_match(year="2000")
    assert not h.is_match(year="1990")


def test_hypothesis_set_iterates():
    a = _hyp("a")
    hs = Hypothesis_Set([a])
    assert list(hs) == [a]


test_proof_depth_on_a_leaf()
test_proof_depth_on_a_tree()
test_hypothesis_set_copy_does_not_share_cached_data()
test_add_hypotheses_accepts_a_tuple()
test_find_hypothesis_returns_none_quietly()
test_is_match_accepts_year_as_string()
test_hypothesis_set_iterates()
