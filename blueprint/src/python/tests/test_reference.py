import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_TESTS_DIR, os.path.dirname(_TESTS_DIR)]

from reference import Reference, Reference_Manager


def test_find_all_skips_entries_without_an_author():
    mgr = Reference_Manager("/dev/null")
    mgr.refs["a"] = Reference("a", "article", {"year": 1999})
    mgr.refs["b"] = Reference("b", "article", {"author": "Tao", "year": 1999})
    assert mgr.find_all(author="Tao") == [mgr.refs["b"]]


def test_find_all_matches_year_via_year_method():
    mgr = Reference_Manager("/dev/null")
    mgr.refs["a"] = Reference("a", "article", {"author": "Tao", "year": 1999})
    assert mgr.find_all(year=1999) == [mgr.refs["a"]]
    assert mgr.find_all(year="1999") == [mgr.refs["a"]]
    assert mgr.find_all(year=2000) == []


def test_max_year_coerces_numeric_strings():
    r1 = Reference.make("A", "1990")
    r2 = Reference.make("B", 2001)
    assert Reference.max_year((r1, r2)) == 2001


def test_max_year_empty_is_minus_one():
    assert Reference.max_year(()) == -1


def test_max_year_unknown_date_is_contagious():
    known = Reference.make("A", 1990)
    unknown = Reference("u", "article", {})
    assert unknown.year() == "Unknown date"
    assert Reference.max_year((known, unknown)) == "Unknown date"


test_find_all_skips_entries_without_an_author()
test_find_all_matches_year_via_year_method()
test_max_year_empty_is_minus_one()
test_max_year_unknown_date_is_contagious()
