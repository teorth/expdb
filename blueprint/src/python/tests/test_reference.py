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


def test_trivial_reference_is_not_literature():
    assert Reference.trivial().is_literature() is False
    assert Reference.classical().is_literature() is False
    assert Reference.conjectured().is_literature() is False
    assert Reference.make("Tao", 1999).is_literature() is True


def test_max_year_of_a_conjecture_is_not_a_crash():
    paper = Reference.make("A", 1990)
    assert Reference.max_year((paper, Reference.conjectured())) == "Unknown date"


def test_bare_and_quoted_bibtex_fields_are_kept():
    import tempfile
    bib = """
@article{bareyear,
  author = {Doe, A.},
  title = {Bare year},
  year = 2016
}
@article{quoteyear,
  author = "Doe, A.",
  title = "{Quoted year}",
  year = "2019"
}
"""
    with tempfile.NamedTemporaryFile("w", suffix=".bib", delete=False) as fh:
        fh.write(bib)
        path = fh.name
    mgr = Reference_Manager(path)
    mgr.load()
    assert mgr.get("bareyear").year() == 2016
    assert mgr.get("quoteyear").year() == 2019
    assert "Quoted year" in mgr.get("quoteyear").entries["title"]


def test_repository_bib_keeps_the_dodos_year():
    path = os.path.join(
        os.path.dirname(__file__), "..", "..", "references.bib"
    )
    mgr = Reference_Manager(os.path.abspath(path))
    mgr.load()
    assert mgr.get("dodos").year() == 2016
    assert mgr.get("MR-II").year() == 2020
    assert mgr.get("atkinson_divisor_1941").year() == 1941


test_find_all_skips_entries_without_an_author()
test_find_all_matches_year_via_year_method()
test_max_year_empty_is_minus_one()
test_max_year_unknown_date_is_contagious()
test_trivial_reference_is_not_literature()
test_max_year_of_a_conjecture_is_not_a_crash()
test_bare_and_quoted_bibtex_fields_are_kept()
test_repository_bib_keeps_the_dodos_year()
