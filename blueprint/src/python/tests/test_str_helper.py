import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_TESTS_DIR, os.path.dirname(_TESTS_DIR)]

from helpers.str_helper import Str_Helper


def test_zero_expression_is_printed():
    assert Str_Helper.format([0, 0], ["x"]) == "0"
    assert Str_Helper.format([0], ["x"]) == "0"


def test_nonzero_expressions_keep_their_signs():
    assert Str_Helper.format([0, 1], ["x"]) == "x"
    assert Str_Helper.format([1, -2], ["x"]) == "1 - 2x"
    assert Str_Helper.format([0, -1], ["x"]) == "-x"


test_zero_expression_is_printed()
test_nonzero_expressions_keep_their_signs()
