# Make both this directory and its parent (the directory holding the modules under test)
# importable, so that this file works when invoked as `python tests/test_all.py` and as
# `python -m tests.test_all`.
import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_TESTS_DIR, os.path.dirname(_TESTS_DIR)]

# Tests will run as part of these imports
import test_affine2

print("All Affine2 test cases passed.")

import test_exppair

print("All exponent pair test cases passed.")

import test_interval

print("All interval test cases passed.")

import test_polynomial

print("All Polynomial test cases passed.")

import test_polytope

print("All Polytope test cases passed.")

import test_rationalfunction

print("All RationalFunction test cases passed.")

import test_region

print("All Region test cases passed.")

import test_ep_to_zd

print("All ep_to_zd regression test cases passed.")

print("All test cases passed.")
