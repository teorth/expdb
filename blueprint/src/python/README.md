# Python documentation

This subdirectory contains Python functions referenced in the [web blueprint](https://teorth.github.io/expdb/blueprint/). The basic object in this library is a [`Hypothesis`](https://github.com/teorth/expdb/blob/main/blueprint/src/python/hypotheses.py) object, representing either a known theorem or a conjecture. A `Hypothesis` may depend on one or more other `Hypothesis`, and a proof is represented as the dependency tree of `Hypothesis`, starting with the axioms. Important theorems from the literature may be either represented as individual `Hypothesis` (recorded in [`literature.py`](https://github.com/teorth/expdb/blob/main/blueprint/src/python/literature.py)) or as functions that map one or more `Hypothesis` to a `Hypothesis`.

# Installation

The third-party dependencies are listed in [`requirements.txt`](https://github.com/teorth/expdb/blob/main/blueprint/src/python/requirements.txt), and can be installed with
```
pip install -r requirements.txt
```
Note that `pycddlib` must be a 2.x release: the polytope code is written against the 2.x API, which `pycddlib` 3.x renamed. Note also that `pycddlib` 2.x ships no prebuilt wheel for recent Python versions, so installing it may require the `cddlib` and GMP development headers (for instance `libcdd-dev` and `libgmp-dev` on Debian or Ubuntu).

# Quickstart guide

The modules in this directory import each other by plain (unqualified) module name, so all of the commands below should be run from within this directory.

There are several entry points to the project. Some immediate examples can be found by running
```
python examples.py
```
For more examples of derived proofs, including proofs of new results not previously found in the literature, one can run
```
python derived.py
```
To run all unit tests, run
```
python tests/test_all.py
```
An individual test module can be run on its own with, for instance,
```
python -m tests.test_polytope
```
The tests are plain `assert` statements that are executed when the test module is imported, rather than a `pytest` suite.
