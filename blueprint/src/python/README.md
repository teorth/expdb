# Python documentation

This subdirectory contains Python functions referenced in the [web blueprint](https://teorth.github.io/expdb/blueprint/). The basic object in this library is a [`Hypothesis`](https://github.com/teorth/expdb/blob/main/blueprint/src/python/hypotheses.py) object, representing either a known theorem or a conjecture. A `Hypothesis` may depend on one or more other `Hypothesis`, and a proof is represented as the dependency tree of `Hypothesis`, starting with the axioms. Important theorems from the literature may be either represented as individual `Hypothesis` (recorded in [`literature.py`](https://github.com/teorth/expdb/blob/main/blueprint/src/python/literature.py)) or as functions that map one or more `Hypothesis` to a `Hypothesis`. 

# Quickstart guide

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
