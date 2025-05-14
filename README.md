



# Analytic Number Theory Exponent Database

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-lightblue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Code Style: Black](https://img.shields.io/badge/Code%20Style-Black-000000.svg)](https://github.com/psf/black)

This is the repository for the Analytic Number Theory Exponent Database (ANTEDB), an ongoing project to systematically record known theorems for various exponents appearing in analytic number theory, as well as the relationships between them. Currently, the database is recording facts and conjectures on [exponent pairs](https://teorth.github.io/expdb/blueprint/exponent-pairs-chapter.html), [exponential sum bounds](https://teorth.github.io/expdb/blueprint/beta-chapter.html), [zero-density](https://teorth.github.io/expdb/blueprint/zero-density-chapter.html) and [moment bounds](https://teorth.github.io/expdb/blueprint/zeta-moment-chapter.html) for the Riemann zeta-function, [large value estimates](https://teorth.github.io/expdb/blueprint/largevalue-chapter.html), [additive energy bounds](https://teorth.github.io/expdb/blueprint/energy-chapter.html) and exponents related to [prime distributions](https://teorth.github.io/expdb/blueprint/primes-sec.html), amongst other results (the full list of exponents can be found [here](https://teorth.github.io/expdb/blueprint/intro-chapter.html)). The database aims to organise and record theorems and relationships as both [human-readable proofs](https://teorth.github.io/expdb/blueprint/intro-chapter.html) and as [executable python code](https://github.com/teorth/expdb/tree/main/blueprint/src/python). This builds upon the tables of exponent pairs and other related exponents that were collected in [this paper of Trudgian and Yang](https://arxiv.org/abs/2306.05599), and aims to fulfil the vision set out in [this blog post of Terence Tao](https://terrytao.wordpress.com/2024/07/07/a-computation-outsourced-discussion-of-zero-density-theorems-for-the-riemann-zeta-function/).

This is intended to be a living database; additional corrections, updates, and contributions to the ANTEDB are welcome.  Instructions for contributors can be [found here](https://github.com/teorth/expdb/blob/main/CONTRIBUTING.md).

- [Main web page](https://teorth.github.io/expdb/) for this project.
- [Web blueprint](https://teorth.github.io/expdb/blueprint/), containing a human-readable version of the database that is accessible online (see also the [introduction](https://teorth.github.io/expdb/blueprint/intro-chapter.html) here).
- [Dependency graph](https://teorth.github.io/expdb/blueprint/dep_graph_document.html) of the theorems in the database.
- [PDF form of blueprint](https://teorth.github.io/expdb/blueprint.pdf), a downloadable, human-readable version of the database stored as a single file.
- [Python code repository](https://github.com/teorth/expdb/tree/main/blueprint/src/python), containing a programmable version of this database and various optimization routines.

## External links
- T. Tao, T. Trudgian, A. Yang, "[New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach](https://arxiv.org/abs/2501.16779)": A paper describing the project, and several of the new results that already arose from it
- [A blog post](https://terrytao.wordpress.com/2025/01/28/new-exponent-pairs-zero-density-estimates-and-zero-additive-energy-estimates-a-systematic-approach/) by Terence Tao describing this paper and project
- [Mathbases](https://github.com/MathBases/MathBases) - A directory of math databases (including this one)

## Python database and proof automation
Each theorem and conjecture in the database is represented both in [natural language](https://teorth.github.io/expdb/blueprint/intro-chapter.html) and also in Python using the `Hypothesis` object. Each `Hypothesis` may either contain a reference to a result from the literature or depend on one or more other `Hypothesis` objects. The special predefined set `literature` contains all hypotheses of the first kind. For example, Ingham's 1940 [zero-density estimate](https://teorth.github.io/expdb/blueprint/zero-density-chapter.html) for the Riemann zeta function 
```math
N(\sigma, T) \ll T^{3(1-\sigma)/(2-\sigma) + o(1)}\qquad (1/2 \le \sigma \le 1, T \to \infty)
```
is represented as a `Hypothesis` in `literature`, and can be retrieved using
```
import literature

h = literature.find_hypothesis(hypothesis_type="Zero density estimate", keywords="Ingham")
print("Hypothesis name:", h)
print("Hypothesis data:", h.data)
print("Hypothesis reference:", h.reference.title(), h.reference.author(), h.reference.year())
print("Hypothesis dependencies:", h.dependencies)
``` 
Console output
```
Hypothesis name: Ingham (1940) zero density estimate
Hypothesis data: A(x) \leq 3/(2 - x) on [1/2,1)
Hypothesis reference: {On} the estimation of ${N}(\sigma, T)$ Ingham 1940
Hypothesis dependencies: set()
```
There are no dependencies because this `Hypothesis` object directly references the literature. Alternatively, theorems may also be represented as a `Hypothesis` that (recursively) depend on other hypotheses. The dependency structure can be represented as a tree whose root is the theorem and whose leaves are other known theorems (either trivial or proved in the literature). This tree also represents a proof of the `Hypothesis`. For instance, we may also derive Ingham's estimate using other hypotheses in the database, such as [large value estimates](https://teorth.github.io/expdb/blueprint/largevalue-chapter.html). This returns a `Hypothesis` object representing the same result but containing a proof of Ingham's result. represented as a dependency tree. 
```
import derived

h = prove_zero_density_ingham_1940_v2()
h.recursively_list_proofs()
``` 
Console output
```
- [Derived zero density estimate]  i.e. A(x) \leq \frac{3 \left(x - 1\right)}{x - 2} on [1/2,1). Follows from computed large value estimates and zeta large value estimates with τ0 = 2 - σ for σ \in [1/2,1)). Dependencies:
        - [Derived large value estimate]  i.e. (σ,τ,ρ) in Disjoint union of {['-4/3 + 2/3σ + τ >= 0', '2 - σ - τ >= 0', '3/4 - σ >= 0', 'ρ >= 0', '-1/2 + σ >= 0', '2 - 2σ - ρ >= 0'], ['2 - σ - τ >= 0', '3/4 - σ >= 0', '-1/2 + σ >= 0', '1 - 2σ + τ - ρ >= 0', '-2 + 2σ + ρ >= 0']}. Follows from 1 large value estimates. Dependencies:
                - [Classical large value estimate]  i.e. ρ <= max(2 - 2σ, 1 - 2σ + τ). Classical.
        - [Derived zeta large value estimate]  i.e. (σ,τ,ρ) in Disjoint union of {['8/3 - 4/3σ + τ >= 0', '-2 + τ >= 0']}. Follows from 1 zeta large value estimates. Dependencies:
                - [Classical large value estimate]  i.e. ρ <= max(2 - 2σ, 1 - 2σ + τ). Classical.
```
The database also implements a limited form of automated theorem proving (ATP) using known relationships between certain exponents in analytic number theory. A number of functions can be used to automatically find proofs given a desired result and a set of assumed `Hypothesis`. For example, to find a proof of the [exponent pair](https://teorth.github.io/expdb/blueprint/exponent-pairs-chapter.html) $(\frac{4}{18}, \frac{11}{18}) = BABA^2B(0, 1)$, one can use
```
from exponent_pair import *

h = best_proof_of_exponent_pair(frac(4,18), frac(11,18))
h.recursively_list_proofs()
```
Console output
```
- [Derived exponent pair (2/9, 11/18)]  i.e. The exponent pair (2/9, 11/18). Follows from "Derived exponent pair (1/9, 13/18)" and taking the van der Corput B transform. Dependencies:
        - [van der Corput B transform]  i.e. van der Corput B transform. See [van der Corput, 1920]. 
        - [Derived exponent pair (1/9, 13/18)]  i.e. The exponent pair (1/9, 13/18). Follows from "Derived exponent pair (2/7, 4/7)" and taking the van der Corput A transform. Dependencies:
                - [van der Corput A transform]  i.e. van der Corput A transform. See [van der Corput, 1920].
                - [Derived exponent pair (2/7, 4/7)]  i.e. The exponent pair (2/7, 4/7). Follows from "Derived exponent pair (1/14, 11/14)" and taking the van der Corput B transform. Dependencies:
                        - [van der Corput B transform]  i.e. van der Corput B transform. See [van der Corput, 1920].
                        - [Derived exponent pair (1/14, 11/14)]  i.e. The exponent pair (1/14, 11/14). Follows from "Derived exponent pair (1/6, 2/3)" and taking the van der Corput A transform. Dependencies:
                                - [van der Corput A transform]  i.e. van der Corput A transform. See [van der Corput, 1920].
                                - [Derived exponent pair (1/6, 2/3)]  i.e. The exponent pair (1/6, 2/3). Follows from "Derived exponent pair (1/2, 1/2)" and taking the van der Corput A transform. Dependencies:
                                        - [van der Corput A transform]  i.e. van der Corput A transform. See [van der Corput, 1920].
                                        - [Derived exponent pair (1/2, 1/2)] i.e. The exponent pair (1/2, 1/2). Follows from "Trivial exponent pair (0, 1)" and taking the van der Corput B transform. Dependencies:
                                                - [van der Corput B transform]  i.e. van der Corput B transform. See [van der Corput, 1920].
                                                - [Trivial exponent pair (0, 1)]  i.e. The exponent pair (0, 1). Triangle inequality.
```