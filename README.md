
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
