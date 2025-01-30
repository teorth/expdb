# Contributing to Analytic Number Theory Exponent Database

Contributions to the ANTEDB are very welcome!  For short contributions (such as correcting a typo, or adding a reference or citation), one can simply email the correction to one of the maintainers of the ANTEDB ([Terence Tao](mailto:tao@math.ucla.edu), [Timothy Trudgian](mailto:timothy.trudgian@unsw.edu.au), or [Andrew Yang](mailto:andrew.yang1@unsw.edu.au)).  Alternatively, if one has familiarity with Github, one can submit a [Github Pull Request](https://github.com/teorth/expdb/pulls) (PR) to the [Github repository](https://github.com/teorth/expdb) with the suggested changes, or raise a [Github issue](https://github.com/teorth/expdb/issues).

Some specific information about types of contributions:

## Contributing a new reference to the blueprint

The references are stored at https://github.com/teorth/expdb/blob/main/blueprint/src/references.bib in BibTeX format.  If possible, please try to maintain the alphabetical ordering of the references by first order.

## Contributing to a chapter of the blueprint

The chapters are stored at https://github.com/teorth/expdb/tree/main/blueprint/src/chapter .  One can use standard LaTeX commands in these chapters.  One can also use the \uses{} macro to indicate which results depend on which other ones, for the purpose of filling out the dependency graph, though at the current stage of the project we are not making heavy use of this feature.

It is fine to suggest incomplete contributions, for instance stating a result with only a very sketchy proof or reference.  For an extremely incomplete contribution (e.g., a vague statement of a result without a reference, or proposing a new direction for the ANTEDB without contributing significant content), consider opening up a Github issue instead of a pull request.

If you propose to add a new chapter, or rearrange existing ones, one will also have to modify https://github.com/teorth/expdb/blob/main/blueprint/src/content.tex .  It is also recommended that each chapter be labeled, in order for the blueprint to be able to assign a stable name to the web page for that chapter.

## Contributing to Python code

Python code is stored at https://github.com/teorth/expdb/tree/main/blueprint/src/python .  References to python code within the blueprint can be made using the (somewhat crude) `\code{}`, `\python{}`, `\literature{}`, and `\derived{}` macros, defined in https://github.com/teorth/expdb/blob/main/blueprint/src/macros/common.tex .

## Contributing to Lean code

At present we do not have any Lean code in this database, although we have structured it to allow for future expansion in this direction.  If you have some ideas for such an expansion, please do not hesitate to discuss them with the maintainers of the ANTEDB, or raise a suitable [issue on the Github repository](https://github.com/teorth/expdb/issues); if there is sufficient interest, we may start a broader discussion to take more steps in this direction.

## A note on notational conventions

The ANTEDB is intended to collect results on analytic number theory exponents in the literature in such a way that they can be combined as readily as possible.  Because of this, we have made efforts to select standardized notations for the ANTEDB, and to express results from the literature in this standardized notation even if they are expressed in a somewhat different notation in the original source material.

The focus on the ANTEDB is thus on the exponent bounds in the literature, rather than the notational choices from that literature, and the addition of notational conventions to the ANTEDB that duplicate ones that are already present in the database is therefore discouraged in order to maximize interoperability of results.  However, it is conceivable that there is a case for globally refactoring the ANTEDB to switch from a current notation convention to new one; if so, one should set up a [Github issue](https://github.com/teorth/expdb/issues) to discuss the pros and cons of such a global switch.
