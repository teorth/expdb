
# Basic code for handling hypotheses in the exponent database.  The code here is
# intended to be as flexible and general as possible, anticipating that different
# types of hypotheses will be manipulated in different ways.

########################################################################################
# A Hypothesis is a mathematical assertion that was in the literature, conjectured,
# or derived from other hypotheses.
# It consists of some data (under the 'data' attribute) together with metadata.
# More precisely, it contains the following information:
#
# - name : the short name of the hypothesis, e.g., 'Bourgain bound on \mu(1/2)'
# - hypothesis_type : the type of the hypothesis.  Current supported types:
#    - 'Upper bound on beta'
#    - 'Exponent pair'
#    - 'Exponent pair transform'
#    - 'Upper bound on mu'
#    - 'Large value estimate'
#    - 'Zeta large value estimate'
#    - 'Zero density estimate'
#    - 'Zero density energy estimate'
#    - 'Large value energy region'
#    - 'Large value energy region transform'
# - description : a longer description of the hypothesis, e.g., 'The bound
#       $\mu(1/2) \leq 13/84$'.  Defaults to the description of the underlying data.
# - citation : a reference to where the hypothesis was discovered/proven/conjectured,
#       e.g., 'A paper of Bourgain'.  Use 'Derived' if the hypothesis was generated
#       by the exponent database as a consequence of other existing hypotheses.
# - proof : a human-readable proof of the hypothesis.  (If from the literature,
#       the proof will simply cite the literature.  If a conjecture, states the
#       hypothesis as a conjecture).
# - data : the function or class that computes whatever the hypothesis is., e.g.
#       Bound_mu(1/2, 13/84).  The type of `data` depends on `hypothesis_type`
# - year : the year in which the result was established.  Use 'Derived' if the
#       hypothesis is derived within the database, 'Conjectured' if it is a
#       conjecture rather than a proven result, or 'Classical' if it is too standard
#       to assign a year to.
# - dependencies : the set of hypotheses that the current hypothesis directly
#       depends on (defaults to the empty set).

# TODO: add `Lean proof` attribute


class Hypothesis:
    def __init__(self, name, hypothesis_type, data, proof, reference):
        self.name = name
        self.hypothesis_type = hypothesis_type
        self.description = data.__repr__()
        self.proof = proof
        self.data = data
        self.reference = reference
        self.dependencies = set()

    # By default, just returns the name of the hypothesis.
    def __repr__(self):
        return self.name

    # Returns the name and description.
    def desc(self):
        return f"{self.name} ({self.description})"

    # Returns the name, description, human-readable proof of the hypothesis.
    def desc_with_proof(self):
        return f"{self.name}: {self.description}.  {self.proof}."

    # Adds a single hypothesis to the set of dependencies.
    def add_dependency(self, fact):
        self.dependencies.add(fact)

    # internal method to add all dependencies of a hypothesis to a Hypothesis_Set
    def add_recursive_dependencies(self, dependency_set):
        dependency_set.add_hypotheses(self.dependencies)
        for subhypothesis in list(self.dependencies):
            subhypothesis.add_recursive_dependencies(dependency_set)

    # Returns the set of all hypotheses that the current hypothesis depends on, recursively.
    def recursive_dependencies(self):
        dependency_set = Hypothesis_Set()
        self.add_recursive_dependencies(dependency_set)
        return dependency_set

    def recursively_list_proofs(self, indentation=0):
        if len(self.dependencies) > 0:
            print(
                "\t" * indentation
                + f"- [{self}]  i.e. {self.data}. {self.proof}. Dependencies:"
            )
            for d in self.dependencies:
                d.recursively_list_proofs(indentation + 1)
        else:
            print("\t" * indentation + f"- [{self}]  i.e. {self.data}. {self.proof}.")

    # Returns True if the hypothesis is of the given type and was established by the given year.
    def is_match(self, hypothesis_type="Any", year="Any"):
        if hypothesis_type == "Any" or hypothesis_type == self.hypothesis_type:
            if (
                year == "Any"
                or self.reference.label == "Derived"
                or self.reference.label == "Classical"
            ):
                return True
            if self.reference.label == "Conjectured":
                return False
            if (
                self.reference.year() != "Unknown date"
                and self.reference.year() <= year
            ):
                return True
        return False

    # Returns the complexity of the proof of this hypothesis
    def proof_complexity(self):
        return sum(h.proof_complexity() for h in self.dependencies) + 1

    # Returns the date of the latest dependency
    def proof_date(self):
        year = self.reference.year()
        if year == "Unknown date":
            year = -1
        return max([year] + [h.proof_date() for h in self.dependencies])


########################################################################################
# A Hypothesis_Set is, as the name suggests, a set of hypotheses.  This will be an input parameter in various exponent calculation routines in the database, which seek to answer questions like "What is the best zero-density bound known for sigma = 3/4 assuming [a given list of hypotheses]"

# In the future, we can allow lists of hypotheses to store precomputed data to speed up calculations.  For instance, if the list of hypotheses contains exponent pairs, we can permit the list to store the convex hull of these pairs.  A routine that wants to compute this convex hull can first look to see if the convex hull has already been computed, and use it if available, otherwise it will compute it and store it in the list.  This way, such computations only need to be performed once.

# Similarly, we can develop methods to self-improve lists of hypotheses.  For instance, if a list of hypotheses contains exponent pairs, we can create methods to discover all zeta function moment bounds that can be derived from those pairs, add those to the list of hypotheses, and remove duplicate or redundant hypotheses in the combined list.


class Hypothesis_Set:
    def __init__(self, hypotheses=set()):
        self.hypotheses = set()
        self.add_hypotheses(hypotheses)
        self.data = {}
        self.data_valid = False  # set to false whenever data needs to be recomputed

    def __repr__(self):
        return f'Set of {len(self.hypotheses)} hypotheses: [{",".join(h.name for h in self.hypotheses)}]'

    # Shallow copy, the hypothesis objects are not cloned
    def __copy__(self):
        copy = Hypothesis_Set(self.hypotheses)
        copy.data = self.data
        copy.data_valid = self.data_valid
        return copy

    def __iter__(self):
        return self.hypotheses.__iter__()

    def __next__(self):
        return self.hypotheses.__next__()

    def __len__(self):
        return len(self.hypotheses)

    def to_list(self):
        return list(self.hypotheses)

    def list_proofs(self):
        for hypothesis in self:
            print(hypothesis.desc_with_proof())

    # Adds a single hypothesis to the hypothesis set.
    # If invalidate_data is true, then the precomputed, cached data (e.g. a convex hull)
    # is invalidated and will be recomputed at the next proof function call. In certain
    # cases calls of this function are guaranteed to not require recomputation of cached
    # data (e.g. in prove_mu_bound where we insert a valid bound on \mu(\sigma) derived
    # from other bounds; such a bound is guaranteed to not enlarge the precomputed convex
    # hull)
    def add_hypothesis(self, hypothesis, invalidate_data=True):
        self.hypotheses.add(hypothesis)
        if invalidate_data:
            self.data_valid = False

    # Adds a set or list of hypotheses, an individual hypothesis, or a Hypothesis_Set
    def add_hypotheses(self, new_hypotheses, invalidate_data=True):
        if isinstance(new_hypotheses, Hypothesis):
            self.add_hypothesis(new_hypotheses, invalidate_data)
        elif isinstance(new_hypotheses, Hypothesis_Set):
            self.add_hypotheses(new_hypotheses.hypotheses, invalidate_data)
        elif isinstance(new_hypotheses, list):
            self.add_hypotheses(set(new_hypotheses), invalidate_data)
        else:
            self.hypotheses.update(new_hypotheses)
            if invalidate_data:
                self.data_valid = False

    # return all hypotheses of a given type, and (optionally) up to a given year.  Note: is now returning a Hypothesis_Set rather than a list
    def list_hypotheses(self, hypothesis_type="Any", year="Any"):
        return [
            hypothesis
            for hypothesis in self
            if hypothesis.is_match(hypothesis_type, year)
        ]

    
    def find_hypothesis(
        self, hypothesis_type="Any", data="Any", name="Any", keywords="Any", year="Any"
    ) -> Hypothesis | None:
        """
        Returns the first instance of a Hypothesis in the set that matches the 
        specified requirements.

        Parameters
        ----------
        hypothesis_type : str, optional
            The type of the hypothesis, e.g. "Exponent pair" (default is "Any").
        data : str or object, optional
            The data that the hypothesis contains, e.g. an object of type Exp_pair
            (default is "Any").
        name : str, optional
            The full name of the hypothesis, e.g. "Jutila large value theorem with k = 3"
            (default is "Any").
        keywords : str, optional
            A comma-separated list of keywords to search for in the name of the 
            hypothesis (default is "Any").
        year : str or int, optional
            The year of the hypothesis (default is "Any").
        
        Returns
        -------
        Hypothesis or None
            The first hypothesis that matches all conditions, or None if no such 
            hypothesis exists. 
        """
        for h in self:
            if hypothesis_type == "Any" or h.hypothesis_type == hypothesis_type:
                if data == "Any" or h.data == data:
                    if name == "Any" or h.name == name:
                        if keywords == "Any" or all(
                            k.strip() in h.name for k in keywords.split(",")
                        ):
                            if year == "Any" or h.reference.year() == year:
                                return h
        print("ERROR: No matching hypothesis found")
        return None
