# Contains some helper functions for the Region class
import itertools

class Region_Helper:

    # Simplify a union of polytopes 
    # polys: list of Polytope
    # union_fn: function taking a list of Polytope and returns Polytope if union is 
    # a convex polytope, None otherwise. 
    def simplify_union_of_polys(polys, union_fn):
        
        # The list of polytopes will be expanded, without removal
        # Instead we keep track of the indexes of the polytopes included in the union
        # in this set 
        index = set(i for i in range(len(polys)))

        # Keep track of the tuples (int, int) of indexes which do not form a pairwise 
        # union
        non_unions = set()
        changed = True
        while changed:

            changed = False
            for c in itertools.combinations(list(index), 2):
                (i, j) = c
                # First check if i, j have not been removed from a previous iteration
                if i not in index or j not in index:
                    continue

                # Skip checking the expensive try_union operation if the combination
                # is known to be non-unionable
                if (i, j) in non_unions or (j, i) in non_unions: 
                    continue

                union = union_fn([polys[i], polys[j]])
                if union is not None:
                    # Union is successful: remove indexes in c and add new index
                    index.remove(i)
                    index.remove(j)
                    index.add(len(polys))
                    polys.append(union)
                    changed = True
                else:
                    non_unions.add((i, j))

        return [polys[i] for i in index]
    