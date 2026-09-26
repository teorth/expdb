from polytope import Polytope

def run_regressions():
    strip = Polytope([[0, 1, 0], [0, 0, 1], [1, 0, -1]])
    ray = strip.project({0})
    assert ray.contains([2]) and not ray.contains([-1])
    bounded = strip.project({1})
    assert bounded.contains([0]) and bounded.contains([1])
    assert not bounded.contains([2])
    whole_line = Polytope([[0, 1, 0], [1, -1, 0]]).project({1})
    assert whole_line.contains([-100]) and whole_line.contains([100])
    negative = Polytope([[0, -1, 0], [0, 0, 1]], linear=False).project({0})
    assert negative.contains([-2]) and not negative.contains([1])
    assert strip.project(set()).contains([])
    assert Polytope([[-1, 1], [0, -1]]).project({0}) is None
    for dims in ({-1}, {2}, {"x"}):
        try:
            strip.project(dims)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid projection coordinate accepted")

run_regressions()
