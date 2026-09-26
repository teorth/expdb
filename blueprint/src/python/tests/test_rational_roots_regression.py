from functions import RationalFunction as R, Interval

def run_regressions():
    interval = Interval(-3, 3, True, True)
    assert R([1, -1], [1, 2]).roots(interval) == [1]
    assert R([1, -1], [1, -1]).roots(interval) == []
    assert R([7]).roots(interval) == []
    assert R([1, -1], [1, 2]).roots(Interval(1, 2, False, True)) == []
    # Cross multiplication vanishes at the shared pole, where neither exists.
    assert R([1], [1, 0]).intersections(R([2], [1, 0]), interval) == []
    assert R([1], [1, 0]).intersections(R([1]), interval) == [1]
    assert R([1]).intersections(R([2]), interval) == []

run_regressions()
