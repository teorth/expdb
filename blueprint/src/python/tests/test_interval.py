import parent

from functions import *

def run_all():
    i1 = Interval(0, 2, True, False)  # [0, 2)
    i2 = Interval(0, 2)
    assert i1 == i2
    # [0, 2) intersect (1, 3) = (1, 2)
    assert i1.intersect(Interval(1, 3, False, False)) == Interval(1, 2, False, False)
    # [0, 2) intersect [1, 3) = [1, 2)
    assert i1.intersect(Interval(1, 3, True, False)) == Interval(1, 2, True, False)
    assert i1.intersect(Interval(2, 3, False, False)).is_empty()
    assert not i1.intersect(Interval(-1, 0, False, True)).is_empty()

run_all()