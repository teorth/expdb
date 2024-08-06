# This file contains some unit tests 
from functions import *
from hypotheses import *
from literature import *
from exponent_pair import *
from test_polytope import *

def run_exp_pair_transform_tests():
    
    transforms = literature.list_hypotheses(hypothesis_type='Exponent pair transform')
    A_process = next(t for t in transforms if t.name == 'van der Corput A transform')
    B_process = next(t for t in transforms if t.name == 'van der Corput B transform')
    E = literature_exp_pair(frac(1,6), frac(2,3), 'Test', 2024)
    
    # B(1/6, 2/3) = (1/6, 2/3)
    BE = B_process.data.transform(E)
    assert BE.data.k == frac(1,6) and BE.data.l == frac(2,3)
    
    # A(1/6, 2/3) = (1/14, 11/14)
    AE = A_process.data.transform(E)
    assert AE.data.k == frac(1,14) and AE.data.l == frac(11,14)

    print('All exponent pair test cases finished')



def run_functions_tests():
    i1 = Interval(0, 2, True, False) # [0, 2)
    i2 = Interval(0, 2)
    assert (i1 == i2)
    # [0, 2) intersect (1, 3) = (1, 2)
    assert (i1.intersect(Interval(1, 3, False, False)) == Interval(1, 2, False, False))
    # [0, 2) intersect [1, 3) = [1, 2)
    assert (i1.intersect(Interval(1, 3, True, False)) == Interval(1, 2, True, False))
    assert (i1.intersect(Interval(2, 3, False, False)).is_empty())
    assert (not i1.intersect(Interval(-1, 0, False, True)).is_empty())
    
    print('All functions test cases finished')
    
    
def run_all():
    #run_functions_tests()
    #run_exp_pair_transform_tests()
    run_polytope_tests()
    
run_all()