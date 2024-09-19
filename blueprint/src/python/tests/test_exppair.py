import parent

from functions import *
from hypotheses import *
from literature import *
from exponent_pair import *
from reference import *


def run_exp_pair_transform_tests():

    transforms = literature.list_hypotheses(hypothesis_type="Exponent pair transform")
    A_process = next(t for t in transforms if t.name == "van der Corput A transform")
    B_process = next(t for t in transforms if t.name == "van der Corput B transform")
    E = literature_exp_pair(frac(1, 6), frac(2, 3), Reference.make("Test", 2024))

    # B(1/6, 2/3) = (1/6, 2/3)
    BE = B_process.data.transform(E)
    assert BE.data.k == frac(1, 6) and BE.data.l == frac(2, 3)

    # A(1/6, 2/3) = (1/14, 11/14)
    AE = A_process.data.transform(E)
    assert AE.data.k == frac(1, 14) and AE.data.l == frac(11, 14)

run_exp_pair_transform_tests()