# This file contains all bounds that can be directly drawn from the literature.
# For bounds that are derived from other bounds, used derived.py instead

import additive_energy as ad
from constants import *
from fractions import Fraction as frac
from hypotheses import *
import bound_beta as bbeta
import exponent_pair as ep
from bound_mu import *
import large_values as lv
import zeta_large_values as zlv
import zero_density_estimate as zd
import zero_density_energy_estimate as ze
from functions import RationalFunction as RF, Interval as Itvl
import os.path
from reference import *
from region import Region, Region_Type
from transform import Transform

# literature will contain all the bounds that are directly drawn from the literature
literature = Hypothesis_Set()

# Set up reference manager
rm = Reference_Manager(os.path.dirname(__file__) + "/../references.bib")
rm.load()

########################################################################################
# Bounds on \beta(\alpha), sorted in chronological order


# Watt (1989), "Exponential sums and the Riemann Zeta Function II" Theorem 5
def add_beta_bound_watt_1989():
    bbeta.add_beta_bound(
        literature,
        [Affine(frac(1, 2), frac(89, 560), Interval("[3/7, 1/2]"))],
        rm.get("watt_exponential_1989"),
    )


add_beta_bound_watt_1989()


# Huxley-Kolesnik (1991), "Exponential Sums and the Riemann Zeta Function III" Theorem 3
def add_beta_bound_huxley_kolesnik_1991():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(8, 22), frac(1, 22), Interval("[2/5, 1/2]")).max_with(
            [
                Affine(frac(112, 158), frac(11, 158), Interval("[2/5, 1/2]")),
                Affine(frac(17, 22), frac(1, 22), Interval("[2/5, 1/2]")),
            ]
        ),
        rm.get("huxley_exponential_1991"),
    )


add_beta_bound_huxley_kolesnik_1991()


# Huxley (1993) "Exponential sums and the Riemann Zeta Function IV" Theorem 1
def add_beta_bound_huxley_1993():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(7, 20), frac(13, 60), Interval("[0, 49/144]")).max_with(
            [Affine(frac(13, 20), frac(11, 120), Interval("[0, 49/144]"))]
        ),
        rm.get("huxley_exponential_1993"),
    )
    bbeta.add_beta_bound(
        literature,
        [Affine(frac(1, 2), frac(89, 570), Interval("[49/144, 1/2]"))],
        rm.get("huxley_exponential_1993"),
    )


add_beta_bound_huxley_1993()


# Huxley (1993) "Exponential sums and the Riemann Zeta Function IV" Theorem 3
def add_beta_bound_huxley_1993_3():
    bbeta.add_beta_bound(
        literature,
        [
            Affine(frac(94, 146), frac(13, 146), Interval("[0, 87/275]")),
            Affine(frac(191, 244), frac(11, 244), Interval("[87/275, 423/1295]")),
            Affine(frac(908, 1282), frac(89, 1282), Interval("[423/1295, 227/601]")),
            Affine(frac(173, 280), frac(29, 280), Interval("[227/601, 12/31]")),
            Affine(frac(103, 128), frac(4, 128), Interval("[12/31, 1]")),
        ],
        rm.get("huxley_exponential_1993"),
    )


add_beta_bound_huxley_1993_3()


# Huxley (1996) "Area, Lattice points and Exponential sums" Table 17.1
def add_beta_bound_huxley_1996():
    bbeta.add_beta_bound(
        literature,
        [
            Affine(frac(253, 318), frac(13, 318), Interval("[2848/12173, 161/646)")),
            Affine(frac(428, 492), frac(11, 492), Interval("[161/646, 19/74)")),
            Affine(frac(2243, 2706), frac(89, 2706), Interval("[19/74, 199/716)")),
            Affine(frac(464, 600), frac(29, 600), Interval("[199/716, 967/3428)")),
            Affine(frac(1351, 1614), frac(49, 1614), Interval("[967/3428, 120/419)")),
            Affine(frac(235, 264), frac(4, 264), Interval("[120/419, 1424/4747)")),
            Affine(frac(94, 146), frac(13, 146), Interval("[1424/4747, 87/275)")),
            Affine(frac(191, 244), frac(11, 244), Interval("[87/275, 423/1295)")),
            Affine(frac(908, 1282), frac(89, 1282), Interval("[423/1295, 227/601)")),
            Affine(frac(173, 280), frac(29, 280), Interval("[227/601, 12/31)")),
            Affine(frac(103, 128), frac(4, 128), Interval("[12/31, 356/873)")),
            Affine(frac(21, 60), frac(13, 60), Interval("[356/873, 5/12)")),
            Affine(frac(78, 120), frac(11, 120), Interval("[5/12, 49/114)")),
            Affine(frac(285, 570), frac(89, 570), Interval("[49/114, 65/114)")),
            Affine(frac(42, 120), frac(29, 120), Interval("[65/114, 7/12)")),
            Affine(frac(39, 60), frac(4, 60), Interval("[7/12, 517/873]")),
        ],
        rm.get("huxley_area_1996"),
    )


add_beta_bound_huxley_1996()


# Sargos (1995) Theorem 2.4 and Lemma 2.6
def add_beta_bound_sargos_1995():
    bbeta.add_beta_bound(
        literature,
        Affine(1 - 3*frac(4,40), frac(3,40), Interval("[0, 1]")).max_with([
            Affine(frac(7,8), 0, Interval("[0, 1]")),
            Affine(frac(1,3) + frac(4,6), -frac(1,6), Interval("[0, 1]")),
            Affine(frac(0), frac(0), Interval("[0, 1]"))
            ]),
        rm.get("sargos_points_1995")
        )
    bbeta.add_beta_bound(
        literature,
        Affine(1 - frac(4,14), frac(1,14), Interval("[0, 1]")).max_with([
            Affine(frac(5,6), 0, Interval("[0, 1]")),
            Affine(frac(1,3) + frac(4,6), -frac(1,6), Interval("[0, 1]")),
            Affine(frac(0), frac(0), Interval("[0, 1]"))
            ]),
        rm.get("sargos_points_1995")
        )

add_beta_bound_sargos_1995()



# Huxley (1996) "Area, Lattice points and Exponential sums" Table 19.2
def add_beta_bound_huxley_1996_2():
    bbeta.add_beta_bound(
        literature,
        [
            Affine(frac(103, 128), frac(4, 128), Interval("[12/31, 68682/171139)")),
            Affine(
                frac(2484, 6410), frac(1273, 6410), Interval("[68682/171139, 307/761)")
            ),
            Affine(frac(1053, 2800), frac(569, 2800), Interval("[307/761, 143/349)")),
            Affine(frac(3624, 5530), frac(491, 5530), Interval("[143/349, 263/638)")),
            Affine(
                frac(897, 1345), frac(113, 1345), Interval("[263/638, 156527/370694)")
            ),
            Affine(
                frac(88442, 134680),
                frac(11897, 134680),
                Interval("[156527/370694, 699371/1647930)"),
            ),
            Affine(
                frac(19177, 29855),
                frac(2819, 29855),
                Interval("[699371/1647930, 675/1574)"),
            ),
            Affine(
                frac(17972, 27290),
                frac(2387, 27290),
                Interval("[675/1574, 106822/246639)"),
            ),
            Affine(
                frac(285, 570),
                frac(89, 570),
                Interval("[106822/246639, 139817/246639]"),
            ),
        ],
        rm.get("huxley_area_1996"),
    )


add_beta_bound_huxley_1996_2()


# Huxley--Kolesnik (2001) "Exponential sums with a large second derivative" Theorem 1
def add_beta_bound_huxley_kolesnik_2001():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(79, 120), frac(7, 80), Interval("[2/5, 1/2]")).max_with(
            [
                Affine(frac(103, 160), frac(3, 32), Interval("[2/5, 1/2]")),
                Affine(frac(13, 40), frac(9, 40), Interval("[2/5, 1/2]")),
            ]
        ),
        rm.get("huxley_exponential_2001"),
    )


add_beta_bound_huxley_kolesnik_2001()


# Robert-Sargos (2002) "A Fourth Derivative Test for Exponential Sums
#    $$ \beta(\alpha) \leq \max\left( \alpha + \frac{1-4\alpha}{13}, -\frac{7(1-4\alpha)}{13}\right).$$
def add_beta_bound_robert_sargos_2002():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(9, 13), frac(1, 13), Interval("[0, 1]")).max_with(
            [Affine(frac(28, 13), frac(-7, 13), Interval("[0, 1]"))]
        ),
        rm.get("robert_fourth_2002"),
    )


add_beta_bound_robert_sargos_2002()


# Sargos (2003) "An analog of van der Corput’s A4-process for exponential sums"
def add_beta_bound_sargos_2003():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(204 - 8, 204), frac(1, 204), Interval("[0, 1]")).max_with(
            [Affine(frac(95 * 8, 204), frac(-95, 204), Interval("[0, 1]"))]
        ),
        rm.refs["sargos_analog_2003"],
    )

    bbeta.add_beta_bound(
        literature,
        Affine(frac(2640 - 7 * 9, 2640), frac(7, 2640), Interval("[0, 1]")).max_with(
            [Affine(frac(1001 * 9, 2640), frac(-1001, 2640), Interval("[0, 1]"))]
        ),
        rm.refs["sargos_analog_2003"],
    )


add_beta_bound_sargos_2003()


# Huxley 2005 "Exponential sums and the Riemann Zeta Function V" Proposition 1 + Theorem 1
def add_beta_bound_huxley_2005():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(59, 170), frac(37, 170), Interval("[1/3, 1/2]")).max_with(
            [Affine(frac(449, 690), frac(21, 230), Interval("[1/3, 1/2]"))]
        ),
        rm.refs["huxley_exponential_2005"],
    )


add_beta_bound_huxley_2005()


# Robert (2016) "On the fourth derivative test for exponential sums" Theorem 1
def add_beta_bound_robert_2016():
    bbeta.add_beta_bound(
        literature,
        Affine(frac(2, 3), frac(1, 12), Interval("[0, 3/7]")).max_with(
            [Affine(frac(11, 12), 0, Interval("[0, 3/7]"))]
        ),
        rm.refs["robert_fourth_2016"],
    )


add_beta_bound_robert_2016()


# Robert (2016) "On van der Corput's k-th derivative test for exponential sums" Theorem 10
def add_beta_bound_robert_2016_2(K):
    for k in range(4, K):
        bbeta.add_beta_bound(
            literature,
            Affine(
                1 - frac(k, 2 * (k - 1) * (k - 2)),
                frac(1, 2 * (k - 1) * (k - 2)),
                Interval(0, frac(k - 1, k * (k - 1) - (2 * k - 3)), True, True),
            ).max_with(
                [
                    Affine(
                        1 - frac(1, 2 * (k - 1) * (k - 2)),
                        0,
                        Interval(0, frac(k - 1, k * (k - 1) - (2 * k - 3)), True, True),
                    )
                ]
            ),
            rm.refs["robert_2016"],
        )


add_beta_bound_robert_2016_2(Constants.BETA_TRUNCATION)


# First few beta bounds in the sequence due to Heath-Brown (2017) (consistent with exponent pairs)
def add_beta_bound_heath_brown_2017(K):
    for k in range(3, K):
        bbeta.add_beta_bound(
            literature,
            Affine(
                frac(k - 2, k - 1), frac(1, k * (k - 1)), Interval("[0, 1]")
            ).max_with(
                [
                    Affine(1 - frac(1, k * (k - 1)), 0, Interval("[0, 1]")),
                    Affine(1, frac(-2, k * k * (k - 1)), Interval("[0, 1]")),
                ]
            ),
            rm.get("heathbrown_new_2017"),
        )


add_beta_bound_heath_brown_2017(Constants.BETA_TRUNCATION)


# Bourgain (2017) eqn. (3.18)
def add_beta_bound_bourgain_2017():
    bbeta.add_beta_bound(
        literature,
        [
            Affine(frac(1, 3), frac(2, 9), Interval("[1/3, 5/12)")),
            Affine(frac(2, 3), frac(1, 12), Interval("[5/12, 3/7)")),
            Affine(frac(1, 2), frac(13, 84), Interval("[3/7, 1/2]")),
        ],
        rm.get("bourgain_decoupling_2017"),
    )


add_beta_bound_bourgain_2017()



def add_beta_bound_trudgian_yang_2024():
    # Other bounds on beta are stated in the LaTeX blueprint, however they have
    # already been added to the beta bounds literature
    bbeta.add_beta_bound(
        literature,
        [
            Affine(frac(359, 414), frac(13, 414), Interval("[0, 2848/12173)")),
            Affine(frac(139, 194), frac(13, 194), Interval("[1328/4447, 104/343)")),
            Affine(frac(521, 796), frac(18, 199), Interval("[1508/3825, 62831/155153]")),
        ],
        rm.get("trudgian-yang"),
    )

add_beta_bound_trudgian_yang_2024()


########################################################################################
# Known exponent pairs in the literature. As a convention we only record bounds
# (k, l) satisfying l >= k + 1/2, since if (k, l) is an exponent pair, the
# van der Corput B process implies that (l - 1/2, k + 1/2) is also an exponent pair.
# The (k, l) values are recorded without epsilons.


def add_literature_exponent_pairs():
    literature.add_hypotheses(
        [
            # Exponent pairs on the line of symmetry l = k + 1/2
            ep.literature_exp_pair(
                frac(9, 56), frac(37, 56), rm.get("huxley_exponential_1988")
            ),
            ep.literature_exp_pair(
                frac(89, 560), frac(369, 560), rm.get("watt_exponential_1989")
            ),
            ep.literature_exp_pair(
                frac(17, 108), frac(71, 108), rm.get("huxley_exponential_1991")
            ),
            ep.literature_exp_pair(
                frac(89, 570), frac(187, 285), rm.get("huxley_exponential_1993")
            ),
            ep.literature_exp_pair(
                frac(32, 205), frac(269, 410), rm.get("huxley_exponential_2005")
            ),
            ep.literature_exp_pair(
                frac(13, 84), frac(55, 84), rm.get("bourgain_decoupling_2017")
            ),
            # Exponent pairs (k, l) with l > k + 1/2.
            # Exponent pair derived from Bombieri-Iwaniec method
            ep.literature_exp_pair(
                frac(2, 13), frac(35, 52), rm.get("huxley_watt_exponential_1990")
            ),
            ep.literature_exp_pair(
                frac(6299, 43860), frac(29507, 43860), rm.get("huxley_area_1996")
            ),
            ep.literature_exp_pair(
                frac(771, 8116), frac(1499, 2029), rm.get("sargos_points_1995")
            ),
            ep.literature_exp_pair(
                frac(21, 232), frac(173, 232), rm.get("sargos_points_1995")
            ),
            ep.literature_exp_pair(
                frac(1959, 21656), frac(16135, 21656), rm.get("sargos_points_1995")
            ),
            ep.literature_exp_pair(
                frac(516247, 6629696),
                frac(5080955, 6629696),
                rm.get("huxley_exponential_2001"),
            ),
            # Exponent pairs derived from k-th derivative tests
            ep.literature_exp_pair(
                frac(1, 13), 1 - frac(3, 13), rm.get("robert_fourth_2002")
            ),
            ep.literature_exp_pair(
                frac(1, 204), 1 - frac(7, 204), rm.get("sargos_analog_2003")
            ),
            ep.literature_exp_pair(
                frac(1, 360), 1 - frac(8, 360), rm.get("robert_2002")
            ),
            ep.literature_exp_pair(
                frac(7, 2640), 1 - frac(56, 2640), rm.get("sargos_analog_2003")
            ),
            ep.literature_exp_pair(
                frac(1, 716), 1 - frac(9, 716), rm.get("sargos_analog_2003")
            ),
            ep.literature_exp_pair(
                frac(1, 649), 1 - frac(9, 649), rm.get("Robert_Sargos_2001")
            ),
            ep.literature_exp_pair(
                frac(7, 4540), 1 - frac(63, 4540), rm.get("robert_2002")
            ),
            ep.literature_exp_pair(
                frac(1, 615), 1 - frac(9, 615), rm.get("robert_2002b")
            ),
            ep.literature_exp_pair(
                frac(1, 915), 1 - frac(10, 915), rm.get("robert_2002b")
            ),
        ]
    )


add_literature_exponent_pairs()


# Exponent pairs in the sequence due to Huxley Ch. 17 (1996)
def add_huxley_exponent_pairs(M):
    for m in range(1, M):
        literature.add_hypothesis(
            literature_exp_pair(
                frac(169, 1424 * pow(2, m) - 338),
                1 - frac(169, 1424 * pow(2, m) - 338) * (m + frac(1577, 712)),
                rm.get("huxley_area_1996"),
            )
        )


add_huxley_exponent_pairs(Constants.EXP_PAIR_TRUNCATION)

# Exponent pairs in the sequence due to Heath-Brown (1986), which
# appear in the end-of-chapter notes in Titchmarsh Ch. 6
# TODO: These are removed for now because Fraction does not support arithmetic with mpf
# objects.  one could round the exponents up to suitable precision (taking care to stay inside the triangle)

# for m in range(3, Constants.EXP_PAIR_TRUNCATION):
#     literature.add_hypothesis(literature_exp_pair(1/(25 * m * m * mp.log(m)),
#                                          1 - 1/(25 * m * m * mp.log(m)),
#                                          rm.get('titchmarsh_theory_1986')))


# Exponent pairs in the sequence due to Heath-Brown (2017)
def add_heath_brown_exponent_pairs(M):
    for m in range(3, M):
        literature.add_hypothesis(
            literature_exp_pair(
                frac(2, (m - 1) * (m - 1) * (m + 2)),
                1 - frac(3 * m - 2, m * (m - 1) * (m + 2)),
                rm.get("heathbrown_new_2017"),
            )
        )


add_heath_brown_exponent_pairs(Constants.EXP_PAIR_TRUNCATION)


def add_exp_pairs_all(hypothesis_list):
    hypothesis_list.add_hypothesis(trivial_exp_pair)
    hypothesis_list.add_hypotheses(literature.list_hypotheses("Exponent pair"))


def add_exp_pairs_up_to(hypothesis_list, year):
    hypothesis_list.add_hypothesis(trivial_exp_pair)
    hypothesis_list.add_hypotheses(literature.list_hypotheses("Exponent pair", year))


###############################################################################
# List of exponent pair transforms (maps of exponent pair -> exponent pair) in
# the literature, e.g. the van der Corput A/B transforms.

# van der Corput A transform (weyl-van der Corput inequality)
def A_transform_function(hypothesis):
    pair = hypothesis.data
    return derived_exp_pair(
        pair.k / (2 * pair.k + 2),
        frac(1, 2) + pair.l / (2 * pair.k + 2),
        f'Follows from "{hypothesis.name}" and taking the van der Corput A transform',
        {hypothesis, A_transform_hypothesis},
    )
A_transform_hypothesis = Hypothesis(
        "van der Corput A transform",
        "Exponent pair transform",
        Transform("van der Corput A transform", A_transform_function),
        "See [van der Corput, 1920]",
        Reference.make("Weyl--van der Corput", 1920),
    )


def B_transform_function(hypothesis):  # van der Corput B transform (Poisson summation)
    pair = hypothesis.data
    return derived_exp_pair(
        pair.l - frac(1, 2),
        pair.k + frac(1, 2),
        f'Follows from "{hypothesis.name}" and taking the van der Corput B transform',
        {hypothesis, B_transform_hypothesis},
    )
B_transform_hypothesis = Hypothesis(
        "van der Corput B transform",
        "Exponent pair transform",
        Transform("van der Corput B transform", B_transform_function),
        "See [van der Corput, 1920]",
        Reference.make("van der Corput", 1920),
    )

def C_transform_function(hypothesis):  # Sargos 2003 transform
    pair = hypothesis.data
    return derived_exp_pair(
        pair.k / (12 * (1 + 4 * pair.k)),
        (11 * (1 + 4 * pair.k) + pair.l) / (12 * (1 + 4 * pair.k)),
        f'Follows from "{hypothesis.name}" and taking the Sargos C transform',
        {hypothesis, C_transform_hypothesis},
    )
C_transform_hypothesis = Hypothesis(
        "Sargos C transform",
        "Exponent pair transform",
        Transform("Sargos C transform", C_transform_function),
        "See [Sargos, 2003]",
        rm.get("sargos_analog_2003"),
    )

# Sargos 1995 transform: exp pair -> list of beta bounds
def D_transform_function(hypothesis):

    if hypothesis.hypothesis_type != "Exponent pair":
        raise ValueError("Parameter hypothesis must be of type Exponent pair")

    k = hypothesis.data.k
    l = hypothesis.data.l
    domain = Interval(0, frac(1,2), True, True)
    pieces = Affine(
        (6 * k + 5 * l + 2) / (2 * (5 * k + 3 * l + 2)),
        (5 * k + l + 2) / (8 * (5 * k + 3 * l + 2)),
        domain,
    ).max_with([Affine(frac(2, 3), frac(1, 12), domain)])

    return [
        bbeta.derived_bound_beta(
                p,
                f'Follows from "{hypothesis.name}" and taking the Sargos D transform',
                {hypothesis, D_transform_hypothesis},
            )
        for p in pieces
        ]
D_transform_hypothesis = Hypothesis(
        "Sargos D transform",
        "Exponent pair to beta bound transform",
        Transform("Sargos D transform", D_transform_function),
        'See [Sargos, 1995] Theorem 7.1',
        rm.get('sargos_points_1995'),
        )


literature.add_hypothesis(A_transform_hypothesis)
literature.add_hypothesis(B_transform_hypothesis)
literature.add_hypothesis(C_transform_hypothesis)
literature.add_hypothesis(D_transform_hypothesis)


########################################################################################
# We now list the known upper bounds on $\mu$ in the literature.
# This list can be updated as new bounds are established, or old bounds not currently present in the list are added.


# All bounds in the literature, including older obsolete bounds - listed for historical
# or pedagogical purposes, or for trying to find the minimal inputs required to establish a given bound.
def add_literature_bounds_mu():
    literature.add_hypotheses(
        [
            # Bounds on the critical line
            literature_bound_mu(
                frac(1, 2), frac(1, 6), rm.get("hardy_littlewood_1923")
            ),
            literature_bound_mu(frac(1, 2), frac(193, 988), rm.get("walfisz_1924")),
            literature_bound_mu(
                frac(1, 2), frac(27, 164), rm.get("titchmarsh_van_1931")
            ),
            literature_bound_mu(
                frac(1, 2), frac(229, 1392), rm.get("phillips_zeta_1933")
            ),
            literature_bound_mu(
                frac(1, 2), frac(19, 116), rm.get("titchmarsh_order_1942")
            ),
            literature_bound_mu(frac(1, 2), frac(15, 92), rm.get("min_on_1949")),
            literature_bound_mu(
                frac(1, 2), frac(6, 37), rm.get("haneke_verscharfung_1963")
            ),
            literature_bound_mu(
                frac(1, 2), frac(173, 1067), Reference.make("Kolesnik", 1973)
            ),
            literature_bound_mu(
                frac(1, 2), frac(35, 216), rm.get("kolesnik_order_1982")
            ),
            literature_bound_mu(
                frac(1, 2), frac(139, 858), Reference.make("Kolesnik", 1985)
            ),
            literature_bound_mu(frac(1, 2), frac(9, 56), rm.get("bombieri_order_1986")),
            literature_bound_mu(
                frac(1, 2), frac(89, 560), rm.get("watt_exponential_1989")
            ),
            literature_bound_mu(
                frac(1, 2), frac(17, 108), rm.get("huxley_exponential_1991")
            ),
            literature_bound_mu(
                frac(1, 2), frac(89, 570), rm.get("huxley_exponential_1993")
            ),
            literature_bound_mu(
                frac(1, 2), frac(32, 205), rm.get("huxley_exponential_2005")
            ),
            literature_bound_mu(
                frac(1, 2), frac(13, 84), rm.get("bourgain_decoupling_2017")
            ),
            # bounds off the critical line
            literature_bound_mu(
                frac(1934, 3655), frac(6299, 43860), rm.get("huxley_area_1996")
            ),
            literature_bound_mu(
                frac(49, 51), frac(1, 204), rm.get("sargos_analog_2003")
            ),
            literature_bound_mu(
                frac(361, 370), frac(1, 370), rm.get("sargos_analog_2003")
            ),
            literature_bound_mu(
                frac(3, 5), frac(1409, 12170), rm.get("Lelechenko_linear_2014")
            ),
            literature_bound_mu(
                frac(4, 5), frac(3, 71), rm.get("Lelechenko_linear_2014")
            ),
            literature_bound_mu(
                frac(11, 15), frac(1, 15), rm.get("demeter_small_2020")
            ),
        ]
    )


add_literature_bounds_mu()

# van der Corput sequence of mu bounds
literature.add_hypotheses(
    literature_bound_mu(
        1 - frac(n, pow(2, n) - 2),
        frac(1, pow(2, n) - 2),
        Reference.make("van der Corput", 1920),
    )
    for n in range(4, Constants.EXP_PAIR_TRUNCATION)
)

# Hardy-Littlewood sequence of mu bounds
literature.add_hypotheses(
    literature_bound_mu(
        1 - frac(1, pow(2, n - 1)),
        frac(1, (n + 1) * pow(2, n - 1)),
        Reference.make("Hardy--Littlewood", "Unknown date"),
    )  # TODO: find the exact year
    for n in range(4, Constants.EXP_PAIR_TRUNCATION)
)

########################################################################################
# Methods for adding bounds on to an ambient hypothesis set.  May be better to just call the Hypothesis_Set methods directly.


# use this method if you just want to add the "classic" bounds on mu
def add_bounds_mu_classic(hypothesis_set):
    hypothesis_set.add_hypothesis(
        literature.find_hypothesis(data=Bound_mu(frac(1, 2), frac(1, 6)))
    )


# use this method if you want to add *all* the (unconditional) bounds on mu
def add_bounds_mu_all(hypothesis_set):
    hypothesis_set.add_hypotheses(literature.list_hypotheses("Upper bound on mu"))


# use this method to add all known bounds on mu up to a given year
def add_bounds_mu_as_of(hypothesis_set, year):
    hypothesis_set.add_hypotheses(literature.list_hypotheses("Upper bound on mu", year))


########################################################################################
# Chapter 7: List of large value estimates in the literature

# Huxley large values theorem: LV(s, t) \leq max(2 - 2s, 4 - 6s + t)
def add_huxley_large_values_estimate():
    literature.add_hypothesis(
        lv.literature_bound_LV_max([[2, -2, 0], [4, -6, 1]], rm.get("Huxley"))
    )

# Heath-Brown large values theorem: LV(s, t) \leq max(2 - 2s, 10 - 13s + t)
def add_heath_brown_large_values_estimate():
    literature.add_hypothesis(
        lv.literature_bound_LV_max(
            [[2, -2, 0], [10, -13, 1]], rm.get("heathbrown_zero_1979")
        )
    )

# Jutila large values theorem:
# LV(s, t) \leq \max(2 - 2s, (4 - 2/k) - (6 - 2/k)s + t, 6k - 8ks + t)
def add_jutila_large_values_estimate(K):
    for k in range(1, K):
        literature.add_hypothesis(
            lv.literature_bound_LV_max(
                [[2, -2, 0], [4 - frac(2, k), -(6 - frac(2, k)), 1], [6 * k, -8 * k, 1]],
                rm.get("jutila_zero_density_1977"),
                params=f" with k = {k}",
            )
        )

# Bourgain large-values theorem with optimal choices of \alpha_1, \alpha_2 given by
# Given in Table 7.1 of the LaTeX blueprint
def add_bourgain_large_values_estimate():
    literature.add_hypotheses(
        lv.get_optimized_bourgain_lv_estimate(rm.get("bourgain_large_2000"))
    )

# Guth-Maynard (2024) large values theorem: LV(s, t) \leq max(2 - 2s, 18/5 - 4s, t + 12/5 - 4s)
def add_guth_maynard_large_values_estimate():
    literature.add_hypothesis(
        lv.literature_bound_LV_max(
            [[2, -2, 0], [frac(18, 5), -4, 0], [frac(12, 5), -4, 1]], rm.get("guth-maynard")
        )
    )

add_huxley_large_values_estimate()
add_heath_brown_large_values_estimate()
add_jutila_large_values_estimate(Constants.LARGE_VALUES_TRUNCATION)
add_bourgain_large_values_estimate()
add_guth_maynard_large_values_estimate()

########################################################################################
# Chapter 8: List of zeta large value estimates in the literature
literature.add_hypothesis(
    zlv.literature_bound_ZLV_Max([[6, -12, 2]], rm.get("heathbrown_twelfth_1978"))
)


#################################################################
# Chapter 10: List of large value energy region theorems from the literature

def add_lver_heath_brown_1979():
    region = ad.union_of_halfplanes(
        [
            [2, 0, 0, 1, 0, -1],                    # 2 + rho - s >= 0
            [1, 0, 0, 2, 0, -1],                    # 1 + 2 * rho - s >= 0
            [1, 0, frac(1,2), frac(5,4), 0, -1],    # 1 + 1/2 * tau + 5/4 * rho - s >= 0
        ],
        ad.Large_Value_Energy_Region.default_constraints()
    )
    literature.add_hypothesis(
        ad.literature_large_value_energy_region(
            region,
            rm.get("heathbrown_large_1979"),
            params= " 1"
        )
    )
add_lver_heath_brown_1979()


def add_lver_ivic_1985():
    # divide into two cases: tau <= 1 and tau >= 1
    rect1 = ad.Large_Value_Energy_Region.default_constraints()
    rect1.append([1, 0, -1, 0, 0, 0]) # tau <= 1
    region = ad.union_of_halfplanes(
        [
            [2, 0, 0, 1, 0, -1],    # 2 + rho - s >= 0
            [1, 0, 0, 2, 0, -1]     # 1 + 2 * rho - s >= 0
        ],
        rect1
    )

    rect2 = ad.Large_Value_Energy_Region.default_constraints()
    rect2.append([-1, 0, 1, 0, 0, 0]) # tau >= 1
    region.child.append(Region(Region_Type.POLYTOPE, Polytope(rect2)))

    literature.add_hypothesis(
        ad.literature_large_value_energy_region(
            region,
            rm.get("ivic")
        )
    )
add_lver_ivic_1985()

# Implementation of "Simplified Heath-Brown relation" Corollary 10.19
def add_lver_heath_brown_1979b():
    # divide into two cases
    rect1 = ad.Large_Value_Energy_Region.default_constraints()
    rect1.append([frac(3,2), 0, -1, 0, 0, 0]) # tau <= 3/2
    region = ad.union_of_halfplanes(
        [
            [1, -2, 0, 3, -1, 0],                   # 1 - 2sigma + 3rho - rho* >= 0
            [4, -4, 0, 1, -1, 0],                   # 4 - 4sigma + rho - rho* >= 0
            [frac(3,2), -2, 0, frac(5,2), -1, 0]    # 3/2 - 2sigma + 5/2rho - rho* >= 0
        ],
        rect1
    )

    rect2 = ad.Large_Value_Energy_Region.default_constraints()
    rect2.append([-frac(3,2), 0, 1, 0, 0, 0]) # tau >= 3/2
    region.child.append(Region(Region_Type.POLYTOPE, Polytope(rect2)))

    literature.add_hypothesis(
        ad.literature_large_value_energy_region(
            region,
            rm.get("heathbrown_zero_1979"),
            params=" 2"
        )
    )
add_lver_heath_brown_1979b()

# Implementation of second Heath-Brown relation (Lemma 10.20)
def add_lver_heath_brown_1979c(K):
    rect = ad.Large_Value_Energy_Region.default_constraints()
    for k in range(1, K):
        region = ad.union_of_halfplanes(
            [
                [2, -2, 0, -1, 0, 0],                           # 2 - 2sigma - rho >= 0
                [frac(3*k,4), -k, frac(1,4), -1, frac(1,4), 0], # 3k/4 - k sigma + tau/4 - rho + rho*/4 >= 0
                [frac(k,2), -k, frac(k,4), -1, frac(1,4), 0]    # k/2 - k sigma + k/4 tau - rho + rho*/4 >= 0
            ],
            rect
        )
        literature.add_hypothesis(
            ad.literature_large_value_energy_region(
                region,
                rm.get("heathbrown_zero_1979"),
                params=f" 3 with k = {k}"
            )
        )
add_lver_heath_brown_1979c(5)

def add_lver_guth_maynard_2024a(K):
    rect = ad.Large_Value_Energy_Region.default_constraints()
    for k in range(1, K):
        polys = []
        # 2 - 2sigma - rho >= 0
        polys.append(Polytope(rect + [[2, -2, 0, -1, 0, 0]]))
        # 1 - 2sigma + S_1/3 - rho >= 0 
        # -7/3 - 2sigma - rho >= 0
        polys.append(Polytope(rect + [[-frac(7,3), -2, 0, -1, 0, 0]]))
        
        # TODO: complete the list of constraints
        raise NotImplementedError()
#add_lver_guth_maynard_2024a(10)

def add_lver_guth_maynard_2024b():
    rect = ad.Large_Value_Energy_Region.default_constraints()
    region = ad.union_of_halfplanes(
        [
            [0, -2, 0, 1, -1, 1]    # -2sigma + rho - rho* + s >= 0
        ],
        rect
    )
    literature.add_hypothesis(
        ad.literature_large_value_energy_region(
            region,
            rm.get("guth-maynard"),
            params=" 2"
        )
    )
add_lver_guth_maynard_2024b()

def add_lver_guth_maynard_2024c():
    rect1 = ad.Large_Value_Energy_Region.default_constraints()
    rect1.append([frac(4,3), 0, -1, 0, 0, 0]) # 1 <= tau <= 4/3
    rect1.append([-1, 0, 1, 0, 0, 0]) 
    region = ad.union_of_halfplanes(
        [
            [4, -4, 0, 1, -1, 0],                   # 4 - 4sigma + rho - rho* >= 0
            [1, -2, frac(1,4), frac(21,8), -1, 0],  # 1 - 2sigma + tau/4 + 21/8 rho - rho* >= 0 
            [1, -2, 0, 3, -1, 0]                    # 1 - 2sigma + 3rho - rho* >= 0
        ],
        rect1
    )
    rect2 = ad.Large_Value_Energy_Region.default_constraints()
    rect2.append([-frac(4,3), 0, 1, 0, 0, 0]) # tau >= 4/3
    region.child.append(Region(Region_Type.POLYTOPE, Polytope(rect2)))

    rect3 = ad.Large_Value_Energy_Region.default_constraints()
    rect3.append([1, 0, -1, 0, 0, 0]) # tau <= 1
    region.child.append(Region(Region_Type.POLYTOPE, Polytope(rect3)))
    
    literature.add_hypothesis(
        ad.literature_large_value_energy_region(
            region,
            rm.get("guth-maynard"),
            params=" 3"
        )
    )
add_lver_guth_maynard_2024c()


########################################################################################
# Chapter 11: List of zero-density estimates in the literature, in chronological order

# Carlson (1921) Uber die Nullstellen der Dirichletschen Reihen und der Riemannschen ζ-Funktion
def add_zero_density_carlson_1921():
    zd.add_zero_density(
        literature, "4*x", Itvl(frac(1, 2), 1), rm.get("carlson_uber_1921")
    )
add_zero_density_carlson_1921()

# Ingham (1940) On the estimation of N(σ,T)
def add_zero_density_ingham_1940():
    zd.add_zero_density(
        literature, "3/(2 - x)", Itvl(frac(1, 2), 1), rm.get("ingham_estimation_1940")
    )
add_zero_density_ingham_1940()

# Montgomery (1971) Topics in Multiplicative Number Theory
def add_zero_density_montgomery_1971():
    zd.add_zero_density(
        literature,
        "1600 * (1 - x) ** (1/2)",
        Itvl(frac(1, 2), 1),
        Reference.make("Montgomery", 1971),
    )
add_zero_density_montgomery_1971()

# M. N. Huxley (1972) On the difference between consecutive primes, Invent. Math., 15, pages 164--170
def add_zero_density_huxley_1972():
    zd.add_zero_density(
        literature, "3/(3 * x - 1)", Itvl(frac(1, 2), 1), rm.get("Huxley")
    )
add_zero_density_huxley_1972()

# M. N. Huxley (1973) Large values of Dirichlet polynomials, Acta Arithmetica Volume: 24, pages 329--346
zd.add_zero_density(
    literature,
    "39/(115 * x - 75)",
    Itvl(frac(55, 67), frac(189, 230)),
    rm.get("huxley_large_1973"),
)
zd.add_zero_density(
    literature, "2", Itvl(frac(189, 230), frac(75, 89)), rm.get("huxley_large_1973")
)
zd.add_zero_density(
    literature,
    "48/(37 * (2 * x - 1))",
    Itvl(frac(61, 74), frac(75, 89)),
    rm.get("huxley_large_1973"),
)

# M. N. Huxley (1975) Large values of Dirichlet polynomials II, Acta Arithmetica Volume: 27, Issue: 1, pages 159--170
zd.add_zero_density(
    literature, "3/(2 * x)", Itvl(frac(37, 42), 1), rm.get("huxley_large_1975a")
)
zd.add_zero_density(
    literature,
    "48/(37 * (2 * x - 1))",
    Itvl(frac(61, 74), frac(37, 42)),
    rm.get("huxley_large_1975a"),
)
zd.add_zero_density(
    literature, "2", Itvl(frac("0.80119"), 1), rm.get("huxley_large_1975a")
)

# M. N. Huxley (1975) Large values of Dirichlet polynomials III, Acta Arithmetica Volume: 26, pages 435--444
zd.add_zero_density(literature, "2", Itvl(frac(4, 5), 1), rm.get("huxley_large_1975b"))

# D. R. Heath-Brown (1979) Zero Density Estimates for the Riemann Zeta-Function and Dirichlet L-Functions, J. Lond. Math. Soc. s2-19, pages 221--232
def add_zero_density_heathbrown_1979():
    zd.add_zero_density(
        literature, "4/(4 * x - 1)", Itvl("[25/28, 1]"), rm.get("heathbrown_zero_1979")
    )
    zd.add_zero_density(
        literature, "3/(10 * x - 7)", Itvl("[3/4, 25/28)"), rm.get("heathbrown_zero_1979"),
    )
    zd.add_zero_density(
        literature, "9/(7 * x - 1)", Itvl("[11/14, 1)"), rm.get("heathbrown_zero_1979")
    )
add_zero_density_heathbrown_1979()

# A. Ivic (1979) A note on the zero-density estimates for the zeta function, Archiv der Mathematik Volume: 33, pages 155--164
zd.add_zero_density(
    literature, "6/(5 * x - 1)", Itvl(frac(67, 87), 1), rm.get("ivic_note_1979")
)
zd.add_zero_density(
    literature,
    "3/(34 * x - 25)",
    Itvl(frac(28, 37), frac(74, 95)),
    rm.get("ivic_note_1979"),
)
zd.add_zero_density(
    literature, "9/(7 * x - 1)", Itvl(frac(74, 95), 1), rm.get("ivic_note_1979")
)
zd.add_zero_density(
    literature, "3/(2 * x)", Itvl(frac(4, 5), 1), rm.get("ivic_note_1979")
)
zd.add_zero_density(
    literature, "68/(98 * x - 47)", Itvl(frac(115, 166), 1), rm.get("ivic_note_1979")
)

# A. Ivic (1980) Exponent pairs and the zeta function of Riemann, Studia Sci. Math. Hung. Volume: 15, pages 157--181
def add_zero_density_ivic_1980():
    zd.add_zero_density(
        literature, "3/(2 * x)", Itvl(frac(3831, 4791), 1), Reference.make("Ivic", 1980)
    )
    zd.add_zero_density(
        literature, "2", Itvl(frac(11, 14), 1), Reference.make("Ivic", 1980)
    )
    zd.add_zero_density(
        literature, "9/(7 * x - 1)", Itvl(frac(41, 53), 1), Reference.make("Ivic", 1980)
    )
    zd.add_zero_density(
        literature, "6/(5 * x - 1)", Itvl(frac(13, 17), 1), Reference.make("Ivic", 1980)
    )
    zd.add_zero_density(
        literature, "4/(2 * x + 1)", Itvl(frac(17, 18), 1), Reference.make("Ivic", 1980)
    )
    zd.add_zero_density(
        literature, "24/(30 * x - 11)", Itvl(frac(155, 174), 1), Reference.make("Ivic", 1980),
    )
add_zero_density_ivic_1980()

# M. Jutila (1982) Zeros of the zeta-function near the critical line, Studies of Pure Mathematics, to the Memory of Paul Tur\'an, Birkha\"user Verlag, Basel-Stuttgart
zd.add_zero_density(
    literature,
    "(3/2 - x) / (1 - x)",
    Itvl(frac(155, 174), 1),
    rm.get("jutila_zeroes_1983"),
)

# A. Ivic (1984) A zero-density theorem for the Riemann zeta-function, Тр. МИАН СССP, том 163, pages 85--89
zd.add_zero_density(
    literature,
    "3/(7 * x - 4)",
    Itvl(frac(3, 4), frac(10, 13)),
    rm.get("ivic_zero_1984"),
)
zd.add_zero_density(
    literature, "9/(8 * x - 2)", Itvl(frac(10, 13), 1), rm.get("ivic_zero_1984")
)

# TODO add reference
zd.add_zero_density(
    literature,
    "15/(22 * x - 10)",
    Itvl(frac(10, 13), frac(5, 6)),
    Reference.make("Ivic", 1984),
)

# A. Ivic (1984) The Riemann zeta-function (11.76, 11.77)
# For k = 2, the estimate is already contained in 
# (A. Ivic (1980) Exponent pairs and the zeta function of Riemann, Studia Sci. Math. Hung. Volume: 15, pages 157--181)
# For k > 2, the result depends on our choice of exponent pair. Based on 
# our current knowledge, in the case of k = 3 the lower limit on sigma is 41/53, 
# and for all higher k the lower limit is given by 
# (9 * k**2 - 4 * k + 2)/(12 * k**2 - 6 * k + 2)
def add_zero_density_ivic_1984():
    for k in range(3, 100):
        if k == 3:
            sigma_lower = frac(41,53)
        else:
            sigma_lower = min(
                frac(6 * k**2 - 5 * k + 2, 8 * k**2 - 7 * k + 2),
                frac(9 * k**2 - 4 * k + 2, 12 * k**2 - 6 * k + 2)
            )
        if sigma_lower > Constants.ZERO_DENSITY_SIGMA_LIMIT:
            break
        zd.add_zero_density(
            literature, f"{3*k}/({3*k-2} * x + {2-k})", Itvl(sigma_lower, 1), rm.get("ivic"),
        )
add_zero_density_ivic_1984()

# B. Conrey (1989)
zd.add_zero_density(
    literature,
    "(1 + (1 - 2 * x) * 4/7) / (1 - x)",
    Itvl(frac(1, 2), 1),
    rm.get("conrey_at_1989"),
)

# J. Bourgain (1995) Remarks on Halasz–Montgomery type inequalities, in: Geometric aspects of functional analysis (Israel, 1992–1994),
# Oper. Theory Adv. Appl. 77, Birkh¨auser, Basel, pages 25--39.
def add_zero_density_bourgain_1995():
    zd.add_zero_density(
        literature, "4 / (30 * x - 25)", Itvl(frac(15, 16), 1), rm.get("bourgain_remarks_1995"),
    )
    zd.add_zero_density(
        literature, "2 / (7 * x - 5)", Itvl(frac(17, 19), 1), rm.get("bourgain_remarks_1995"),
    )
add_zero_density_bourgain_1995()

# J. Bourgain (2000) On large values estimates for Dirichlet polynomials and the density hypothesis for the Riemann zeta function,
# Internat. Math. Res. Notices, no. 3, pages 133--146.
def add_zero_density_bourgain_2000():
    zd.add_zero_density(
        literature, "2", Itvl(frac(25, 32), 1), rm.get("bourgain_large_2000")
    )
add_zero_density_bourgain_2000()

# K. Ford (2002) Vinogradov’s integral and bounds for the Riemann zeta function, Proc. London Math. Soc. (3) 85, no. 2, pages 565--633
zd.add_zero_density(
    literature, "58.05 * (1 - x) ** (1/2)", Itvl(frac(1, 2), 1), rm.get("FordZeta")
)

# J. Bourgain (2002) On the distribution of Dirichlet sums II, in: Number theory for the millennium, I (Urbana, IL, 2000), pages 87--109, A K Peters, Ltd., Natick, MA
def add_zero_density_bourgain_2002():
    zd.add_zero_density(
        literature, "3/(2 * x)", Itvl(frac(3734, 4694), 1), Reference.make("Bourgain", 2002)
    )
add_zero_density_bourgain_2002()

# D. R. Heath-Brown (2017) Proc. Steklov Inst. Math. 296, no. 1, pages 88--103
zd.add_zero_density(
    literature,
    "6.42 * (1 - x) ** (1/2)",
    Itvl(frac(9, 10), 1),
    rm.get("heathbrown_new_2017"),
)

#  J. Pintz (2023) Density theorems for Riemann’s zeta-function near the line \Re s = 1. Acta Arithmetica 208.1, pages 1--13
def add_zero_density_pintz_2023():
    zd.add_zero_density(
        literature,
        "3/(24 * x - 20)",
        Itvl(frac(23, 24), frac(39, 40)),
        rm.get("pintz_density_2023"),
    )
    zd.add_zero_density(
        literature,
        "2/(15 * x - 12)",
        Itvl(frac(39, 40), frac(41, 42)),
        rm.get("pintz_density_2023"),
    )
    zd.add_zero_density(
        literature,
        "3/(40 * x - 35)",
        Itvl(frac(41, 42), frac(59, 60)),
        rm.get("pintz_density_2023"),
    )
    for k in range(6, 100):
        sigma_lower = 1 - frac(1, 2 * k * (k - 1))
        if sigma_lower > Constants.ZERO_DENSITY_SIGMA_LIMIT:
            break
        zd.add_zero_density(
            literature,
            f"3 / ({k} * (1 - 2 * ({k} - 1) * (1 - x)))",
            Itvl(sigma_lower, 1 - frac(1, 2 * k * (k + 1))),
            rm.get("pintz_density_2023"),
        )
add_zero_density_pintz_2023()

# B. Chen, G. Debruyne and J. Vindas (2024) On the density hypothesis for L-functions associated with holomorphic cusp forms. Revista Matematica Iberoamericana
def add_zero_density_chen_debruyne_vindas_2024():
    zd.add_zero_density(
        literature,
        "24 / (30 * x - 11)",
        Itvl(frac(279, 314), frac(17, 18)),
        rm.get("chen_debruyne_vindas_density_2024"),
    )
add_zero_density_chen_debruyne_vindas_2024()

# L. Guth, J. Maynard (2024) New large value estimates for Dirichlet polynomials
def add_zero_density_guth_maynard_2024():
    zd.add_zero_density(
        literature, "15 / (3 + 5 * x)", Itvl(frac(1, 2), 1), rm.get("guth-maynard")
    )
add_zero_density_guth_maynard_2024()


########################################################################################
# Chapter 12: List of zero-density energy estimates in the literature
def add_zero_density_energy_heath_brown_1979():
    ref = rm.get("heath_brown_consecutive_II")
    literature.add_hypothesis(
        ze.literature_zero_density_energy_estimate(
            "(10 - 11 * x) / ((2 - x) * (1 - x))",
            Interval(frac(1,2), frac(2,3), True, True),
            ref
        )
    )
    literature.add_hypothesis(
        ze.literature_zero_density_energy_estimate(
            "(18 - 19 * x) / ((4 - 2 * x) * (1 - x))",
            Interval(frac(2,3), frac(3,4), True, True),
            ref
        )
    )
    literature.add_hypothesis(
        ze.literature_zero_density_energy_estimate(
            "12 / (4 * x - 1)",
            Interval(frac(3,4), frac(1), True, True),
            ref
        )
    )
add_zero_density_energy_heath_brown_1979()
