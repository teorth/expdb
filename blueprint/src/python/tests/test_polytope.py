import parent
from functions import *

def get_unit_square():
    return Polytope(
        [
            [1, -1, 0],  # x <= 1
            [0, 1, 0],  # x >= 0
            [1, 0, -1],  # y <= 1
            [0, 0, 1],  # y >= 0
        ]
    )

def run_polytope_tests():

    p = get_unit_square()
    assert p.polyhedron.rep_type == cdd.RepType.INEQUALITY

    assert p.contains([frac(1, 2), frac(1, 2)])

    assert p.get_centroid() == [frac(1, 2), frac(1, 2)]

    # Intersects with the line y = 1 - x, i.e. -1 + x + y = 0
    assert p.intersects(Hyperplane([-1, 1, 1]))

    # Does not intersect with the line y = 3 - x
    assert not p.intersects(Hyperplane([-3, 1, 1]))

    # Polytope is not empty
    assert not p.is_empty()

    # An empty polytope
    p2 = Polytope(
        [
            [1, -1, 0],  # x <= 1
            [0, 0, 1],  # y >= 0
            [-2, 1, -1],  # y <= x - 2
        ]
    )
    assert p2.is_empty()

    # A single point (x, y) = (1, 0)
    p3 = Polytope(
        [
            [1, -1, 0],  # x <= 1
            [0, 0, 1],  # y >= 0
            [-1, 1, -1],  # y <= x - 1
        ]
    )
    assert not p3.is_empty()
    # however, if we are not including the boundary then p3 is empty
    assert p3.is_empty(include_boundary=False)

    p4 = Polytope(
        [
            [2, -1, 0],  # x <= 2
            [-1, 1, 0],  # x >= 1
            [1, 0, -1],  # y <= 1
            [0, 0, 1],  # y >= 0
        ]
    )

    # Intersection operation
    assert set(p.intersect(p4).get_vertices()) == {(1, 1), (1, 0)}

    p_degen = Polytope([[1, -1, 0], [-1, 1, 0], [0, 0, 1]])  # x = 1, y >= 0
    p_degen.compute_V_rep()
    assert set(p_degen.vertices) == {(1, 0)}  # Single vertex at the point (1, 0)
    assert set(p_degen.rays) == {
        (0, 1)
    }  # Ray going to infinity in the direction (0, 1)
    assert p4.intersect(p_degen) == p_degen.intersect(p4)

    # Union operation
    assert Polytope.try_union([p, p4]) == Polytope(
        [
            [2, -1, 0],  # x <= 2
            [0, 1, 0],  # x >= 0
            [1, 0, -1],  # y <= 1
            [0, 0, 1],  # y >= 0
        ]
    )

    assert Polytope.try_union([p, p3]) == p

    p5 = Polytope(
        [
            [1, -1, -1],  # 1 - x - y >= 0
            [1, 1, -1],  # 1 + x - y >= 0
            [1, -1, 1],  # 1 - x + y >= 0
            [1, 1, 1],  # 1 + x + y >= 0
        ]
    )
    assert set(Polytope.try_union([p5, p]).get_vertices()) == {
        (-1, 0),
        (0, 1),
        (1, 1),
        (1, 0),
        (0, -1),
    }

def run_emptiness_tests():
    # Square
    p1 = Polytope.rect((0, 1), (0, 1))
    assert not p1.is_empty(include_boundary=True)
    assert not p1.is_empty(include_boundary=False)

    # Line
    p2 = Polytope.rect((0, 0), (0, 1))
    assert not p2.is_empty(include_boundary=True)
    assert p2.is_empty(include_boundary=False)

    # Trapezium
    p3 = Polytope([
        [0, 0, 1],
        [2, 0, -1],
        [0, 1, 0],
        [1, -1, 0],
        [1, 1, -1]
    ])
    assert not p3.is_empty(include_boundary=True)
    assert not p3.is_empty(include_boundary=False)

    # Empty region
    p4 = Polytope([
        [0, 0, 1],  # y >= 0
        [-1, 0, -1] # y <= -1
    ])
    assert p4.is_empty(include_boundary=True)
    assert p4.is_empty(include_boundary=False)

    # Point 
    p5 = Polytope.rect((0, 0), (0, 0))
    assert not p5.is_empty(include_boundary=True)
    assert p5.is_empty(include_boundary=False)

    # Half-infinite rectangle 
    p6 = Polytope([
        [0, 1, 0],  # x >= 0
        [0, 0, 1], [-1, 0, 1] # 0 <= y <= 1
    ])
    assert not p6.is_empty(include_boundary=True)
    assert not p6.is_empty(include_boundary=False)

    # A halfplane in R^3
    p7 = Polytope([[1, 1, 1, 1]])
    assert not p7.is_empty(include_boundary=True)
    assert not p7.is_empty(include_boundary=False)

    # A polygon embedded in 3 dimensions
    p8 = Polytope([
        [0, 1, 0, 0], [1, -1, 0, 0], # 0 <= x <= 1
        [0, 0, 1, 0], [1, 0, -1, 0], # 0 <= y <= 1
        [0, 0, 0, 1], [0, 0, 0, -1]  # z = 0
    ])
    assert not p8.is_empty(include_boundary=True)
    assert p8.is_empty(include_boundary=False)

def run_V_init_test():
    verts = [[0, 0], [0, 1], [1, 1], [1, 0]]
    p = Polytope.from_V_rep(verts)
    assert set(tuple(v) for v in verts) == set(tuple(v) for v in p.get_vertices())

def run_union_test():
    # '-2 + 6x - 5/3y >= 0', '6 - 10x + 7/6y >= 0', '-25/32 + x >= 0'
    p1 = Polytope([[-2, 6, -frac(5,3)], [6, -10, frac(7,6)], [-frac(25,32), 1, 0]])
    # '2 - 6x + 5/3y >= 0', '-25/32 + x >= 0', '15/2 - 21/2x + 1/2y >= 0', '12 - 12x - y >= 0'
    p2 = Polytope([[2, -6, frac(5,3)], [-frac(25,32), 1, 0], [frac(15,2), -frac(21,2), frac(1,2)], [12, -12, -1]])

    union = Polytope.try_union([p1, p2])
    assert union is not None
    
    p1 = Polytope([
        [-19, 24, 0],
        [12, -12, -1],
        [16, -24, 2]
        ])
    p2 = Polytope([
        [19, -24, 0],
        [12, -12, -1],
        [-frac(1,2), 1, 0],
        [3, 0, -1],
        [16, -24, 2],
        [-1, 0, 1]
        ])
    
    union = Polytope([
        [12, -12, -1],
        [-frac(1,2), 1, 0],
        [3, 0, -1],
        [16, -24, 2],
        [-1, 0, 1]
        ])
    assert Polytope.try_union([p1, p2]) == union
    
    # Another test case - with a vertex at infinity (this fails)
    p1 = Polytope([
        [18, -22, 0],
        [12, -12, -1],
        #[3, 0, -1],
        #[-1, 0, 1],
        #[-frac(1,2), 1, 0],
        #[1, -1, 0]
        ])
    p2 = Polytope([
        [13, -15, -frac(1,2)],
        [-7, 5, frac(3,2)],
        [12, -12, -1],
        [3, 0, -1]
        ])
    U = Polytope.try_union([p1, p2])
    
    p1 = Polytope([
        [-3, 0, 2],
        [4, -4, -frac(1,2)],
        [-4, 4, frac(2,3)],
        [3, 0, -1],
        [6, -10, frac(7,6)],
        [13, -14, -frac(4,3)]
        ])
    p2 = Polytope([
        [13, -14, -frac(4,3)],
        [-4, 4, frac(1,2)],
        [3, 0, -1]
        ])
    p3 = Polytope([
        [-3, 0, 2],
        [-4, 4, frac(2,3)],
        [3, 0, -1],
        [6, -10, frac(7,6)],
        [13, -14, -frac(4,3)]
        ])
    assert Polytope.try_union([p1, p2]) == p3

    # Test edge cases - union between line and square 
    square = Polytope.rect((0, 1), (0, 1))
    line = Polytope.rect((frac(1,2), frac(1,2)), (0, 1))
    assert Polytope.try_union([square, line]) == square

    # Test unsuccessful union
    long_line = Polytope.rect((frac(1,2), frac(1,2)), (-1, 2))
    assert Polytope.try_union([square, long_line]) is None

def run_union_3d_test():
    # Higher dimensional polytope test
    # ['1 - x >= 0', '8 - 14x + 5/2y - z >= 0', '8 - 4x - z >= 0', '3 - y >= 0', '-34/3 + 32/3x + z >= 0']
    p1 = Polytope([
        [1, -1, 0, 0], 
        [8, -14, frac(5,2), -1],
        [8, -4, 0, -1],
        [3, 0, -1, 0],
        [-frac(34,3), frac(32,3), 0, 1]
    ])

    # ['z >= 0', '1 - x >= 0', '-1 + 2x >= 0', '1000000 - z >= 0', '1000000 - y >= 0'] ->
    p2 = Polytope([
        [-1, 2, 0, 0],
        [1, -1, 0, 0],
        [0, 0, 1, 0],
        [100, 0, -1, 0],
        [0, 0, 0, 1],
        [100, 0, 0, -1]
    ])
    assert Polytope.try_union([p1, p2]) == p2


def run_union_linset_test():
    # Test union edge cases which involve the use of lin_set
    A = Polytope([
        [-1, 0, 1, 0], [frac(3,2), 0, -1, 0],
        [0, 0, 0, 1], [3, 0, 0, -1],
        [-frac(1,2), 1, 0, 0], [frac(1,2), -1, 0, 0] # the linear constraint
    ], canonicalize=True)

    assert Polytope.try_union([A, A]) == A

def run_3_way_union_test():
    p1 = Polytope([[-8, 11, -1], [1, -1, 0], [-2, 0, 1], [frac(31, 12), 0, -1]])
    p2 = Polytope(
        [
            [8, -11, 1],
            [-7, 8, 0],
            [-2, 0, 1],
            [frac(31, 12), 0, -1],
            [-frac(103, 128), 1, -frac(1, 32)],
        ]
    )
    p3 = Polytope([[-2, 0, 1], [7, -8, 0], [-frac(103, 128), 1, -frac(1, 32)]])

    f = Piecewise(
        [Affine2([0, 0, 0], p1), Affine2([0, 0, 0], p2), Affine2([0, 0, 0], p3)]
    )
    f.simplify()

    assert len(f.pieces) == 1
    assert f.pieces[0].domain == Polytope(
        [
            [-frac(103, 128), 1, -frac(1, 32)],
            [1, -1, 0],
            [-2, 0, 1],
            [frac(31, 12), 0, -1],
        ]
    )

def run_setminus_test():
    A = Polytope.rect((0, 5), (0, 5))
    B = Polytope.rect((1, 2), (1, 2))
    
    A_minus_B = A.set_minus(B)
    
    assert not any(p.contains((1.5, 1.5)) for p in A_minus_B)
    assert any(p.contains((3, 3)) for p in A_minus_B)
    
def run_subs_test():
    A = Polytope.rect((1, 2), (1, 2), (1, 2))
    assert A.substitute({1: frac(3,2)}) == Polytope.rect((1, 2), (1, 2))

    B = Polytope([
            [0, 1, 0], # x >= 0
            [0, 0, 1], # y >= 0
            [1, -1, -1] # x + y <= 1
        ])
    assert B.substitute({1: frac(1,2)}) == Polytope.rect((0, frac(1,2)))

def run_lift_test():
    A = Polytope.rect((1, 2))
    assert A.lift([0, (3, 4)]) == Polytope.rect((1, 2), (3, 4))
    assert A.lift([(3, 4), 0]) == Polytope.rect((3, 4), (1, 2))
    
    B = Polytope.rect((-1, 1), (4, 9))
    assert B.lift([(2, 3), 0, (3, 6), 1]) == Polytope.rect((2, 3), (-1, 1), (3, 6), (4, 9))
    assert B.lift([0, 1]) == B
    
    C = Polytope([
            [0, 1, 0], # x >= 0
            [0, 0, 1], # y >= 0
            [1, -1, -1] # x + y <= 1
        ])
    assert C.lift([0, (5, 10), 1]) == Polytope([
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [1, -1, 0, -1],
            [-5, 0, 1, 0],
            [10, 0, -1, 0]
        ])

def run_scale_test():
    A = Polytope.rect((1, 2), (3, 4))
    assert A.scale_all([frac(2), frac(3)]) == Polytope.rect((2, 4), (9, 12))

    B = Polytope.from_V_rep([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert B.scale_all([frac(2), frac(3), frac(4)]) == \
        Polytope.from_V_rep([[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]])


run_polytope_tests()
run_emptiness_tests()
run_V_init_test()
run_union_test()
run_union_3d_test()
run_union_linset_test()
run_3_way_union_test()
run_setminus_test()
run_subs_test()
run_lift_test()
run_scale_test()