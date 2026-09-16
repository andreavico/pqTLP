from sage.all import EllipticCurve


def theta_null_point_to_montgomery_curve(O0):
    """
    Given a level 2 theta null point (a:b), compute a Montgomery curve equation.
    We use the model where the 4-torsion point (1:0) above (a:-b) is sent to
    (1:1) in Montgomery coordinates.

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    a, b = O0

    aa = a**2
    bb = b**2

    T1 = aa + bb
    T2 = aa - bb

    # Montgomery coefficient
    A = -(T1**2 + T2**2) / (T1 * T2)

    # Construct curve
    F = a.parent()
    E = EllipticCurve(F, [0, A, 0, 1, 0])
    return E


def theta_point_to_montgomery_point(O0, O):
    """
    Given an elliptic curve in Montgomery form
       E : y^2 = x^3 + Ax^2 + x
    and the theta null point with coordinates
       O0 = (a, b)

    Converts a theta point O = (U : V) on θ_(a,b)
    to a point P = (X : Z) on E

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    a, b = O0
    U, V = O
    X = a * V + b * U
    Z = a * V - b * U

    return (X, Z)
