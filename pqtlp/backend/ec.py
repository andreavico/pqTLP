"""Torsion bases and quaternion actions on elliptic curves."""

from sage.all import ZZ, Zmod, matrix, pari

from . import params


def PointCoords(X, P, Q, zeta=None, e=None):
    """
    Given a point X in E[2^e] and a basis (P, Q) of E[2^f] compute the vector
    (a, b) such that X = a*P + b*Q.
    Input:
    - X: a point in E[2^e]
    - P, Q: a basis for E[2^f]
    - zeta: (optional) the scaled tate pairing between P and Q
    - e: (optional) the order of X; if None, params.f is used; if e !=
      params.f, the output (a, b) is such that
        (2**(params.f - e)) * (aP + bQ) = X
    Output:
    - a, b such that X = a*P + b*Q
    """
    if not e:
        e = params.f
    E = X.curve()

    cof = (params.p**2 - 1) // 2**e
    if not zeta:
        zeta = pari.elltatepairing(E, P, Q, 2**params.f) ** cof

    zeta1 = pari.elltatepairing(E, X, Q, 2**e) ** cof
    zeta2 = pari.elltatepairing(E, X, -P, 2**e) ** cof

    a = ZZ(pari.fflog(zeta1, zeta, 2**e))
    b = ZZ(pari.fflog(zeta2, zeta, 2**e))

    assert (2 ** (params.f - e)) * (a * P + b * Q) == X

    return a, b


def ChangeOfBasis(B1, B2, e=None):
    """
    Compute the matrix M so that M * B1 = B2 where B1 = (P1, P2)^T and B2 =
    (Q1, Q2)^T are two basis of the 2^f torsion [Algorithm 2.5].
    Input:
    - B1 = (P1, P2) a basis for E[2^f]
    - B2 = (Q1, Q2) a basis for E[2^e]
    - e: optional, the order of the points; if None params.f is used
    Output:
    - A matrix M =(x_i) so that Q1 = x1*P1 + x2*P2 and Q2 = x3*P1 + x4*P2
    """
    if not e:
        e = params.f

    cof = (params.p**2 - 1) // 2**e
    P1, P2 = B1
    E = P1.curve()
    zeta = pari.elltatepairing(E, P1, P2, 2**params.f) ** cof

    Q1, Q2 = B2

    x1, x2 = PointCoords(Q1, P1, P2, zeta, e)
    x3, x4 = PointCoords(Q2, P1, P2, zeta, e)

    return matrix(Zmod(2**e), 2, [x1, x2, x3, x4])


def EvalMatrix(M, basis=None):
    """
    Given a matrix M and two points (P, Q) evaluate M * (P, Q)^T
    Input:
    - M: a matrix
    - basis = (P, Q): two points; if None, (P0, Q0) is used
    Output:
    - M * (P, Q)^T
    """
    if basis:
        P, Q = basis
    else:
        P, Q = params.P0, params.Q0

    x1, x2, x3, x4 = M.list()
    return x1 * P + x2 * Q, x3 * P + x4 * Q


def EvalBasisPrecomp(theta, basis=None):
    """
    Given a quaternion theta in O0 and a basis (P, Q) of E[2^f] evaluate theta
    on the points using the precomputed action of a basis of O0 on the
    2-torsion.
    Input:
    - theta: an element of O0
    - basis: a basis (P, Q) of E[2^f]; if None, the standard (P0, Q0) is used
    Output:
    - theta(P), theta(Q)
    """
    # Express theta wrt <1, i, (i+j)/2, (1+k)/2>
    a, b, c, d = theta
    x, y, z, w = map(ZZ, (a - d, b - c, 2 * c, 2 * d))

    M = x * params.mat_1 + y * params.mat_i + z * params.mat_ij2 + w * params.mat_1k2
    return EvalMatrix(M, basis)


def TorsionBasis(E_A, e=None):
    """
    Given a montgomery curve E_A and an integer e <= f, computes a
    deterministic basis (R, S) of E_A[2^e] such that S is above (0, 0).
    Modified version of [Algorithm 2.1] working with (x,y)-points and no hint.
    Input:
    - A: montgomery coefficient
    - e: 2-order of the torsion; if not provided, f is used
    Output:
    - R, S: a basis of E_A[2^e] such that S is above (0, 0)
    """
    A = E_A.a_invariants()[1]
    assert E_A.a_invariants() == (0, A, 0, 1, 0), "not a Montgomery model"

    h = params.Fp(0)
    _i = params.Fp2_i
    if A.is_square():
        while True:
            h += 1
            if (1 + h**2).is_square():  # Norm of (1 + ih)
                continue
            xR = -A / (1 + _i * h)  # Now this is non-square
            assert not xR.is_square()
            if E_A.is_x_coord(xR):
                break
    else:
        while True:
            h += 1
            xR = h * A
            if E_A.is_x_coord(xR):
                break

    xRS = -xR - A
    R = E_A.lift_x(xR)
    RS = E_A.lift_x(xRS)
    S = R - RS

    # Clear cofactors
    if not e:
        e = params.f

    cof = (params.p + 1) // 2**e
    R *= cof
    S *= cof

    assert (2 ** (e - 1) * S).xy() == (0, 0)
    return R, S
