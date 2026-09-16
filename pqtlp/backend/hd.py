"""Product-isogeny evaluation through the bundled theta arithmetic."""

from sage.all import pari

from . import params
from .theta.theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
from .theta.theta_structures.couple_point import CouplePoint


def Dim2Iso(K, n):
    """
    Compute the 2D (2^n,2^n)-isogeny with kernel K. Wrapper for the theta isogeny
    code. Does not handle split first steps.
    Input:
    - K: two pair of points on E1 x E2, the kernel of the isogeny
    - n: log 2 of the (polarized) degree (or the length of the chain)
    Output:
    - phi: the corresponding isogeny
    """

    (P1, P2), (Q1, Q2) = K
    P = CouplePoint(P1, P2)
    Q = CouplePoint(Q1, Q2)
    kernel = (P, Q)

    return EllipticProductIsogenySqrt(kernel, n)


def EvalTwoPoints(phi, X1, X2):
    """
    Given two 2D-points X1 and X2 on E1xE2 and phi: E1xE2 -> E3xE4, evaluate X1
    and X2 consistently fixing the sign.
    Input:
    - phi: a 2D isogeny
    - X1, X2: CouplePoints on the domain of phi
    Output:
    - phi(X1), phi(X2)
    """
    X3 = X1 - X2

    Y1 = phi(X1)
    Y2 = phi(X2)
    Y3 = phi(X3)

    images = []
    for P, Q, difference in zip(Y1, Y2, Y3):
        signs = {
            P - Q: (P, Q),
            P + Q: (P, -Q),
            -P - Q: (-P, Q),
            -P + Q: (-P, -Q),
        }
        images.append(signs[difference])
    (P1, Q1), (P2, Q2) = images
    return CouplePoint(P1, P2), CouplePoint(Q1, Q2)


def EmbeddedIsogeny(phi, n, d, basis=None):
    """
    Given the (2^n, 2^n)-isogeny phi embedding a 1D isogeny psi of degree d,
    recover the embedded isogeny and evaluate it on a provided basis (P, Q).
    Input:
    - phi: the 2D-isogeny
    - n: (polarized) 2-degree of phi
    - d: degree of psi
    - basis = (P, Q): basis of the 2^n torsion on the starting curve; if not
      provided, the standard basis (P0, Q0) on E0 is used
    Output:
    - E_psi: the codomain of psi
    - P_psi, Q_psi = psi(P), psi(Q)
    """
    if not basis:
        P = params.P0
        Q = params.Q0
    else:
        P, Q = basis

    T = CouplePoint(P, phi.E2(0))
    S = CouplePoint(Q, phi.E2(0))

    # Evaluate 2D iso on the points
    phi_T, phi_S = EvalTwoPoints(phi, T, S)

    # Detect the correct side with pairings
    cof = (params.p**2 - 1) // 2**n
    tPQ = pari.elltatepairing(phi.E1, P, Q, 2**n) ** (cof * d)
    t1 = pari.elltatepairing(phi.codomain()[0], phi_T[0], phi_S[0], 2**n) ** cof

    if tPQ == t1:
        return phi.codomain()[0], phi_T[0], phi_S[0]

    assert (
        pari.elltatepairing(phi.codomain()[1], phi_T[1], phi_S[1], 2**n) ** cof == tPQ
    )
    return phi.codomain()[1], phi_T[1], phi_S[1]
