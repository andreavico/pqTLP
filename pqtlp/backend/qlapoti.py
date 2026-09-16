"""Qlapoti ideal-to-isogeny conversion."""

from sage.all import QQ, ZZ, Matrix, gcd, inverse_mod, round, vector

from . import ec, hd, params, sos
from . import quaternions as qt


def qlapoti(J):
    """
    Given an ideal J, compute two equivalent ideals such that n(I1) + n(I2) =
    2^e
    Input:
    - J: the target left O0-ideal
    Output:
    - mu1, mu2: the two quaternions whose reduced norm sum to 2^e in the ideal
      I equivalent to J
    - N: the norm of I
    - I: an ideal equivalent to J
    - beta_ij: the element in J whose corresponding equivalent ideal is I
    """
    # Smallest equivalent ideal
    I, beta_ij = qt.ReducedIdeal(J)
    I_basis = qt.LLLReducedBasis(I)

    N = ZZ(I.norm())

    keep_alpha = False
    lam = 1
    alpha = alpha_0 = None

    while True:
        # Bound on the scaling on alpha
        if lam > 20000:
            keep_alpha = False

        # Scaling alpha
        if keep_alpha:
            alpha += alpha_0
            lam += 1
            while gcd(lam, N) > 1:
                # Need 2*alpha invertible
                alpha += alpha_0
                lam += 1

        # New generator
        else:
            lam = 1
            alpha_0 = qt.SmallishGenerator(I, N, I_basis=I_basis)
            alpha = alpha_0
            alpha_0_norm = alpha.reduced_norm()

        a_alpha, b_alpha, _, _ = alpha.coefficient_tuple()

        M = ZZ(2**params.f - 2 * lam**2 * alpha_0_norm / N)
        assert M > 0

        if keep_alpha:
            a_alpha_inv2 = inverse_mod(ZZ(2 * a_alpha), N)
            T = (M * a_alpha_inv2) % N
        else:
            a_alpha_inv2 = inverse_mod(ZZ(2 * a_alpha), N)
            x = ZZ((2 * b_alpha * a_alpha_inv2) % N)
            L = Matrix(ZZ, [[N - x, 1], [N, 0]])
            L = L.LLL()

            T = ZZ((M * a_alpha_inv2) % N)
            L_inv = Matrix(QQ, L).inverse()
            if qt._succ_min(L)[1] < 1 and gcd(2 * a_alpha, 2 * b_alpha) == 1:
                # Check for reasonable lattices
                keep_alpha = True

        v_target = vector(ZZ, [-T, 0])
        AB_vec = vector([round(c) for c in v_target * L_inv])
        v_close = AB_vec * L

        A, B = v_close - v_target

        assert (2 * (a_alpha * A + b_alpha * B)) % N == M % N

        M2 = M - 2 * a_alpha * A - 2 * b_alpha * B
        M2 = ZZ(M2 / N)

        # Complete the square
        M4 = 2 * M2 - A**2 - B**2
        if M4 < 0:
            continue

        # Unsolvable cases
        if M4 % 8 == 0:
            continue
        if A % 2 == B % 2 == 0:
            if M4 % 4 != 0:
                continue
        elif A % 2 == B % 2 == 1:
            if M4 % 4 != 2:
                continue
        else:
            if M4 % 4 != 1:
                continue

        ab = sos.sum_of_squares(M4)
        if not ab:
            continue

        # Back substitutions
        ad1, bd1 = ab
        if ad1 % 2 != A % 2:
            ad1, bd1 = bd1, ad1

        assert ad1 % 2 == A % 2 and bd1 % 2 == B % 2

        a1, b1 = ZZ((ad1 + A) / 2), ZZ((bd1 + B) / 2)
        a2, b2 = A - a1, B - b1

        gamma1 = params.B([a1, b1, 0, 0])
        gamma2 = params.B([a2, b2, 0, 0])

        mu1 = N * gamma1 + alpha
        mu2 = N * gamma2 + alpha

        theta = mu2 * mu1.conjugate() / N

        tx, _, tz, _ = 2 * theta
        if not ((tx % 2 != tz % 2) and ([c % 4 for c in 2 * theta].count(2) == 1)):
            keep_alpha = False
            continue

        break

    assert mu1 in I
    assert mu2 in I
    assert mu1.reduced_norm() + mu2.reduced_norm() == 2**params.f * N
    assert theta == mu2 * mu1.conjugate() / N

    return mu1, mu2, theta, N, I, beta_ij


def IdealToIsogeny(J):
    """
    Ideal to isogeny algorithm from qlapoti (2025/1604)
    Code adapted from github.com/KULeuven-COSIC/Qlapoti
    Input:
    - J: a left O0-ideal
    Output:
    - EJ: the codomain of phi_J
    - PJ, QJ: the images under phi_J of the standard basis (P0, Q0)
    """
    # Solve the norm equation
    mu1, _, theta, N, _, beta_ij = qlapoti(J)

    # Given the basis (P0, Q0) the kernel is given by
    # ((P0, d1^-1 * theta(P0), (Q0, d1^-1 * theta(Q0)))
    # -> compute d1^-1 * theta(P) using precomputations [Sec 3.2.1.1]
    d1 = mu1.reduced_norm() / N
    d_inv = inverse_mod(d1, 2**params.f)
    theta_P, theta_Q = ec.EvalBasisPrecomp(d_inv * theta)
    K = ((params.P0, theta_P), (params.Q0, theta_Q))

    # Compute the 2D isogeny and recover phi_I1
    phi = hd.Dim2Iso(K, params.f)
    EI, PI, QI = hd.EmbeddedIsogeny(phi, params.f, d1)

    # Recover phi_J
    theta_j = (mu1 * beta_ij) / N
    PJ, QJ = ec.EvalBasisPrecomp(d_inv * theta_j, (PI, QI))

    return EI, PJ, QJ
