"""Randomized base ideals, short equivalents, and ideal-side pushforward."""

from typing import Optional

from sage.all import QQ, ZZ, Matrix, gcd, randint

from .backend import params
from .backend import qlapoti as qlpt
from .backend import quaternions as qt
from .parameters import SHORTCUT_SEARCH_ATTEMPTS, VERTICAL_DEGREE


def order_z_plus_nO(O, N: ZZ):
    """Return the quaternion order ``Z + N O``.

    Algorithm 2 samples a short element that is scalar modulo ``N``.  Representing
    ``Z + N O`` as an explicit order lets Sage compute the lattice intersection
    ``I_delay cap (Z + N O)`` used for that search.
    """

    B = O.quaternion_algebra()
    OB = O.basis()
    rows = []
    for b in OB:
        rows.append([QQ(c) * QQ(N) for c in b.coefficient_tuple()])
    rows.append([QQ(1), QQ(0), QQ(0), QQ(0)])
    M = Matrix(QQ, rows)
    d = M.denominator()
    H = (M * d).change_ring(ZZ).echelon_form()
    basis = []
    for row in H.rows()[:4]:
        coeffs = [QQ(c) / QQ(d) for c in row]
        basis.append(B(coeffs))
    return B.quaternion_order(basis)


def short_equivalent_ideal(
    I_delay,
    congruence_modulus: ZZ = VERTICAL_DEGREE,
    degree_coprime_to: Optional[ZZ] = None,
    attempts: int = SHORTCUT_SEARCH_ATTEMPTS,
):
    """Algorithm 2: find a short ideal equivalent to ``I_delay``.

    ``congruence_modulus`` is Algorithm 2's ``N``.  In the puzzle protocol this
    is the odd vertical degree ``q(2^a - q)``, not the ``2^a`` torsion modulus.
    """

    congruence_modulus = ZZ(congruence_modulus)
    if congruence_modulus <= 0:
        raise ValueError("congruence_modulus must be positive")
    if degree_coprime_to is None:
        degree_coprime_to = congruence_modulus
    degree_coprime_to = ZZ(degree_coprime_to)
    if degree_coprime_to <= 0:
        raise ValueError("degree_coprime_to must be positive")

    O = I_delay.left_order()
    scalar_order = order_z_plus_nO(O, congruence_modulus)
    lattice = I_delay.intersection(scalar_order)
    basis = qt.LLLReducedBasis(lattice)
    nI = ZZ(I_delay.norm())

    # Try the reduced basis vectors, then small random combinations.
    candidates = list(basis)

    for bound in [1, 2, 4, 8, 16, 32, 64]:
        for _ in range(attempts // 7):
            coeffs = [randint(-bound, bound) for _ in basis]
            alpha = sum(c * b for c, b in zip(coeffs, basis))
            if alpha:
                candidates.append(alpha)

    best = None
    best_norm = None
    for alpha in candidates:
        if alpha == 0:
            continue
        nr = ZZ(alpha.reduced_norm())
        if nr == 0:
            continue
        if nr % nI != 0:
            continue
        J = I_delay * alpha.conjugate() * (QQ(1) / QQ(nI))
        jn = ZZ(J.norm())
        # The replacement isogeny degree must be invertible where it is used.
        if gcd(jn, degree_coprime_to) != 1:
            continue
        if J.left_order() != I_delay.left_order():
            continue
        if not J.is_left_equivalent(I_delay):
            continue
        if best is None or jn < best_norm:
            best = J
            best_norm = jn

    if best is None:
        raise RuntimeError("failed to find a short equivalent ideal")
    return best


def sample_random_base_curve():
    """Sample the paper's public base curve and retain its secret ring witness.

    A random large odd-degree ideal from the fixed reference order plays the
    role of the setup's stored random walk.  Its codomain is the public curve;
    the ideal and its right order remain confined to ``SetupState``.
    """

    while True:
        I_base = qt.RandomIdealGivenNorm(params.D_mix, True)
        if not I_base:
            raise RuntimeError("failed to sample the randomized base ideal")
        E_base, P_full, Q_full = qlpt.IdealToIsogeny(I_base)
        E_base.set_order((params.p + 1) ** 2)
        if E_base.j_invariant() != params.E0.j_invariant():
            return I_base, E_base, (P_full, Q_full)


def sample_psi_ideal(*, base_ideal):
    """Sample an ideal of the fixed vertical degree ``N``."""

    I_reference = qt.RandomIdealGivenNorm(VERTICAL_DEGREE, False)
    if not I_reference:
        raise RuntimeError("failed to sample the vertical ideal")
    return qt.Pushforward(I_reference, base_ideal)


def ideal_side_pushforward_codomain(
    shortcut_ideal,
    psi_ideal,
    *,
    base_ideal,
    shortcut_push_lattice,
):
    """Compute the generation-side pushforward codomain using ideals."""

    # Setup precomputes the fixed right-order term of the pushforward.
    pushed_ideal = shortcut_ideal.conjugate() * psi_ideal + shortcut_push_lattice
    output_ideal = base_ideal * shortcut_ideal * pushed_ideal
    output_curve, _, _ = qlpt.IdealToIsogeny(output_ideal)
    return output_curve
