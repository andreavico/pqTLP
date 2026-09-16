"""Torsion representations and the public CGL pushforward."""

import hashlib

from sage.all import ZZ, Matrix, Zmod, inverse_mod, pari, proof, vector

from .backend import ec, hd, misc, params
from .backend import qlapoti as qlpt
from .backend.theta.montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
from .backend.theta.montgomery_isogenies.kummer_line import KummerLine
from .parameters import PSI_TORSION_EXPONENT, WALK_SEED_BYTES, Q
from .types import Isogeny2D


def scale_full_basis_point_to(P, torsion_exponent: int):
    """Map a full-basis image down to 2^e torsion."""

    torsion_exponent = int(torsion_exponent)
    if not (1 <= torsion_exponent <= int(params.f)):
        raise ValueError(f"torsion_exponent must satisfy 1 <= e <= {params.f}")
    return (ZZ(2) ** ZZ(params.f - torsion_exponent)) * P


def ideal_to_2d_repr(
    I,
    *,
    base_curve,
    base_full_basis,
    base_change_matrix,
    base_ideal,
) -> Isogeny2D:
    """Compute the 2D representation of an ideal-generated isogeny.

    ``I`` is a left ideal of the secret randomized base order.
    Composing it with ``base_ideal`` lets the Qlapoti
    backend evaluate the isogeny from its fixed reference curve without
    exposing that setup-only ideal in the returned representation.
    """

    evaluation_ideal = base_ideal * I
    E, P_transport_image, Q_transport_image = qlpt.IdealToIsogeny(evaluation_ideal)
    P_full, Q_full = (
        c1 * P_transport_image + c2 * Q_transport_image for c1, c2 in base_change_matrix
    )
    P_image = scale_full_basis_point_to(P_full, PSI_TORSION_EXPONENT)
    Q_image = scale_full_basis_point_to(Q_full, PSI_TORSION_EXPONENT)
    return Isogeny2D(
        domain=base_curve,
        codomain=E,
        full_domain_basis=base_full_basis,
        images=(P_image, Q_image),
    )


def kernel_decomposed_to_ideal_2e(c1, c2, e: int):
    """Return the O0-ideal for ``<c1*P0 + c2*Q0>`` in ``E0[2^e]``.

    This is the decomposed-kernel-to-ideal map with the exponent made explicit
    instead of hard-wired to ``params.f``.
    """

    ker = kernel_generator_2e(c1, c2, e)
    modulus = ZZ(2) ** ZZ(e)
    I = params.O0 * ker + params.O0 * modulus
    if ZZ(I.norm()) != modulus:
        raise RuntimeError(f"kernel ideal has wrong norm {I.norm()} != {modulus}")
    return I


def kernel_generator_2e(c1, c2, e: int):
    """Return the quaternion generator for the decomposed kernel ideal."""

    e = int(e)
    if not (1 <= e <= int(params.f)):
        raise ValueError(f"e must satisfy 1 <= e <= {params.f}")

    modulus = ZZ(2) ** ZZ(e)
    v = vector(Zmod(modulus), (c1, c2))

    M_theta = -params.mat_i + 2 * params.mat_ij2 + params.mat_1k2
    M_theta = M_theta.change_ring(Zmod(modulus))
    mat_i = params.mat_i.change_ring(Zmod(modulus))

    d1, d2 = M_theta.transpose() * v
    M = Matrix(Zmod(modulus), 2, [c1, d1, c2, d2])

    a, b = M**-1 * mat_i.transpose() * v
    a, b = ZZ(a), ZZ(b)

    return params.B((a + b / 2, -1, b, b / 2))


def deterministic_basis_2e(E, e: int):
    """Return the public, deterministic basis of ``E[2^e]``.

    Setup and solving must use the same basis; errors must never trigger a
    random replacement. The reference curve uses its precomputed basis.
    """

    e = int(e)
    if not 1 <= e <= params.f:
        raise ValueError(f"e must satisfy 1 <= e <= {params.f}")
    if E == params.E0:
        scale = ZZ(2) ** (params.f - e)
        return scale * params.P0, scale * params.Q0
    return ec.TorsionBasis(E, e=ZZ(e))


def point_coords_in_basis_2e(X, P, Q, e: int):
    """Return ``a, b`` such that ``X = a P + b Q`` in ``E[2^e]``."""

    e = int(e)
    D = ZZ(2) ** ZZ(e)
    E = P.curve()
    zeta = pari.ellweilpairing(E, P, Q, D)
    a = ZZ(pari.fflog(pari.ellweilpairing(E, X, Q, D), zeta, D))
    b = ZZ(pari.fflog(pari.ellweilpairing(E, X, -P, D), zeta, D))
    if X != a * P + b * Q:
        raise RuntimeError("failed to recover coordinates in 2^e basis")
    return a, b


def evaluate_2d_on_torsion(sigma: Isogeny2D, X):
    """Evaluate ``sigma`` on a point in the represented domain ``2^a`` torsion."""

    P, Q = sigma.torsion_basis
    c1, c2 = point_coords_in_basis_2e(X, P, Q, PSI_TORSION_EXPONENT)
    return c1 * sigma.P_image + c2 * sigma.Q_image


def normalized_torsion_images(sigma: Isogeny2D):
    """Return ``q^-1 * psi(P), q^-1 * psi(Q)`` on the ``2^a`` basis."""

    inv_q = inverse_mod(Q, ZZ(2) ** PSI_TORSION_EXPONENT)
    return inv_q * sigma.P_image, inv_q * sigma.Q_image


def evaluate_basis_from_2d(
    E_start,
    E_sigma,
    B_start_full,
    B_start,
    B_sigma,
):
    """Evaluate a 2D representation on a full ``2^f`` domain basis."""

    a = PSI_TORSION_EXPONENT
    q = Q
    P_start_full, Q_start_full = B_start_full
    P_start, Q_start = B_start
    P_sigma, Q_sigma = B_sigma

    Phi = hd.Dim2Iso(
        ((P_start, P_sigma), (Q_start, Q_sigma)),
        a,
    )

    phi1_P, mpsi2_P = Phi(hd.CouplePoint(P_start_full, E_sigma(0)))
    phi1_Q, mpsi2_Q = Phi(hd.CouplePoint(Q_start_full, E_sigma(0)))

    P_sigma_full, Q_sigma_full = deterministic_basis_2e(E_sigma, int(params.f))
    hphi2_P, psi1_P = Phi(hd.CouplePoint(E_start(0), P_sigma_full))
    hphi2_Q, psi1_Q = Phi(hd.CouplePoint(E_start(0), Q_sigma_full))

    pair_start = pari.ellweilpairing(
        E_start,
        P_start_full,
        Q_start_full,
        ZZ(2) ** ZZ(params.f),
    )
    pair_phi1 = pari.ellweilpairing(
        phi1_P.curve(),
        phi1_P,
        phi1_Q,
        ZZ(2) ** ZZ(params.f),
    )
    pair_start_q = pair_start**q
    pair_start_q_inv = pair_start_q ** (-1)

    if pair_phi1 not in [pair_start_q, pair_start_q_inv]:
        phi1_P, phi1_Q = mpsi2_P, mpsi2_Q
        hphi2_P, hphi2_Q = psi1_P, psi1_Q
        pair_phi1 = pari.ellweilpairing(
            phi1_P.curve(),
            phi1_P,
            phi1_Q,
            ZZ(2) ** ZZ(params.f),
        )

    if pair_phi1 == pair_start_q_inv:
        phi1_Q = -phi1_Q
    elif pair_phi1 != pair_start_q:
        raise ValueError("wrong input degree q for 2D representation")

    pair_sigma = pari.ellweilpairing(
        E_sigma,
        P_sigma_full,
        Q_sigma_full,
        ZZ(2) ** ZZ(params.f),
    )
    pair_hphi2 = pari.ellweilpairing(
        hphi2_P.curve(),
        hphi2_P,
        hphi2_Q,
        ZZ(2) ** ZZ(params.f),
    )
    pair_sigma_cof = pair_sigma ** (ZZ(2) ** ZZ(a) - q)
    pair_sigma_cof_inv = pair_sigma_cof ** (-1)

    if pair_hphi2 == pair_sigma_cof_inv:
        hphi2_Q = -hphi2_Q
    elif pair_hphi2 != pair_sigma_cof:
        raise ValueError("wrong input degree 2^a-q for 2D representation")

    M_dlog = ec.ChangeOfBasis((hphi2_P, hphi2_Q), (phi1_P, phi1_Q))
    c11, c12, c21, c22 = M_dlog.list()
    scalar = ZZ(2) ** ZZ(a) - q
    sigma_P = scalar * (c11 * P_sigma_full + c12 * Q_sigma_full)
    sigma_Q = scalar * (c21 * P_sigma_full + c22 * Q_sigma_full)

    return sigma_P, sigma_Q


def evaluate_2d_on_full_basis(sigma: Isogeny2D):
    """Return ``sigma`` evaluated on its full ``2^f`` domain basis."""

    return evaluate_basis_from_2d(
        sigma.domain,
        sigma.codomain,
        sigma.full_domain_basis,
        sigma.torsion_basis,
        normalized_torsion_images(sigma),
    )


def kernel_from_basis_2e(basis, scalar: ZZ):
    """Return the CGL kernel generator ``P + scalar*Q`` for a ``2^e`` basis."""

    P, Q = basis
    return P + ZZ(scalar) * Q


def seeded_cgl_scalar(
    seed: bytes,
    E,
    chunk_index: int,
    e: int,
    forbidden_backtrack=None,
):
    """Derive one deterministic CGL chunk from the paper's public walk seed.

    The hash includes the current curve, chunk index, and previous dual
    order-two subgroup.  Rejection sampling excludes that subgroup at chunk
    boundaries, so concatenating the cyclic chunks remains non-backtracking.
    """

    seed = bytes(seed)
    chunk_index = int(chunk_index)
    e = int(e)
    if len(seed) != WALK_SEED_BYTES:
        raise ValueError(f"walk seed must be exactly {WALK_SEED_BYTES} bytes")
    if chunk_index < 0 or e <= 0:
        raise ValueError("invalid chunk index or exponent")

    basis = deterministic_basis_2e(E, e)
    forbidden_encoding = b""
    if forbidden_backtrack is not None:
        forbidden_encoding = misc.encode_fp2(forbidden_backtrack[0])

    scalar_bytes = (e + 7) // 8
    for counter in range(256):
        transcript = (
            b"pqtlp-hwalk-v1"
            + seed
            + chunk_index.to_bytes(8, "big")
            + counter.to_bytes(2, "big")
            + misc.encode_curve_j(E)
            + forbidden_encoding
        )
        scalar = ZZ(
            int.from_bytes(hashlib.shake_256(transcript).digest(scalar_bytes), "big")
        ) % (ZZ(2) ** ZZ(e))
        K = kernel_from_basis_2e(basis, scalar)
        if (
            forbidden_backtrack is None
            or (ZZ(2) ** ZZ(e - 1)) * K != forbidden_backtrack
        ):
            return scalar, basis
    raise RuntimeError("failed to derive a non-backtracking CGL chunk")


def cgl_dual_backtrack_subgroup(phi, basis, e: int):
    """Return the dual order-two subgroup forbidden for the next chunk."""

    e = int(e)
    Q_image = phi(basis[1])
    forbidden = (ZZ(2) ** ZZ(e - 1)) * Q_image
    if forbidden.is_zero():
        raise RuntimeError("failed to recover the dual backtracking subgroup")
    return forbidden


def check_kernel_order_2e(K, e: int, label: str = "kernel") -> None:
    """Check that ``K`` has exact order ``2^e`` and annotate it for Sage."""

    e = int(e)
    if e <= 0:
        raise ValueError("kernel exponent must be positive")
    order = ZZ(2) ** e
    if not (order * K).is_zero() or ((order // 2) * K).is_zero():
        raise RuntimeError(f"{label} does not have full order 2^{e}")
    K._order = order


def isogeny_from_kernel_2e_x_only(E, K, e: int):
    """Construct a degree-``2^e`` isogeny from ``K`` using Kummer arithmetic."""

    check_kernel_order_2e(K, e)
    L = KummerLine(E)
    phi = KummerLineIsogeny(L, L(K), ZZ(2) ** ZZ(e))
    E_out = phi.codomain().curve()
    E_out.set_order((params.p + 1) ** 2)
    return phi, E_out


def composite_dual_fast(phi):
    """Build a composite dual without re-proving the fixed characteristic."""

    # Sage constructs a composite dual by dualising every degree-two factor.
    # Its curve-model matching asks whether the (already fixed and verified)
    # base-field characteristic is prime for every factor.  A proved
    # primality test at each step dominates the public solver at full-size
    # parameters.  Probable-prime mode is sufficient for these redundant
    # internal checks and is scoped so callers retain Sage's normal policy.
    with proof.WithProof("arithmetic", False):
        return phi.dual()


def montgomery_isogeny_from_kernel(E, K, e: int):
    """Construct Sage's Montgomery-model isogeny from a full-order kernel."""

    check_kernel_order_2e(K, e, label="pushed chunk kernel")
    E.set_order((params.p + 1) ** 2)
    phi = E.isogeny(K, model="montgomery")
    phi.codomain().set_order((params.p + 1) ** 2)
    return phi


def push_2d_through_cgl_chunk(
    sigma: Isogeny2D,
    scalar: ZZ,
    chunk_size: int,
    *,
    basis=None,
):
    """Push a 2D-represented odd isogeny through one public CGL chunk.

    This implements one step of the Section 5 pushforward computation.  The current
    representation is ``sigma: E_i -> E_i'`` on ``2^a`` torsion.  The CGL chunk
    ``phi_i: E_i -> E_{i+1}`` has kernel ``<P + scalar Q>`` for the public
    deterministic ``2^u`` basis of ``E_i``.  The output represents
    ``[phi_i]_* sigma: E_{i+1} -> E'_{i+1}``.
    """

    E = sigma.domain
    a = PSI_TORSION_EXPONENT
    u = int(chunk_size)
    if u <= 0:
        raise ValueError("chunk_size must be positive")
    if u > a:
        raise NotImplementedError("current 2D evaluator requires chunk_size <= a")
    if u > int(params.f) - a:
        raise ValueError("chunk_size must satisfy u <= f - a")
    chunk_basis = deterministic_basis_2e(E, u) if basis is None else basis
    phi, _, K = cgl_chunk_isogeny(E, scalar, u, basis=chunk_basis)
    next_forbidden = cgl_dual_backtrack_subgroup(phi, chunk_basis, u)

    sigma_K = evaluate_2d_on_torsion(sigma, K)
    phi_prime = montgomery_isogeny_from_kernel(sigma.codomain, sigma_K, u)
    E_sigma_next = phi_prime.codomain()

    E_next = phi.codomain()
    P_next_full, Q_next_full = deterministic_basis_2e(E_next, int(params.f))
    P_old_full, Q_old_full = sigma.full_domain_basis
    sigma_P_full, sigma_Q_full = evaluate_2d_on_full_basis(sigma)

    dual_phi = composite_dual_fast(phi)
    # For u <= f-a, evaluate the dual on 2^(f-a-u) times the full basis.
    # The pushforward identity then yields the action on the desired 2^a
    # basis directly:
    #   phi' sigma dual_phi(2^(f-a-u) R)
    #       = sigma_push(2^(f-a) R).
    refresh_scale = ZZ(2) ** ZZ(int(params.f) - a - u)
    hphi_P = dual_phi(refresh_scale * P_next_full)
    hphi_Q = dual_phi(refresh_scale * Q_next_full)
    pP, pQ = point_coords_in_basis_2e(
        hphi_P,
        P_old_full,
        Q_old_full,
        int(params.f),
    )
    qP, qQ = point_coords_in_basis_2e(
        hphi_Q,
        P_old_full,
        Q_old_full,
        int(params.f),
    )

    push_phi_P = phi_prime(sigma_P_full)
    push_phi_Q = phi_prime(sigma_Q_full)
    # The scaled dual inputs give the action on [2^(f-a)]R directly.
    P_image = pP * push_phi_P + pQ * push_phi_Q
    Q_image = qP * push_phi_P + qQ * push_phi_Q

    out = Isogeny2D(
        domain=E_next,
        codomain=E_sigma_next,
        full_domain_basis=(P_next_full, Q_next_full),
        images=(P_image, Q_image),
    )
    return out, next_forbidden


def cgl_chunk_isogeny(E, scalar: ZZ, e: int, basis=None):
    """Compute the chunk isogeny with kernel ``<P + scalar*Q>``."""

    e = int(e)
    if basis is None:
        P, Q = deterministic_basis_2e(E, e)
    else:
        P, Q = basis
    K = kernel_from_basis_2e((P, Q), scalar)
    phi = montgomery_isogeny_from_kernel(E, K, e)
    return phi, (P, Q), K
