"""PARI Weil pairing used by product-curve points."""

from sage.all import pari


def weil_pairing_pari(P, Q, D, check=False):
    """
    Wrapper around Pari's implementation of the Weil pairing
    Allows the check of whether P,Q are in E[D] to be optional
    """
    if check:
        nP, nQ = D * P, D * Q
        if not nP.is_zero() or not nQ.is_zero():
            raise ValueError("points must both be n-torsion")

    return pari.ellweilpairing(P.curve(), P, Q, D)
