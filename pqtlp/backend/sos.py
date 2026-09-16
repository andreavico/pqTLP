"""Sum-of-two-squares search used by Qlapoti's norm equation."""

from sage.all import (
    ZZ,
    Factorization,
    gcd,
    is_pseudoprime,
    isqrt,
    prime_range,
    prod,
    two_squares,
)
from sage.rings.finite_rings.integer_mod import Mod

# Trial-division primes below 10,000. Generate the table once at import.
small_primes = list(prime_range(10000))
small_squares = {p: two_squares(p) for p in small_primes if p % 4 == 1}
good_prime_prod = prod(small_squares)
bad_prime_prod = prod(p for p in small_primes if p % 4 == 3)


def two_squares_factored(factors):
    """
    This is the function `two_squares` from sage, except we give it the
    factorisation of n already.
    """
    n = factors.expand()
    F = factors
    for p, e in F:
        if e % 2 == 1 and p % 4 == 3:
            raise ValueError(f"{n} is not a sum of two squares")

    if n == 0:
        return (0, 0)
    a = ZZ.one()
    b = ZZ.zero()
    for p, e in F:
        if p == 1:
            continue
        if e >= 2:
            m = p ** (e // 2)
            a *= m
            b *= m
        if e % 2 == 1:
            if p == 2:
                # (a + bi) *= (1 + I)
                a, b = a - b, a + b
            else:  # p = 1 mod 4
                if p in small_squares:
                    r, s = small_squares[p]
                else:
                    # Find a square root of -1 mod p.
                    # If y is a non-square, then y^((p-1)/4) is a square root of -1.
                    y = Mod(2, p)
                    while True:
                        s = y ** ((p - 1) / 4)
                        if not s * s + 1:
                            s = s.lift()
                            break
                        y += 1
                    # Apply Cornacchia's algorithm to write p as r^2 + s^2.
                    r = p
                    while s * s > p:
                        r, s = s, r % s
                    r %= s

                # Multiply (a + bI) by (r + sI)
                a, b = a * r - b * s, b * r + a * s

    a = a.abs()
    b = b.abs()
    assert a * a + b * b == n
    if a <= b:
        return (a, b)
    else:
        return (b, a)


def rep_gcd(a, b):
    """
    Given a and b returns (a1, common_factor) where
    a1*common_factor = a and common_factor contains the factors in common
    between a and b.
    """
    out = 1
    common_factor = gcd(a, b)
    while common_factor != 1:
        out *= common_factor
        a /= common_factor
        common_factor = gcd(a, common_factor)

    return a, out


def sum_of_squares_friendly(n):
    """
    We can write any n = x^2 + y^2 providing that there
    are no prime power factors p^k | n such that
    p = 3 mod 4 and k odd.
    """
    # We consider the odd part of n and try and determine if there are bad factors
    n_val = n.valuation(2)
    n_odd = n // (2**n_val)
    fact = [(2, n_val)]

    if n_odd % 4 == 3:
        return False, fact

    n_odd, bad_cof = rep_gcd(n_odd, bad_prime_prod)

    sbf = isqrt(bad_cof)
    if sbf**2 == bad_cof:
        fact.append((sbf, 2))
    else:
        return False, fact

    # Good primes 1 mod 4
    n_odd, good_cof = rep_gcd(n_odd, good_prime_prod)
    good_cof = ZZ(good_cof)

    if n_odd == 1:
        return True, Factorization([*fact, *good_cof.factor()])

    else:
        return is_pseudoprime(n_odd), Factorization(
            [*fact, *good_cof.factor(), (n_odd, 1)]
        )


def sum_of_squares(n):
    """
    Attempts to compute x,y such that n = x^2 + y^2
    """
    n = ZZ(n)

    if n == 0:
        return ZZ(0), ZZ(0)
    if n < 0:
        return []
    b, fact = sum_of_squares_friendly(n)
    if not b:
        return []
    return two_squares_factored(fact)
