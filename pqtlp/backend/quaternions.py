"""Quaternion ideal arithmetic used by setup and Qlapoti."""

from sage.all import (
    QQ,
    ZZ,
    Matrix,
    ceil,
    gcd,
    is_pseudoprime,
    isqrt,
    kronecker,
    mod,
    randint,
    round,
    sqrt,
    valuation,
)

from . import params


def ModularSQRT(n, m):
    """
    Modular square root [Alg 3.1]
    Input:
    - m: an odd prime
    - n: an integer which is a square mod m
    Output:
    - an integer x s.t. x^2 = n mod m
    """
    if n % m == 0:
        return ZZ(0)
    if m % 4 == 3:
        return ZZ(mod(n, m) ** ((m + 1) / 4))

    # Tonelli-Shanks for m = 1 mod 4
    # https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python
    s = valuation(m - 1, 2)
    q = (m - 1) // 2**s

    z = 2
    while kronecker(z, m) != -1:
        z += 1

    c = mod(z, m) ** q
    n = mod(n, m)
    r = n ** ((q + 1) // 2)
    t = n**q

    e = s
    t2 = 0
    while (t - 1) != 0:
        t2 = t * t
        for i in range(1, m):
            if t2 - 1 == 0:
                break
            t2 = t2 * t2
        b = c ** (1 << (e - i - 1))
        r = r * b
        c = b * b
        t = t * c
        e = i
    return ZZ(r)


def GeneralizedRepresentInteger(M, om, O):
    """
    Given an order O and a target norm M, compute an element gamma in O with
    norm M [Algorithm 3.12]
    Input:
    - M: odd integer > p
    - om: quaternion such that q = -om^2 is a positive integer and q = 1 mod 4
    - O: a special extremal maximal order containing Z[om] + jZ[om] as
      suborder, with j from the standard basis of B_p,inf
    Output:
    - gamma in O with nrd(gamma) = M, or False
    """
    q = ZZ(-(om**2))
    bound = ceil(4 * M / (params.p * sqrt(q)))
    counter = 0

    m = isqrt(4 * M / params.p - q)
    while counter < bound:
        counter += 1
        z = randint(1, m)
        m1 = isqrt((4 * M - params.p * z**2) / (q * params.p))
        t = randint(-m1, m1)

        M1 = ZZ(4 * M - params.p * (z**2 + q * t**2))
        # ``M1`` is a disposable internal candidate.  A proved primality test
        # here dominated repeated puzzle generation at full-size parameters;
        # a strong probable-prime test is enough before Cornacchia, whose
        # returned representation and the resulting ideal norm are checked
        # exactly below.
        if is_pseudoprime(M1):
            res = Cornacchia(q, M1)
            if not res:
                continue
            x, y = res

            j = params.B.gens()[1]
            gamma = x + om * y + j * z + om * j * t
            d = 1
            while True:
                if gamma / (d + 1) in O:
                    d += 1
                else:
                    break

            if d != 2:
                continue
            return gamma / d

    return False


def Cornacchia(q, m):
    """
    Solve the equation x^2 + q*y^2 = m [Algorithm 3.11]
    Input:
    - q: positive integer at most m
    - m: an odd prime, or 2
    Output:
    - x, y such that x^2 + qy^2 = m or False if a solution do not exist
    """
    if not kronecker(-q, m):
        return False
    if m == 2:
        if q == 1:
            return 1, 1
        return False

    r = ModularSQRT(-q % m, m)
    s = m

    while r**2 > m:
        r, s = s % r, r

    x, y = r, isqrt((m - r**2) / q)
    if x**2 + q * y**2 == m:
        return x, y
    return False


def RandomIdealGivenNorm(N, prime):
    """
    Compute a random ideal from given norm [Alg 3.10]. It may fail if N is not
    prime.
    Input:
    - N: positive integer coprime with p, norm of the output ideal
    - prime: a boolean indicating whether N is prime
    Output:
    - a random left ideal J of O0 of norm N, or False if not found
    """
    if prime:
        while True:
            g1, g2, g3 = [randint(0, N - 1) for _ in range(3)]
            gamma = params.B([0, g1, g2, g3])
            n_gamma = gamma.reduced_norm()
            if kronecker(-n_gamma, N) == 1:
                gamma += ModularSQRT(-n_gamma, N)
                break
    else:
        gamma = GeneralizedRepresentInteger(
            params.QUAT_prime_cofactor * N, params.B.gens()[0], params.O0
        )
        if not gamma:
            return False

    # Scaling by a random beta mod N
    while True:
        x, y, z, w = [randint(0, N - 1) for _ in range(4)]
        beta = params.B([x, y, z, w])
        if gcd(beta.reduced_norm(), N) == 1:
            break
    ideal = params.O0 * (gamma * beta) + params.O0 * N
    if ZZ(ideal.norm()) != ZZ(N):
        raise RuntimeError("sampled ideal has the wrong norm")
    return ideal


def LLLReducedBasis(I):
    """
    Given an ideal I, computes an LLL-reduced basis of I. Calls pari on the
    Gram matrix of I.
    Input:
    - I: an ideal
    Output:
    - 4 elements of I forming an LLL-reduced basis
    """
    B = I.basis()
    M = []
    for a in B:
        M.append([QQ(2) * (a * b.conjugate()).reduced_trace() for b in B])
    G = Matrix(QQ, M)
    U = G.LLL_gram().transpose()
    return [sum(c * beta for c, beta in zip(row, B)) for row in U]


def ReducedIdeal(I):
    """
    Given an ideal I compute the ideal equivalent to I with the smallest norm
    Input:
    - I: a left O0-ideal
    Output:
    - J: the ideal equivalent to I with smallest norm
    - beta: the element in I generating J
    """
    B = LLLReducedBasis(I)
    beta = B[0]
    J = I * (beta.conjugate() / I.norm())
    assert J.is_left_equivalent(I)
    return J, beta


def SmallishGenerator(I, N, I_basis):
    """
    Compute a rather small random generator of I
    Input:
    - I: target left O0-ideal
    - N: norm of I
    - I_basis: an LLL-reduced basis of I
    Output:
    - alpha: a generator of I
    """
    while True:
        alpha = sum(randint(1, 100000) * gen for gen in I_basis)

        a_alpha = alpha.coefficient_tuple()[0]
        if gcd(2 * a_alpha, N) == 1 and gcd(alpha.reduced_norm(), N * N) == N:
            break

    return alpha


def _succ_min(L):
    """
    Compute the first two minima of L
    Input:
    - L: a matrix
    Output:
    - lam1, lam2: the two minima of L
    """
    fourth_root_p = round(params.p ** (1 / 4), 10)
    lam1 = round(L.row(0).norm() / fourth_root_p, 10)
    lam2 = round(L.row(1).norm() / fourth_root_p, 10)
    return lam1, lam2


def Pushforward(I, J):
    """
    Compute the pushforward of I under J, i.e. [J]_* I.
    Input:
    - I, J: two left O-ideals with coprime norm
    Output:
    - the ideal [J]_*I
    """
    assert I.left_order() == J.left_order()
    assert gcd(I.norm(), J.norm()) == 1

    O = J.right_order()
    return J.conjugate() * I + O * I.norm()
