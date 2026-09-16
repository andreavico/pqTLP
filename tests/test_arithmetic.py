"""Algebraic regression checks for the cleaned arithmetic helpers."""

import unittest

from sage.all import ZZ, matrix

from pqtlp.backend import params
from pqtlp.backend.sos import sum_of_squares, two_squares_factored
from pqtlp.backend.theta.theta_isogenies.isomorphism import Isomorphism
from pqtlp.backend.theta.theta_structures.dimension_two import ThetaStructure
from pqtlp.backend.theta.utilities.discrete_log import weil_pairing_pari
from pqtlp.curves import check_kernel_order_2e


class ArithmeticTests(unittest.TestCase):
    def test_sum_of_squares(self):
        for n in (0, 1, 2, 5, 25, 45, 65, 125, 10009, 2**20 * 5**4):
            with self.subTest(n=n):
                a, b = sum_of_squares(n)
                self.assertEqual(a * a + b * b, n)
        for n in (-1, 3, 7, 15):
            self.assertFalse(sum_of_squares(n))
        with self.assertRaises(ValueError):
            two_squares_factored(ZZ(3).factor())

    def test_exact_kernel_order(self):
        exponent = 67
        P = 2 ** (params.f - exponent) * params.P0
        check_kernel_order_2e(P, exponent)
        for invalid in (2 * P, params.P0, params.E0(0)):
            with self.subTest(invalid=invalid):
                with self.assertRaises(RuntimeError):
                    check_kernel_order_2e(invalid, exponent)

    def test_pairing_membership_check(self):
        degree = 2**67
        scale = 2 ** (params.f - 67)
        P, Q = scale * params.P0, scale * params.Q0
        self.assertEqual(
            weil_pairing_pari(P, Q, degree, check=True), weil_pairing_pari(P, Q, degree)
        )
        with self.assertRaises(ValueError):
            weil_pairing_pari(params.P0, Q, degree, check=True)

    def test_theta_coordinate_inverse(self):
        domain = ThetaStructure(tuple(params.Fp2(x) for x in (1, 2, 3, 5)))
        point = domain(tuple(params.Fp2(x) for x in (7, 11, 13, 17)))
        change = Isomorphism(
            matrix(params.Fp2, 4, [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0])
        )
        change._domain = domain
        change._codomain = ThetaStructure(change.apply_isomorphism(domain.null_point()))
        self.assertEqual(change.dual()(change(point)), point)
        self.assertEqual(point.double_iter(0), point)
        self.assertEqual(point.double_iter(2), point.double().double())


if __name__ == "__main__":
    unittest.main()
