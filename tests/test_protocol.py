import pickle
import unittest
from unittest.mock import patch

from sage.all import ZZ, next_prime, proof

import pqtlp as tlp
from pqtlp.backend import params
from pqtlp.crypto import unmask_next_puzzle
from pqtlp.curves import seeded_cgl_scalar


class ProtocolTests(unittest.TestCase):
    def test_lambda_128_parameters(self):
        self.assertEqual(tlp.SECURITY_BITS, 128)
        self.assertEqual(tlp.SOLUTION_BYTES, 16)
        self.assertEqual(tlp.P.nbits(), 326)
        self.assertTrue(tlp.P.is_prime())
        self.assertEqual((tlp.P + 1).valuation(2), 324)
        self.assertEqual((tlp.P + 1) // (2**324), 3)
        self.assertEqual(params.D_mix, next_prime(tlp.P**2))

        self.assertEqual(tlp.PSI_TORSION_EXPONENT, 257)
        self.assertEqual(tlp.CHUNK_SIZE, 67)
        self.assertEqual(tlp.Q, next_prime(ZZ(2) ** 256))
        self.assertEqual(
            tlp.VERTICAL_DEGREE,
            tlp.Q * (ZZ(2) ** 257 - tlp.Q),
        )
        self.assertEqual(tlp.VERTICAL_DEGREE.nbits(), 512)

    def test_randomized_base_flat_compiler_round_trip(self):
        solution = b"compiler-test".ljust(tlp.SOLUTION_BYTES, b"\0")
        puzzle = tlp.pgen(
            67,
            solution,
            tlp_count=3,
            workers=2,
        )

        self.assertNotEqual(
            puzzle.first_puzzle.psi.domain.j_invariant(),
            params.E0.j_invariant(),
        )
        self.assertEqual(len(puzzle.first_puzzle.delay_seed), 32)
        self.assertEqual(puzzle.first_puzzle.delay_length, 67)
        self.assertFalse(hasattr(puzzle.first_puzzle, "base_ideal"))
        self.assertEqual(len(puzzle.masked_puzzles), 2)

        current = puzzle.first_puzzle
        codomains = {current.psi.codomain.j_invariant()}
        for masked in puzzle.masked_puzzles:
            key = tlp.psolve_base(current)
            current = pickle.loads(unmask_next_puzzle(key, masked))
            self.assertIsInstance(current, tlp.BasePuzzle)
            codomains.add(current.psi.codomain.j_invariant())
        self.assertEqual(len(codomains), puzzle.tlp_count)
        self.assertEqual(tlp.psolve_base(current), solution)
        self.assertEqual(tlp.psolve(puzzle), solution)

    def test_two_chunk_seeded_delay_round_trip(self):
        solution = b"two-chunk-test".ljust(tlp.SOLUTION_BYTES, b"\0")
        state = tlp.setup_delay(
            134,
            walk_seed=b"S" * 32,
        )
        self.assertEqual(
            state.shortcut_ideal.left_order(),
            state.base_ideal.right_order(),
        )

        public_puzzle = tlp.pgen_base(state, solution)
        self.assertEqual(len(state.base_change_matrix), 2)
        proof_policy = proof.arithmetic()
        self.assertEqual(tlp.psolve_base(public_puzzle), solution)
        self.assertEqual(proof.arithmetic(), proof_policy)

        # Public walk derivation must not depend on either party's random state.
        from sage.all import set_random_seed

        with patch.object(state.base_curve, "random_point", side_effect=AssertionError):
            set_random_seed(1)
            first = seeded_cgl_scalar(
                state.delay_seed, state.base_curve, 0, tlp.CHUNK_SIZE
            )
            set_random_seed(2)
            second = seeded_cgl_scalar(
                state.delay_seed, state.base_curve, 0, tlp.CHUNK_SIZE
            )
        self.assertEqual(first, second)

    def test_invalid_workers_fail_before_setup(self):
        with patch("pqtlp.protocol.setup_delay") as setup:
            with self.assertRaises(ValueError):
                tlp.pgen(tlp.CHUNK_SIZE, bytes(tlp.SOLUTION_BYTES), workers=0)
            setup.assert_not_called()


if __name__ == "__main__":
    unittest.main()
