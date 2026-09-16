"""Check parameter selection against exhaustive feasible choices."""

import unittest

from pqtlp.calibration import choose_parameters, estimate_times


class CalibrationTests(unittest.TestCase):
    def test_minimum_generation_cost(self):
        profile = {
            67: {"setup": 2, "tlp_pgen": 1, "base_solve": 0.5},
            134: {"setup": 6, "tlp_pgen": 1, "base_solve": 3},
            268: {"setup": 12, "tlp_pgen": 2, "base_solve": 8},
        }
        for target in (0.1, 0.5, 3, 8, 17, 60):
            for ratio in (0.9, 1, 1.2):
                with self.subTest(target=target, ratio=ratio):
                    candidates = []
                    for root in profile:
                        for count in range(1, 200):
                            generation, solve = estimate_times(root, count, profile)
                            if solve >= target * ratio:
                                candidates.append(
                                    (generation, abs(solve - target), root, count)
                                )
                    expected = min(candidates)
                    actual = choose_parameters(
                        target, profile=profile, min_solve_ratio=ratio
                    )
                    self.assertEqual(
                        (actual.root_delay, actual.tlp_count), expected[2:]
                    )
                    self.assertEqual(actual.predicted_pgen_seconds, expected[0])

    def test_invalid_durations(self):
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    choose_parameters(value)
                with self.assertRaises(ValueError):
                    choose_parameters(1, min_solve_ratio=value)


if __name__ == "__main__":
    unittest.main()
