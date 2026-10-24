import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from altmin import (all_minima, c_mu, character_polynomial, degree,
                    minimum_fixed_level, minimum_record, number_of_classes,
                    number_of_irreducibles, sn_character,
                    sum_of_squared_degrees)
from altmin.alternating import _split_values
from altmin.partitions import is_even_cycle_type, partitions


class TestSymmetricGroup(unittest.TestCase):
    def test_s4_character_table(self):
        # rows (4),(3,1),(2,2); columns 1^4, 2 1^2, 2^2, 3 1, 4
        cols = [(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)]
        expected = {
            (4,): [1, 1, 1, 1, 1],
            (3, 1): [3, 1, -1, 0, -1],
            (2, 2): [2, 0, 2, -1, 0],
        }
        for lam, row in expected.items():
            self.assertEqual([sn_character(lam, c) for c in cols], row)

    def test_degrees_hook_formula(self):
        self.assertEqual(degree((4, 4)), 14)
        self.assertEqual(degree((3, 2, 1)), 16)


class TestAlternatingGroup(unittest.TestCase):
    def test_a5(self):
        recs = {r.label: r for r in all_minima(5)}
        self.assertEqual(sorted(degree(l) // (2 if r.split else 1) for l, r in recs.items()),
                         [1, 3, 4, 5])
        plus, minus = _split_values((3, 1, 1), (5,))
        self.assertAlmostEqual(plus, (1 + math.sqrt(5)) / 2)
        self.assertAlmostEqual(minus, (1 - math.sqrt(5)) / 2)

    def test_orthogonality_counts(self):
        for n in range(5, 13):
            self.assertEqual(sum_of_squared_degrees(n), math.factorial(n) // 2)
            self.assertEqual(number_of_irreducibles(n), number_of_classes(n))

    def test_threshold_small(self):
        for n in range(9, 15):
            good = sorted(r.label for r in all_minima(n) if r.minimum >= -1)
            self.assertEqual(good, sorted([(n,), (n - 1, 1), (n - 2, 2)]))


class TestCharacterPolynomials(unittest.TestCase):
    def test_known_polynomials(self):
        self.assertEqual(character_polynomial((1, 1)),
                         {(0, 0): 1, (0, 1): -1, (1, 0): -1, (2, 0): 1})
        self.assertEqual(character_polynomial((2, 1)),
                         {(0, 0, 1): -1, (1, 0, 0): 1, (2, 0, 0): -2, (3, 0, 0): 2})

    def test_against_brute_force(self):
        for mu in [(1, 1), (2,), (3,), (2, 1), (1, 1, 1), (4,), (3, 1), (2, 2)]:
            r = sum(mu)
            for n in range(r + mu[0], 17):
                lam = (n - r,) + mu
                brute = min(sn_character(lam, rho) for rho in partitions(n)
                            if is_even_cycle_type(rho))
                self.assertEqual(minimum_fixed_level(mu, n)[0], brute, (mu, n))


class TestHFunction(unittest.TestCase):
    def test_c_mu(self):
        cases = {(2,): 0.0, (1, 1): 1 / 2, (1, 1, 1): 1 / 3, (3, 1): 1 / 8,
                 (2, 1, 1): 1 / 4, (1, 1, 1, 1): 1 / 4, (2, 2): 0.0}
        for mu, value in cases.items():
            self.assertAlmostEqual(c_mu(mu), value, places=9)


if __name__ == "__main__":
    unittest.main()
