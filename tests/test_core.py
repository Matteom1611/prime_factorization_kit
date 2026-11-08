import unittest

from prime_factorization_kit import factorize, is_prime, pollard_rho


class TestIsPrime(unittest.TestCase):
    def test_small_primes(self):
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 97, 997):
            self.assertTrue(is_prime(p), p)

    def test_small_composites(self):
        for c in (0, 1, 4, 6, 8, 9, 15, 21, 25, 49, 100, 1001):
            self.assertFalse(is_prime(c), c)

    def test_negative_numbers_are_not_prime(self):
        for n in (-2, -3, -7, -100):
            self.assertFalse(is_prime(n), n)

    def test_large_primes(self):
        # Well-known primes.
        self.assertTrue(is_prime(2**61 - 1))  # Mersenne prime
        self.assertTrue(is_prime(10**18 + 9))

    def test_large_composites(self):
        # 10**18 + 9 is prime; 10**18 + 11 = 3 * 333333333333333337.
        self.assertFalse(is_prime(10**18 + 11))
        # Product of two large primes.
        p = 1000000007
        q = 1000000009
        self.assertFalse(is_prime(p * q))

    def test_carmichael_number(self):
        # 561 is the smallest Carmichael number; deterministic MR must
        # reject it.
        self.assertFalse(is_prime(561))


class TestPollardRho(unittest.TestCase):
    def test_returns_factor_of_composite(self):
        n = 1000000007 * 1000000009
        d = pollard_rho(n)
        self.assertGreater(d, 1)
        self.assertLess(d, n)
        self.assertEqual(n % d, 0)

    def test_even_number_returns_two(self):
        self.assertEqual(pollard_rho(4), 2)
        self.assertEqual(pollard_rho(100), 2)

    def test_deterministic_for_same_input(self):
        n = 1000000007 * 1000000009
        self.assertEqual(pollard_rho(n), pollard_rho(n))

    def test_rejects_prime(self):
        with self.assertRaises(ValueError):
            pollard_rho(101)

    def test_rejects_too_small(self):
        with self.assertRaises(ValueError):
            pollard_rho(3)


class TestFactorize(unittest.TestCase):
    def test_zero(self):
        self.assertEqual(factorize(0), [0])

    def test_one(self):
        self.assertEqual(factorize(1), [])

    def test_prime(self):
        self.assertEqual(factorize(97), [97])

    def test_powers_of_two(self):
        self.assertEqual(factorize(1024), [2] * 10)

    def test_typical_composite(self):
        self.assertEqual(factorize(360), [2, 2, 2, 3, 3, 5])

    def test_negative_input(self):
        self.assertEqual(factorize(-12), [-1, 2, 2, 3])

    def test_product_matches_input(self):
        for n in (2, 3, 4, 12, 100, 997, 1024, 123456, 999983):
            factors = factorize(n)
            product = 1
            for f in factors:
                product *= f
            self.assertEqual(product, n, n)

    def test_large_semiprime(self):
        p = 1000000007
        q = 1000000009
        self.assertEqual(factorize(p * q), [p, q])

    def test_large_prime(self):
        n = 10**18 + 9
        self.assertEqual(factorize(n), [n])

    def test_sorted_output(self):
        factors = factorize(2 * 3 * 5 * 7 * 11 * 13)
        self.assertEqual(factors, sorted(factors))

    def test_repeated_calls_are_stable(self):
        n = 1000000007 * 1000000009
        self.assertEqual(factorize(n), factorize(n))


if __name__ == "__main__":
    unittest.main()
