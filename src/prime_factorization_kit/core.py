"""Core integer factorization routines.

Design notes:

* ``is_prime`` uses a deterministic Miller-Rabin with the well-known
  witness set {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}, which is
  proven correct for all n < 3.3 * 10**24. That comfortably covers the
  10**18 range this library targets, and avoids the cost of random
  witnesses.

* ``pollard_rho`` uses Brent's variant of Pollard's rho with the
  ``f(x) = x*x + c`` polynomial. Brent's cycle detection roughly halves
  the number of gcd calls compared to Floyd's and is the version used in
  production factoring code.

* ``factorize`` strips small primes by trial division first, then hands
  the remaining cofactor to Pollard's rho. Trial division by the primes
  up to 1000 removes the common small factors cheaply and leaves rho
  with a cofactor that is either prime or has no small factors, which is
  where rho performs best.

All arithmetic uses Python's arbitrary-precision integers, so the only
practical limit is the cost of rho on semiprimes with two ~10**9 factors.
"""

from __future__ import annotations

import math
import random
from typing import List

__all__ = ["factorize", "is_prime", "pollard_rho"]


# Deterministic Miller-Rabin witnesses for all n < 3.3 * 10**24.
_MR_WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

# Small primes used to strip cheap factors before invoking rho.
_SMALL_PRIMES = (
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
    223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379,
    383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
    463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563,
    569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
    647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
    743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829,
    839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937,
    941, 947, 953, 967, 971, 977, 983, 991, 997,
)


def _miller_rabin(n: int, witnesses) -> bool:
    """Return True if n passes Miller-Rabin for every witness.

    ``n`` must be odd and greater than 2. The caller is responsible for
    the small cases (0, 1, 2) and for even n.
    """
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for a in witnesses:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def is_prime(n: int) -> bool:
    """Return True if ``n`` is a prime number.

    Deterministic for all n < 3.3 * 10**24 via a fixed witness set, so
    no randomness is involved and repeated calls agree.

    Negative numbers, 0 and 1 are not prime.
    """
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    # Cheap trial division before the heavier Miller-Rabin loop.
    for p in _SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False
    return _miller_rabin(n, _MR_WITNESSES)


def _pollard_rho_brent(n: int, seed: int) -> int:
    """Brent's variant of Pollard's rho. Returns a non-trivial factor of n.

    ``seed`` controls the polynomial constant and the initial value. The
    caller provides a deterministic seed so the whole library remains
    reproducible.
    """
    if n % 2 == 0:
        return 2

    rng = random.Random(seed)
    y = rng.randrange(1, n)
    c = rng.randrange(1, n)
    m = 128

    g = 1
    r = 1
    q = 1

    while g == 1:
        x = y
        for _ in range(r):
            y = (y * y + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            upper = min(m, r - k)
            for _ in range(upper):
                y = (y * y + c) % n
                q = (q * abs(x - y)) % n
            g = math.gcd(q, n)
            k += m
        r *= 2

    if g == n:
        # Backtrack one step at a time to find the factor that the
        # batched gcd collapsed into n.
        g = 1
        while g == 1:
            ys = (ys * ys + c) % n
            g = math.gcd(abs(x - ys), n)

    return g


def pollard_rho(n: int) -> int:
    """Return a non-trivial factor of composite ``n``.

    Uses Brent's cycle detection with a deterministic seed derived from
    ``n`` itself, so the same input always yields the same factor. This
    keeps the function pure and makes the test suite reproducible.

    Raises ValueError if ``n`` is prime or less than 4, since those have
    no non-trivial factor.
    """
    if n < 4:
        raise ValueError("pollard_rho requires n >= 4")
    if is_prime(n):
        raise ValueError("pollard_rho requires a composite n")
    if n % 2 == 0:
        return 2

    # Deterministic seed: a simple mix of n into a 64-bit-ish value.
    seed = (n * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
    factor = _pollard_rho_brent(n, seed)
    if factor == n:
        # Extremely unlikely, but retry with a different seed rather than
        # returning the trivial factor.
        factor = _pollard_rho_brent(n, seed ^ 0xA5A5A5A5A5A5A5A5)
    if factor == n:
        raise RuntimeError("pollard_rho failed to find a non-trivial factor")
    return factor


def factorize(n: int) -> List[int]:
    """Return the prime factorization of ``n`` as a sorted list.

    The list contains each prime factor with multiplicity, so
    ``factorize(12) == [2, 2, 3]`` and the product of the list equals
    ``n`` for positive ``n``.

    ``factorize(0)`` returns ``[0]`` and ``factorize(1)`` returns ``[]``.
    Negative inputs are factorized as ``-1`` times the factorization of
    their absolute value, with ``-1`` first in the result.
    """
    if n == 0:
        return [0]
    if n == 1:
        return []
    if n < 0:
        result = [-1]
        n = -n
    else:
        result = []

    # Strip small primes first; this is much faster than rho for
    # numbers with small factors and leaves rho a smaller cofactor.
    for p in _SMALL_PRIMES:
        if n % p == 0:
            while n % p == 0:
                result.append(p)
                n //= p
        if n == 1:
            return result

    # Factor the remaining cofactor with rho, recursing on composites.
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            result.append(m)
            continue
        d = pollard_rho(m)
        stack.append(d)
        stack.append(m // d)

    result.sort()
    return result
