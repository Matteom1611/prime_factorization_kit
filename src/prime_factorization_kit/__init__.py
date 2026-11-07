"""Prime factorization utilities for integers up to around 10**18.

Exports:
    factorize: return a sorted list of prime factors (with repetition).
    is_prime: deterministic Miller-Rabin for 64-bit integers.
    pollard_rho: find a non-trivial factor of a composite number.
"""

from .core import factorize, is_prime, pollard_rho

__all__ = ["factorize", "is_prime", "pollard_rho"]
