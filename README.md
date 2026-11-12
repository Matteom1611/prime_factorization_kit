# prime_factorization_kit

Factor integers up to around 10**18 into their prime factors using trial division, Pollard's rho, and deterministic Miller-Rabin primality testing.

## Usage

```python
from prime_factorization_kit import factorize, is_prime, pollard_rho

factorize(360)
# [2, 2, 2, 3, 3, 5]

factorize(1000000007 * 1000000009)
# [1000000007, 1000000009]

is_prime(10**18 + 9)
# True

pollard_rho(1000000007 * 1000000009)
# 1000000007  (or 1000000009; either non-trivial factor is valid)
```

`factorize(n)` returns a sorted list of the prime factors of `n` with
multiplicity, so the product of the list equals `n`. It handles the
awkward cases explicitly:

* `factorize(0)` returns `[0]`.
* `factorize(1)` returns `[]`.
* `factorize(-12)` returns `[-1, 2, 2, 3]`; the sign is factored out as a
  leading `-1` and the magnitude is factorized normally.

`is_prime(n)` is deterministic for all `n < 3.3 * 10**24` (the witness
set `{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}` is proven sufficient
in that range), so it does not depend on randomness and repeated calls
always agree.

`pollard_rho(n)` returns a single non-trivial factor of a composite `n`.
It raises `ValueError` if `n` is prime or less than 4, because those
inputs have no non-trivial factor. The seed is derived from `n` itself,
so the same input always produces the same factor.

## Why this exists

Factoring 64-bit integers is a common need in Project Euler-style
problems, key-generation sanity checks, and hash-table sizing. Python's
standard library has no integer factorization, and pulling in SymPy or
`primefac` adds a large dependency for what is a few hundred lines of
well-understood number theory.

The trade-off is speed versus dependency count. Trial division alone is
hopeless past 10**12, and a full quadratic sieve or general number field
sieve is far too much code for a small library. Pollard's rho with
Brent's cycle detection is the sweet spot: it factors a semiprime with
two ~10**9 factors in well under a second on modern hardware, uses only
`math` and `random` from the standard library, and stays short enough to
read in one sitting. The cost is that rho is not competitive on
adversarial inputs (a semiprime with two ~10**15 factors will take a
long time), so this library is not a replacement for a real factoring
tool on hard instances.

## The awkward edge

`pollard_rho` is not defined for primes. A prime has no non-trivial
factor, so the function raises `ValueError` rather than returning the
input or `1`. Callers who have an arbitrary integer in hand should check
`is_prime` first, or simply call `factorize`, which handles the prime
case internally by returning a single-element list.

## Design notes

The window stores values eagerly rather than keeping running aggregates. Running
sums drift with floating point over long streams, and recomputing from a small
buffer is cheap enough that the drift is not worth the speed.

