"""Character polynomials and exact minima of fixed-level characters.

For a partition mu of r and n >= r + mu_1, the value of chi^{(n-r, mu)} at a
permutation g depends only on X_1(g), ..., X_r(g), the numbers of cycles of
length 1, ..., r. It is given by the character polynomial q_mu. We store q_mu
in the binomial basis

    q_mu = sum_t  c_t * prod_i binom(X_i, t_i),

which makes the coefficients integers and the interpolation triangular.
"""

from functools import lru_cache
from math import comb
from typing import Dict, List, Optional, Tuple

from .mn import sn_character
from .partitions import Partition

Exponents = Tuple[int, ...]


def weighted_tuples(r: int) -> List[Exponents]:
    """All (t_1, ..., t_r) with sum_i i * t_i <= r."""
    out: List[Exponents] = []

    def rec(i: int, remaining: int, current: List[int]) -> None:
        if i > r:
            out.append(tuple(current))
            return
        for t in range(remaining // i + 1):
            rec(i + 1, remaining - i * t, current + [t])

    rec(1, r, [])
    return out


def _cycle_type(x: Exponents, n: int) -> Partition:
    parts = [i for i, xi in enumerate(x, start=1) for _ in range(xi)]
    rest = n - sum(parts)
    if rest:
        parts.append(rest)
    return tuple(sorted(parts, reverse=True))


def _binomial_monomial(t: Exponents, x: Exponents) -> int:
    value = 1
    for ti, xi in zip(t, x):
        if ti:
            value *= comb(xi, ti)
    return value


@lru_cache(maxsize=None)
def character_polynomial(mu: Partition) -> Dict[Exponents, int]:
    """Coefficients c_t of q_mu in the binomial basis (zero terms omitted).

    Interpolation points: for every exponent vector x of weight <= r, evaluate
    chi^{(N - r, mu)} at the cycle type with X_i = x_i (i <= r) plus one long
    cycle of length > r, where N = 4r + 2. Since binom(x_i, t_i) vanishes
    unless t_i <= x_i, the system is unitriangular.
    """
    r = sum(mu)
    big_n = 4 * r + 2
    lam = (big_n - r,) + tuple(mu)
    points = sorted(weighted_tuples(r), key=sum)
    coeff: Dict[Exponents, int] = {}
    for x in points:
        value = sn_character(lam, _cycle_type(x, big_n))
        known = sum(c * _binomial_monomial(t, x) for t, c in coeff.items())
        coeff[x] = value - known
    return {t: c for t, c in coeff.items() if c}


def evaluate(q: Dict[Exponents, int], x: Exponents) -> int:
    return sum(c * _binomial_monomial(t, x) for t, c in q.items())


@lru_cache(maxsize=None)
def _rest_feasible(points: int, min_part: int, even_parts_parity: int) -> bool:
    """Can ``points`` be split into cycles of length >= min_part whose number
    of even-length cycles has the given parity?"""
    if points == 0:
        return even_parts_parity == 0
    for length in range(min_part, points + 1):
        flip = 1 if length % 2 == 0 else 0
        if _rest_feasible(points - length, min_part, (even_parts_parity - flip) % 2):
            return True
    return False


def minimum_fixed_level(mu: Partition, n: int, group: str = "A"
                        ) -> Tuple[int, Exponents]:
    """Exact minimum of chi^{(n - r, mu)} over A_n (group="A") or S_n ("S").

    Requires n >= r + mu_1. Returns (minimum, X) where X = (X_1, ..., X_r)
    describes a minimising element (cycles of length > r are not recorded).
    """
    r = sum(mu)
    if n < r + mu[0]:
        raise ValueError("character polynomial not valid for this n")
    q = character_polynomial(tuple(mu))
    best: Optional[int] = None
    arg: Exponents = ()

    def rec(i: int, remaining: int, x: List[int]) -> None:
        nonlocal best, arg
        if i > r:
            if group == "A":
                even_small = sum(x[j - 1] for j in range(2, r + 1, 2)) % 2
                ok = _rest_feasible(remaining, r + 1, even_small)
            else:
                ok = remaining == 0 or remaining >= r + 1
            if ok:
                value = evaluate(q, tuple(x))
                if best is None or value < best:
                    best, arg = value, tuple(x)
            return
        for xi in range(remaining // i + 1):
            rec(i + 1, remaining - i * xi, x + [xi])

    rec(1, n, [])
    assert best is not None
    return best, arg
