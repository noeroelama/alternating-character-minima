"""The polynomial H_mu and the constant c_mu of the asymptotic theorem.

    H_mu(s) = sum_j chi^mu(1^{r-2j} 2^j) / ((r-2j)! j! 2^j) * s^{r-2j},

i.e. the Schur function s_mu specialised at p_1 = s, p_2 = 1, p_k = 0 (k>=3).
c_mu = -min(0, inf_{s >= 0} H_mu(s)).
"""

from fractions import Fraction
from math import factorial
from typing import List, Tuple

from .mn import sn_character
from .partitions import Partition


def h_coefficients(mu: Partition) -> List[Tuple[int, Fraction]]:
    """List of (exponent, coefficient) pairs of H_mu."""
    r = sum(mu)
    coeffs = []
    for j in range(r // 2 + 1):
        rho = tuple([2] * j + [1] * (r - 2 * j))
        value = sn_character(tuple(mu), rho)
        coeffs.append((r - 2 * j, Fraction(value, factorial(r - 2 * j) * factorial(j) * 2 ** j)))
    return coeffs


def h_value(mu: Partition, s: float) -> float:
    return float(sum(float(c) * s ** e for e, c in h_coefficients(mu)))


def _derivative_value(coeffs, s: float) -> float:
    return sum(float(c) * e * s ** (e - 1) for e, c in coeffs if e > 0)


def c_mu(mu: Partition, tol: float = 1e-12) -> float:
    """c_mu >= 0, computed from the critical points of H_mu on [0, S].

    The leading coefficient is chi^mu(1)/r! > 0, so H_mu is positive for
    large s; S is a bound (Fujiwara) beyond which H_mu is increasing.
    """
    coeffs = h_coefficients(mu)
    lead = float(coeffs[0][1])
    # Fujiwara-type bound for the real roots of H_mu': every root s satisfies
    # |s| <= 2 max_k |a_k / a_deg|^(1/(deg-k)), where a_k are the coefficients
    # of the derivative (Fujiwara's bound has an extra factor 1/2 for the
    # constant term; omitting it only enlarges the bound). The search grid
    # below is relative to this bound, so a needlessly large bound, such as a
    # Cauchy-type bound when chi^mu(1)/r! is small, makes the grid too coarse.
    deg = coeffs[0][0]
    if len(coeffs) > 1 and deg > 1:
        lead_d = lead * deg
        bound = 2.0 * max((abs(float(c)) * e / lead_d) ** (1.0 / (deg - e))
                          for e, c in coeffs[1:] if e > 0) if any(e > 0 for e, _ in coeffs[1:]) else 1.0
    else:
        bound = 1.0
    bound = max(bound, 2.0)
    grid = 20000
    candidates = [0.0, bound]
    prev_s, prev_d = 0.0, _derivative_value(coeffs, 0.0)
    for k in range(1, grid + 1):
        s = bound * k / grid
        d = _derivative_value(coeffs, s)
        if prev_d == 0.0:
            candidates.append(prev_s)
        elif prev_d * d < 0:
            lo, hi = prev_s, s
            for _ in range(200):
                if hi - lo <= tol * max(1.0, abs(hi)):
                    break
                mid = (lo + hi) / 2
                if _derivative_value(coeffs, lo) * _derivative_value(coeffs, mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            candidates.append((lo + hi) / 2)
        prev_s, prev_d = s, d
    minimum = min(h_value(mu, s) for s in candidates)
    return max(0.0, -minimum)
