"""Murnaghan--Nakayama rule for irreducible characters of the symmetric group.

Rim hooks are removed using beta-sets (first column hook lengths): removing
a rim hook of length k corresponds to replacing an element b of the beta-set
by b - k, provided b - k >= 0 is not already in the set. The sign is
(-1) ** (number of beta-set elements strictly between b - k and b).
"""

from functools import lru_cache

from .partitions import Partition


def _beta_set(lam: Partition):
    length = len(lam)
    return [lam[i] + (length - 1 - i) for i in range(length)]


def _from_beta_set(beta) -> Partition:
    beta = sorted(beta, reverse=True)
    length = len(beta)
    parts = (beta[i] - (length - 1 - i) for i in range(length))
    return tuple(p for p in parts if p > 0)


@lru_cache(maxsize=None)
def sn_character(lam: Partition, rho: Partition) -> int:
    """Value of the irreducible character chi^lam of S_n on cycle type ``rho``.

    Both arguments are partitions of the same integer n. The cycles of ``rho``
    are removed in the given order (any order gives the same value).
    """
    if not rho:
        return 1 if not lam else 0
    k, rest = rho[0], rho[1:]
    beta = _beta_set(lam)
    members = set(beta)
    total = 0
    for b in beta:
        c = b - k
        if c >= 0 and c not in members:
            sign = -1 if sum(1 for x in beta if c < x < b) % 2 else 1
            new_lam = _from_beta_set([x for x in beta if x != b] + [c])
            total += sign * sn_character(new_lam, rest)
    return total


def degree(lam: Partition) -> int:
    """Degree chi^lam(1)."""
    return sn_character(lam, (1,) * sum(lam))
