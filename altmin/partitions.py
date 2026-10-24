"""Integer partitions and cycle types.

A partition is a weakly decreasing tuple of positive integers. A cycle type
of a permutation of {1, ..., n} is a partition of n.
"""

from typing import Iterator, Optional, Tuple

Partition = Tuple[int, ...]


def partitions(n: int, max_part: Optional[int] = None) -> Iterator[Partition]:
    """Yield all partitions of ``n`` in reverse lexicographic order."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, max_part), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def conjugate(lam: Partition) -> Partition:
    """Return the conjugate (transposed) partition."""
    if not lam:
        return ()
    return tuple(sum(1 for part in lam if part > i) for i in range(lam[0]))


def is_even_cycle_type(rho: Partition) -> bool:
    """True iff permutations of cycle type ``rho`` lie in the alternating group.

    A permutation is even iff it has an even number of cycles of even length.
    """
    return sum(1 for part in rho if part % 2 == 0) % 2 == 0


def diagonal_hooks(lam: Partition) -> Partition:
    """Diagonal hook lengths h_11 > h_22 > ... of ``lam``."""
    lam_c = conjugate(lam)
    d = sum(1 for i, part in enumerate(lam) if part > i)
    return tuple(lam[i] + lam_c[i] - 2 * i - 1 for i in range(d))


def splits_in_alternating_group(rho: Partition) -> bool:
    """True iff the S_n-class of cycle type ``rho`` splits into two A_n-classes.

    This happens iff all parts are odd and pairwise distinct.
    """
    return len(set(rho)) == len(rho) and all(part % 2 == 1 for part in rho)
