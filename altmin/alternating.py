"""Irreducible characters of the alternating group A_n and their minima.

Irreducible characters of A_n are labelled by partitions lam of n up to
conjugation (James--Kerber, Section 2.5):

* if lam != lam', the restriction of chi^lam to A_n is irreducible;
* if lam == lam', the restriction splits as chi^{lam+} + chi^{lam-}. These
  agree with chi^lam / 2 except on the two A_n-classes of cycle type
  h = (h_11, h_22, ...), the diagonal hook lengths of lam, where

      chi^{lam +/-} = (eps +/- sqrt(eps * h_11 * h_22 * ...)) / 2,
      eps = (-1) ** ((n - d) / 2),   d = len(h).

  The values are non-real exactly when eps = -1.

For a character chi we define m(chi) as the smallest *real* value of chi on
A_n. For a split pair the set of real values is the same for chi^{lam+} and
chi^{lam-}, so one record covers both.
"""

import math
from dataclasses import dataclass
from typing import List, Optional

from .mn import degree, sn_character
from .partitions import (Partition, conjugate, diagonal_hooks,
                         is_even_cycle_type, partitions,
                         splits_in_alternating_group)


@dataclass(frozen=True)
class MinimumRecord:
    label: Partition          # representative with label >= its conjugate
    split: bool               # True if lam is self-conjugate
    minimum: float            # smallest real value on A_n
    argmin: Partition         # a cycle type where the minimum is attained
    has_nonreal_values: bool  # only possible for split characters


def irreducible_labels(n: int) -> List[Partition]:
    """One partition from each conjugate pair {lam, lam'} of n."""
    labels = []
    for lam in partitions(n):
        if lam >= conjugate(lam):
            labels.append(lam)
    return labels


def even_cycle_types(n: int) -> List[Partition]:
    return [rho for rho in partitions(n) if is_even_cycle_type(rho)]


def _split_values(lam: Partition, rho: Partition):
    """Values (chi^{lam+}(g), chi^{lam-}(g)) at an element g of type rho.

    Returns None for the non-real case.
    """
    n = sum(lam)
    h = diagonal_hooks(lam)
    if tuple(sorted(rho, reverse=True)) != h:
        v = sn_character(lam, rho) / 2
        return (v, v)
    eps = (-1) ** ((n - len(h)) // 2)
    prod = math.prod(h)
    if eps * prod < 0:
        return None
    root = math.sqrt(eps * prod)
    return ((eps + root) / 2, (eps - root) / 2)


def minimum_record(lam: Partition, classes: Optional[List[Partition]] = None) -> MinimumRecord:
    """Smallest real value of the A_n-character labelled by ``lam``."""
    n = sum(lam)
    if classes is None:
        classes = even_cycle_types(n)
    if lam != conjugate(lam):
        best = min((sn_character(lam, rho), rho) for rho in classes)
        return MinimumRecord(lam, False, best[0], best[1], False)
    best_val, best_rho, nonreal = None, None, False
    for rho in classes:
        vals = _split_values(lam, rho)
        if vals is None:
            nonreal = True
            continue
        v = min(vals)
        if best_val is None or v < best_val:
            best_val, best_rho = v, rho
    return MinimumRecord(lam, True, best_val, best_rho, nonreal)


def all_minima(n: int) -> List[MinimumRecord]:
    """Minimum records for all irreducible characters of A_n."""
    classes = even_cycle_types(n)
    return [minimum_record(lam, classes) for lam in irreducible_labels(n)]


def sum_of_squared_degrees(n: int) -> int:
    """Sum of chi(1)^2 over Irr(A_n); equals n!/2 (sanity check)."""
    total = 0
    for lam in irreducible_labels(n):
        d = degree(lam)
        total += 2 * (d // 2) ** 2 if lam == conjugate(lam) else d * d
    return total


def number_of_irreducibles(n: int) -> int:
    return sum(2 if lam == conjugate(lam) else 1 for lam in irreducible_labels(n))


def number_of_classes(n: int) -> int:
    return sum(2 if splits_in_alternating_group(rho) else 1 for rho in even_cycle_types(n))


def column_orthogonality_defect(n: int) -> float:
    """Largest relative defect in the column orthogonality relations of A_n.

    For every conjugacy class C of A_n the relation
    sum_{chi in Irr(A_n)} |chi(C)|^2 = |C_{A_n}(g)|, g in C, is checked using
    the values described in the module docstring (including the non-real
    values of the split characters). The centraliser order is z_rho / 2 if
    the S_n-class of type rho does not split, and z_rho if it splits.
    Returns max |lhs - rhs| / rhs over all classes; it should be ~1e-15.
    """
    labels = irreducible_labels(n)
    worst = 0.0
    for rho in even_cycle_types(n):
        z = 1
        for part in set(rho):
            mult = rho.count(part)
            z *= part ** mult * math.factorial(mult)
        splits = splits_in_alternating_group(rho)
        target = z if splits else z // 2
        total = 0.0
        for lam in labels:
            if lam != conjugate(lam):
                total += sn_character(lam, rho) ** 2
                continue
            h = diagonal_hooks(lam)
            if tuple(rho) != h:
                total += 2 * (sn_character(lam, rho) / 2) ** 2
                continue
            eps = (-1) ** ((n - len(h)) // 2)
            prod = math.prod(h)
            # |(eps +/- sqrt(eps*prod))/2|^2 summed over the pair = (1 + prod)/2
            if eps > 0:
                root = math.sqrt(prod)
                total += ((eps + root) / 2) ** 2 + ((eps - root) / 2) ** 2
            else:
                total += 2 * (1 + prod) / 4
        worst = max(worst, abs(total - target) / target)
    return worst
