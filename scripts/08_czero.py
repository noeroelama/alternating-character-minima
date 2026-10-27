"""Proposition 1.2(ii), Remark 7.1 and Theorem D (the case c_mu = 0).

1. For all mu with 2 <= |mu| <= R (default 12): if all parts of mu are
   even then c_mu = 0 (Proposition 1.2(ii)); for |mu| even and mu with
   empty 2-core, c_mu > 0 iff mu has an odd part (Remark 7.1).
2. Theorem D(a): the minimum of chi^{(n-4,2,2)} over A_n lies between
   -n^(4/3)/4 - 4n - 1 and the value of the explicit element of the proof.
3. Theorem D(b): the explicit elements of the proof give the stated values
   for (5), (3,2), (2,2,1) and (4,1).

Output: results/czero.txt
"""

import argparse
import math

from _common import label, result_path

from altmin import minimum_fixed_level
from altmin.charpoly import character_polynomial, evaluate
from altmin.hfunction import c_mu
from altmin.partitions import is_even_cycle_type, partitions


def empty_two_core(mu) -> bool:
    beta = [mu[i] + len(mu) - 1 - i for i in range(len(mu))]
    even = sum(1 for b in beta if b % 2 == 0)
    odd = len(beta) - even
    return even == odd if len(beta) % 2 == 0 else even == odd + 1


def counts(cycles, r):
    x = [0] * r
    for c in cycles:
        if c <= r:
            x[c - 1] += 1
    return tuple(x)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rmax", type=int, default=12)
    parser.add_argument("--nmax", type=int, default=150)
    args = parser.parse_args()
    lines, ok = [], True

    for r in range(2, args.rmax + 1):
        bad = []
        for mu in partitions(r):
            positive = c_mu(mu) > 1e-12
            if all(p % 2 == 0 for p in mu) and positive:
                bad.append(label(mu))
            if r % 2 == 0 and empty_two_core(mu) and positive != any(p % 2 for p in mu):
                bad.append(label(mu))
        ok &= not bad
        lines.append(f"r={r}: Hermite criterion " + ("OK" if not bad else f"FAILED {bad}"))

    q22 = character_polynomial((2, 2))
    for n in (216, 266, 316):
        a = math.ceil(n ** (1 / 3)) + 1
        ext = next(e for e in (0, 5, 7) if (n - a - e) % 3 == 0)
        l = (n - a - ext) // 3
        cycles = [1] * a + [3] * l + ([ext] if ext else [])
        assert is_even_cycle_type(tuple(sorted(cycles, reverse=True)))
        value = evaluate(q22, counts(cycles, 4))
        m, _ = minimum_fixed_level((2, 2), n)
        lower = -n ** (4 / 3) / 4 - 4 * n - 1
        good = lower <= m <= value
        ok &= good
        lines.append(f"(n-4,2,2), n={n}: min={m}, element={value}, "
                     f"lower bound -n^(4/3)/4-4n-1={lower:.1f} " + ("OK" if good else "FAILED"))

    for n in range(30, args.nmax + 1):
        if n % 4 == 0:
            k, gamma = n // 2, []
        elif n % 2 == 1:
            k = (n - 7) // 2 if ((n - 7) // 2) % 2 == 0 else (n - 9) // 2
            gamma = [n - 2 * k]
        else:
            k, gamma = (n - 8) // 2, [8]
        expected = {(5,): -k * (k - 1) // 2, (3, 2): -k * (k - 1) // 2,
                    (2, 2, 1): -k * (k - 3) // 2}
        assert is_even_cycle_type(tuple(sorted([2] * k + gamma, reverse=True)))
        for mu, exp in expected.items():
            v = evaluate(character_polynomial(mu), counts([2] * k + gamma, 5))
            ok &= v == exp
        k4 = 2 * (n // 8)
        ext = {0: 0, 1: 7, 2: 11}[(n - 2 * k4) % 3]
        l4 = (n - 2 * k4 - ext) // 3
        cycles4 = [2] * k4 + [3] * l4 + ([ext] if ext else [])
        assert is_even_cycle_type(tuple(sorted(cycles4, reverse=True)))
        v = evaluate(character_polynomial((4, 1)), counts(cycles4, 5))
        ok &= v == -(k4 - 1) * l4
    lines.append(f"Theorem D(b) elements checked for 30 <= n <= {args.nmax}")
    lines.append("ALL OK" if ok else "SOME CHECKS FAILED")
    for line in lines:
        print(line)
    with open(result_path("czero.txt"), "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
