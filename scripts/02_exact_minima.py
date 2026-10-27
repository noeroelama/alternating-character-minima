"""Theorem C: exact minima on A_n for (n-2,1,1), (n-3,3), (n-3,2,1), (n-4,4).

For (n-4,4) the minimum is min(-5, f(n)) for n >= 10.

Compares the closed formulas with the exact minimum computed from the
character polynomial, for every n from the first admissible value up to
--nmax (default 160; 120 for (n-4,4), which is slower).

Output: results/exact_minima.txt
"""

import argparse
from fractions import Fraction as F

from _common import label, result_path

from altmin import minimum_fixed_level


def f_211(n):
    return {0: -F(n - 2, 2), 1: -F(n - 1, 2), 2: -F(n - 2, 2), 3: -F(n - 5, 2)}[n % 4]


def f_33(n):
    return {0: -F(n, 2), 1: -F(n - 5, 2), 2: -F(n - 4, 2), 3: -F(n - 5, 2)}[n % 4]


def f_321(n):
    return {0: -F(n, 3), 1: -F(n - 4, 3), 2: -F(n - 2, 3)}[n % 3]


def f_44(n):
    # For n >= 10 the minimum is min(-5, f(n)); for n >= 20 this equals f(n).
    return min(F(-5), {0: -F(n, 3), 1: -F(n - 7, 3), 2: -F(n - 5, 3)}[n % 3])


CASES = [
    ((1, 1), f_211, 6),
    ((3,), f_33, 6),
    ((2, 1), f_321, 7),
    ((4,), f_44, 10),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmax", type=int, default=160)
    parser.add_argument("--nmax4", type=int, default=120)
    args = parser.parse_args()

    lines = []
    all_ok = True
    for mu, formula, n_from in CASES:
        r = sum(mu)
        n_max = args.nmax4 if mu == (4,) else args.nmax
        failures = []
        for n in range(r + mu[0], n_max + 1):
            value, _ = minimum_fixed_level(mu, n)
            if n >= n_from and value != formula(n):
                failures.append(n)
            if n < n_from and value != formula(n):
                lines.append(f"  note: {label((n - r,) + mu)} at n={n}: m={value}, formula gives {formula(n)}")
        ok = not failures
        all_ok &= ok
        lines.append(f"(n-{r},{','.join(map(str, mu))}): formula checked for {n_from} <= n <= {n_max}: "
                     + ("OK" if ok else f"FAILED at {failures}"))
        print(lines[-1], flush=True)
    lines.append("ALL OK" if all_ok else "SOME CHECKS FAILED")
    with open(result_path("exact_minima.txt"), "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(lines[-1])


if __name__ == "__main__":
    main()
