"""Supplementary Section S5, cross-check: character-polynomial minimisation versus complete tables.

For every mu with |mu| <= 4 and every admissible n <= --nmax (default 24),
the minimum of chi^{(n-r,mu)} over A_n computed from the character
polynomial must equal the minimum over all even cycle types computed with
the Murnaghan--Nakayama rule.

Output: results/validation.txt
"""

import argparse

from _common import label, result_path

from altmin import minimum_fixed_level, sn_character
from altmin.partitions import is_even_cycle_type, partitions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmax", type=int, default=24)
    args = parser.parse_args()

    checked = 0
    failures = []
    for r in range(1, 5):
        for mu in partitions(r):
            for n in range(r + mu[0], args.nmax + 1):
                lam = (n - r,) + mu
                brute = min(sn_character(lam, rho) for rho in partitions(n)
                            if is_even_cycle_type(rho))
                fast, _ = minimum_fixed_level(mu, n)
                checked += 1
                if brute != fast:
                    failures.append((label(mu), n, brute, fast))
    msg = f"{checked} (mu, n) pairs checked, {len(failures)} mismatches"
    with open(result_path("validation.txt"), "w", newline="\n") as fh:
        fh.write(msg + "\n")
        for f in failures:
            fh.write(f"{f}\n")
    print(msg)


if __name__ == "__main__":
    main()
