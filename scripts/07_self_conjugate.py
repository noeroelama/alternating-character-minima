"""Lemma 3.2 and Proposition 3.3 (split characters), Table S2.

For every self-conjugate partition lambda of n, 9 <= n <= N (default 40),
check that the smaller of chi^lambda(3-cycle) and chi^lambda(2^2) is <= -3,
and that the degree exceeds (n+1)(n-2) for n >= 10. Also list the
self-conjugate partitions with 10 <= n <= 18 and their degrees.

Output: results/self_conjugate.txt
"""

import argparse

from _common import label, result_path

from altmin.mn import degree, sn_character
from altmin.partitions import conjugate, partitions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmax", type=int, default=40)
    args = parser.parse_args()

    lines = ["Self-conjugate partitions with 10 <= n <= 18 (Table S2):",
             "n | lambda | degree | (n+1)(n-2)"]
    all_ok = True
    for n in range(9, args.nmax + 1):
        count = 0
        for lam in partitions(n):
            if tuple(lam) != tuple(conjugate(lam)):
                continue
            count += 1
            f = degree(lam)
            v3 = sn_character(lam, (3,) + (1,) * (n - 3))
            v22 = sn_character(lam, (2, 2) + (1,) * (n - 4))
            ok = min(v3, v22) <= -3 and (n < 10 or f > (n + 1) * (n - 2))
            all_ok &= ok
            if 10 <= n <= 18:
                lines.append(f"{n} | {label(lam)} | {f} | {(n + 1) * (n - 2)}")
            if not ok:
                lines.append(f"FAILED: {label(lam)}")
        print(f"n={n}: {count} self-conjugate partitions checked", flush=True)
    lines.append(f"Proposition 3.3 checked for 9 <= n <= {args.nmax}")
    lines.append("ALL OK" if all_ok else "SOME CHECKS FAILED")
    with open(result_path("self_conjugate.txt"), "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(lines[-1])


if __name__ == "__main__":
    main()
