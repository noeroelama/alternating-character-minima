"""Theorem A, Proposition 3.1 and Remark 3.2.

For 5 <= n <= N (default 28) list the irreducible characters of A_n whose
smallest real value is >= -1, check the expected answer for n >= 9, and run
the sanity checks sum chi(1)^2 = n!/2 and #Irr = #classes.

Output: results/threshold.txt
"""

import argparse
import math

from _common import label, result_path

from altmin import (all_minima, number_of_classes, number_of_irreducibles,
                    sum_of_squared_degrees)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmax", type=int, default=28)
    args = parser.parse_args()

    lines = ["n | labels with m(chi) >= -1 | split characters with non-real values | checks"]
    all_ok = True
    for n in range(5, args.nmax + 1):
        records = all_minima(n)
        good = sorted((r for r in records if r.minimum >= -1 - 1e-12),
                      key=lambda r: r.label, reverse=True)
        nonreal = sum(1 for r in records if r.has_nonreal_values)
        checks = (sum_of_squared_degrees(n) == math.factorial(n) // 2
                  and number_of_irreducibles(n) == number_of_classes(n))
        labels = ", ".join(label(r.label) + ("+-" if r.split else "") for r in good)
        if n >= 9:
            expected = sorted([(n,), (n - 1, 1), (n - 2, 2)], reverse=True)
            ok = [r.label for r in good] == expected
            all_ok &= ok and checks
            status = "OK" if ok and checks else "FAIL"
        else:
            status = "sanity OK" if checks else "sanity FAIL"
        lines.append(f"{n} | {labels} | {nonreal} | {status}")
        print(lines[-1], flush=True)
    lines.append("ALL OK" if all_ok else "SOME CHECKS FAILED")
    with open(result_path("threshold.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(lines[-1])


if __name__ == "__main__":
    main()
