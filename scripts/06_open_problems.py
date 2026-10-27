"""Data for Section 7: growth of m(chi^{(n-r,mu)}) when c_mu = 0.

For each mu, report m at two values of n together with the empirical
exponent log(m(n2)/m(n1)) / log(n2/n1).

Output: results/open_problems.csv
"""

import argparse
import csv
import math

from _common import label, result_path

from altmin import c_mu, minimum_fixed_level

MUS = [(3,), (2, 1), (4,), (2, 2), (5,), (4, 1), (3, 2), (2, 2, 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n1", type=int, default=60)
    parser.add_argument("--n2", type=int, default=120)
    args = parser.parse_args()

    rows = []
    for mu in MUS:
        m1, x1 = minimum_fixed_level(mu, args.n1)
        m2, x2 = minimum_fixed_level(mu, args.n2)
        exponent = math.log(m2 / m1) / math.log(args.n2 / args.n1)
        rows.append([label(mu), f"{c_mu(mu):.4f}", m1, label(x1), m2, label(x2), f"{exponent:.3f}"])
        print(rows[-1], flush=True)
    with open(result_path("open_problems.csv"), "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["mu", "c_mu", f"m(n={args.n1})", "argmin X", f"m(n={args.n2})",
                         "argmin X", "empirical exponent"])
        writer.writerows(rows)


if __name__ == "__main__":
    main()
