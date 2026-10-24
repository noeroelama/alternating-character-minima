"""Table 2: m(chi^{(n-r,mu)}) on A_n and the ratio m / n^(r/2) versus -c_mu.

Output: results/table_ratio.csv
"""

import argparse
import csv

from _common import label, result_path

from altmin import c_mu, minimum_fixed_level

MUS = [(1, 1), (1, 1, 1), (3, 1), (2, 1, 1), (1, 1, 1, 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, nargs="+", default=[40, 80, 160])
    args = parser.parse_args()

    rows = []
    for mu in MUS:
        r = sum(mu)
        row = [label(mu), f"{-c_mu(mu):.4f}"]
        for n in args.n:
            value, x = minimum_fixed_level(mu, n)
            row += [value, f"{value / n ** (r / 2):.4f}", "X=" + label(x)]
        rows.append(row)
        print(row, flush=True)
    header = ["mu", "-c_mu"]
    for n in args.n:
        header += [f"m(n={n})", f"ratio(n={n})", f"argmin(n={n})"]
    with open(result_path("table_ratio.csv"), "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
