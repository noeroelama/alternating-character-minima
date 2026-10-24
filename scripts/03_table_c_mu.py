"""Table 1: the constant c_mu of Theorem C for 2 <= |mu| <= 8.

Also records chi^mu on the involution class 1^eps 2^floor(r/2)
(Corollary 1.2) and the number of mu with c_mu > 0 for each r.

Output: results/table_c_mu.csv
"""

import csv

from _common import label, result_path

from altmin import c_mu, sn_character
from altmin.partitions import partitions


def main() -> None:
    rows = []
    summary = []
    for r in range(2, 9):
        positive = 0
        total = 0
        for mu in partitions(r):
            rho = tuple([2] * (r // 2) + [1] * (r % 2))
            c = c_mu(mu)
            rows.append([r, label(mu), f"{c:.6f}", sn_character(mu, rho)])
            positive += c > 1e-12
            total += 1
        summary.append(f"r={r}: {positive}/{total} partitions with c_mu > 0")
    with open(result_path("table_c_mu.csv"), "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["r", "mu", "c_mu", "chi_mu_on_1^eps_2^floor(r/2)"])
        writer.writerows(rows)
    with open(result_path("table_c_mu_summary.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
