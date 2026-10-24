# alternating-character-minima

Code and data accompanying the manuscript

> A. Nurwahid, *Minimum values of irreducible characters of the alternating groups* (in preparation).

For an irreducible character χ of the alternating group A<sub>n</sub>, let m(χ) be the smallest real value of χ. The scripts in this repository reproduce every table and every computer verification in the manuscript.

## Requirements

* Python 3.9 or newer
* No third-party packages (standard library only)

## Reproducing the results

```bash
python run_all.py
```

This runs the unit tests and then the scripts below in order. The total running time is about 3–4 minutes on a laptop. Output files are written to `results/`; the copies committed in this repository were produced with the same command.

| Script | Manuscript | Output |
|---|---|---|
| `scripts/01_threshold.py` | Theorem A, Proposition 3.1, Remark 3.2 | `results/threshold.txt` |
| `scripts/02_theorem_b.py` | Theorem B (exact minima, checked for n ≤ 160, resp. n ≤ 120) | `results/theorem_b.txt` |
| `scripts/03_table_c_mu.py` | Table 1 (constants c<sub>μ</sub>, \|μ\| ≤ 8) | `results/table_c_mu.csv`, `results/table_c_mu_summary.txt` |
| `scripts/04_asymptotics.py` | Table 2 (convergence in Theorem C) | `results/table_ratio.csv` |
| `scripts/05_validate.py` | Section 6.1 (cross-check of the two algorithms) | `results/validation.txt` |
| `scripts/06_open_problems.py` | Section 6.3 (cases with c<sub>μ</sub> = 0) | `results/open_problems.csv` |

Each script accepts `--help`; the ranges of n can be changed from the command line.

## Package overview

* `altmin/partitions.py` – partitions, conjugation, cycle types in A<sub>n</sub>
* `altmin/mn.py` – Murnaghan–Nakayama rule on β-sets for characters of S<sub>n</sub>
* `altmin/alternating.py` – irreducible characters of A<sub>n</sub>, including the split characters of self-conjugate partitions, and their minima
* `altmin/charpoly.py` – character polynomials q<sub>μ</sub> and exact minimisation of χ<sup>(n−r,μ)</sup> over A<sub>n</sub> or S<sub>n</sub>
* `altmin/hfunction.py` – the polynomial H<sub>μ</sub> and the constant c<sub>μ</sub>

Two independent methods are used. For complete character tables (n ≤ 28) every character is evaluated on every even cycle type. For characters of fixed level r the value depends only on the numbers X<sub>1</sub>, …, X<sub>r</sub> of short cycles. The minimum is then a finite optimisation over these numbers, subject to the condition that the remaining points form cycles of length > r with the correct parity. `scripts/05_validate.py` checks that the two methods agree.

## Tests

```bash
python -m unittest discover -s tests -v
```

## License

MIT, see `LICENSE`.

## Citation

See `CITATION.cff`.
