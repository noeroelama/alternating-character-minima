# alternating-character-minima

Code and data accompanying the manuscript

> A. Nurwahid, *Minimum values of irreducible characters of alternating groups* (manuscript).

For an irreducible character χ of the alternating group A<sub>n</sub>, let m(χ) be the smallest real value of χ. The scripts in this repository reproduce every table and every computer verification in the manuscript.

## Requirements

* Python 3.9 or newer
* No third-party packages (standard library only)

## Reproducing the results

```bash
python run_all.py
```

This runs the unit tests and then the scripts below in order. The total running time is about 6 minutes on a laptop, and script 01 needs about 1.5 GB of memory. Output files are written to `results/`; the copies committed in this repository were produced with the same command.

| Script | Manuscript | Output |
|---|---|---|
| `scripts/01_threshold.py` | Theorem A(d), Remark 3.4, Table S1 | `results/threshold.txt` |
| `scripts/02_exact_minima.py` | Theorem C (exact minima, checked for n ≤ 160, resp. n ≤ 120) | `results/exact_minima.txt` |
| `scripts/03_table_c_mu.py` | Section 7 (counts for \|μ\| ≤ 8) and Table S3 (constants c<sub>μ</sub>, \|μ\| ≤ 6) | `results/table_c_mu.csv`, `results/table_c_mu_summary.txt` |
| `scripts/04_asymptotics.py` | Table S4 (convergence in Theorem B) | `results/table_ratio.csv` |
| `scripts/05_validate.py` | Supplementary Section S5 (cross-check of the two algorithms) | `results/validation.txt` |
| `scripts/06_open_problems.py` | Section 7 (cases with c<sub>μ</sub> = 0) | `results/open_problems.csv` |
| `scripts/07_self_conjugate.py` | Lemma 3.2, Proposition 3.3, Table S2 (split characters) | `results/self_conjugate.txt` |
| `scripts/08_czero.py` | Proposition 1.2(ii), Remark 7.1, Theorem D (case c<sub>μ</sub> = 0) | `results/czero.txt` |

All scripts except 03 accept `--help`; the ranges of n can be changed from the command line.

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
