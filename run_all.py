"""Run the test suite and every reproduction script in order.

Usage:  python run_all.py
"""

import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    "01_threshold.py",
    "02_exact_minima.py",
    "03_table_c_mu.py",
    "04_asymptotics.py",
    "05_validate.py",
    "06_open_problems.py",
]


def run(cmd) -> None:
    start = time.time()
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)
    print(f"  ({time.time() - start:.1f} s)\n", flush=True)


def main() -> None:
    run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"])
    for script in SCRIPTS:
        run([sys.executable, os.path.join("scripts", script)])


if __name__ == "__main__":
    main()
