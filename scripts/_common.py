"""Shared helpers for the reproduction scripts."""

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESULTS = os.path.join(ROOT, "results")
sys.path.insert(0, ROOT)


def result_path(name: str) -> str:
    os.makedirs(RESULTS, exist_ok=True)
    return os.path.join(RESULTS, name)


def label(lam) -> str:
    return "(" + ",".join(str(p) for p in lam) + ")"
