"""Minimum values of irreducible characters of the alternating groups."""

from .alternating import (MinimumRecord, all_minima, irreducible_labels,
                          minimum_record, number_of_classes,
                          number_of_irreducibles, sum_of_squared_degrees)
from .charpoly import character_polynomial, evaluate, minimum_fixed_level
from .hfunction import c_mu, h_coefficients, h_value
from .mn import degree, sn_character
from .partitions import conjugate, partitions

__all__ = [
    "MinimumRecord", "all_minima", "irreducible_labels", "minimum_record",
    "number_of_classes", "number_of_irreducibles", "sum_of_squared_degrees",
    "character_polynomial", "evaluate", "minimum_fixed_level",
    "c_mu", "h_coefficients", "h_value", "degree", "sn_character",
    "conjugate", "partitions",
]

__version__ = "0.2.0"
