"""Trace diagnostic and benchmark analysis helpers.

The package is intentionally dependency free.  It reads the JSON form of the
Rust trace diagnostic models and uses the same zero-based, half-open interval
convention as the production normalization scripts.
"""

from .diagnostics import (
    FAILURE_CATEGORIES,
    attribute_report,
    load_reports,
    summarize_reports,
)
from .intervals import (
    IntervalError,
    interval_length,
    interval_metrics,
    subtract_intervals,
    union_intervals,
)

__all__ = [
    "FAILURE_CATEGORIES",
    "IntervalError",
    "attribute_report",
    "interval_length",
    "interval_metrics",
    "load_reports",
    "subtract_intervals",
    "summarize_reports",
    "union_intervals",
]
