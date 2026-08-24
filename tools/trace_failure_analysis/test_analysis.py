#!/usr/bin/env python3
"""Focused tests for diagnostic attribution and shared interval projection."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.trace_failure_analysis.diagnostics import (
    DiagnosticError,
    attribute_report,
    load_reports,
    summarize_reports,
)
from tools.trace_failure_analysis.intervals import interval_metrics, project_alignment


class DiagnosticAttributionTests(unittest.TestCase):
    def test_fixture_assigns_every_category_once(self) -> None:
        fixture = ROOT / "tests/data/archive_seeds/diagnostic_categories.json"
        report = load_reports(fixture)[0]
        result = attribute_report(report)
        self.assertEqual(result["truth_bases"], 130)
        self.assertEqual(result["supported_bases"], 0)
        self.assertEqual(result["missing_truth_bases"], 130)
        self.assertIsNone(result["metrics"]["base_precision"])
        self.assertIsNone(result["metrics"]["interval_precision"])
        summary = summarize_reports([report])
        self.assertIsNone(summary["base_precision"])
        self.assertIsNone(summary["interval_precision"])
        self.assertEqual(set(result["failure_counts"]), {
            "candidate_miss",
            "no_retained_seed",
            "no_matching_seed",
            "repetitive_seed_suppression",
            "anchor_cap",
            "no_valid_chain",
            "chain_limit",
            "sequence_budget",
            "alignment_band_failure",
            "alignment_rejection",
            "alignment_cap",
            "mosaic_selection",
            "other",
        })
        self.assertTrue(all(value == 10 for value in result["failure_counts"].values()))
        self.assertEqual(result["stage_totals"]["bytes_read"], 128)

    def test_coordinate_mode_unions_overlapping_truth(self) -> None:
        report = {
            "metagenome_id": "overlap",
            "truth_intervals": [
                {"truth_interval": {"start": 0, "end": 10}, "candidate_selected": True, "mosaic_supported_bases": 0, "failure_category": "other"},
                {"truth_interval": {"start": 5, "end": 15}, "candidate_selected": True, "mosaic_supported_bases": 0, "failure_category": "other"},
            ],
            "supported_intervals": [{"start": 0, "end": 8}],
            "stages": [],
        }
        result = attribute_report(report)
        self.assertEqual(result["accounting_mode"], "coordinate")
        self.assertEqual(result["truth_bases"], 15)
        self.assertEqual(result["supported_bases"], 8)
        self.assertEqual(result["missing_truth_bases"], 7)

        false_positive = dict(report)
        false_positive["supported_intervals"] = [{"start": 0, "end": 20}]
        summary = summarize_reports([false_positive])
        self.assertEqual(summary["observed_bases"], 20)
        self.assertEqual(summary["intersection_bases"], 15)
        self.assertEqual(summary["base_precision"], 0.75)

    def test_coordinate_attribution_splits_boundaries_and_conflicts(self) -> None:
        report = {
            "metagenome_id": "boundary-conflict",
            "truth_intervals": [
                {
                    "truth_interval": {"start": 0, "end": 10},
                    "candidate_selected": False,
                    "mosaic_supported_bases": 0,
                    "failure_category": "candidate_miss",
                },
                {
                    "truth_interval": {"start": 5, "end": 15},
                    "candidate_selected": True,
                    "mosaic_supported_bases": 0,
                    "failure_category": "no_matching_seed",
                },
            ],
            "supported_intervals": [{"start": 0, "end": 2}],
            "stages": [],
        }
        result = attribute_report(report)
        self.assertEqual(
            [(item["interval"], item["category"]) for item in result["attributions"]],
            [
                ({"start": 2, "end": 5}, "candidate_miss"),
                ({"start": 5, "end": 10}, "other"),
                ({"start": 10, "end": 15}, "no_matching_seed"),
            ],
        )

    def test_malformed_boolean_is_rejected(self) -> None:
        for malformed in ("unknown", "2", 2, [], {}):
            with self.subTest(value=malformed), self.assertRaises(DiagnosticError):
                attribute_report(
                    {
                        "metagenome_id": "bad-bool",
                        "truth_intervals": [
                            {
                                "truth_interval": {"start": 0, "end": 1},
                                "candidate_selected": malformed,
                            }
                        ],
                        "stages": [],
                    }
                )

    def test_gap_error_is_symmetric_difference(self) -> None:
        self.assertEqual(
            interval_metrics([(0, 5)], [(0, 6)])["gap_error_bases"],
            1,
        )
        self.assertEqual(
            interval_metrics([(0, 5), (8, 10)], [(0, 6)])["gap_error_bases"],
            3,
        )

    def test_count_mode_rejects_overlapping_truth(self) -> None:
        report = {
            "metagenome_id": "bad-overlap",
            "truth_intervals": [
                {"truth_interval": {"start": 0, "end": 10}, "mosaic_supported_bases": 0},
                {"truth_interval": {"start": 5, "end": 15}, "mosaic_supported_bases": 0},
            ],
            "stages": [],
        }
        with self.assertRaises(DiagnosticError):
            attribute_report(report)

    def test_projection_reuses_query_coordinates(self) -> None:
        result = project_alignment(
            [{"start": 0, "end": 10}],
            10,
            cigar="4=2I4=",
            query_deletion_operation="I",
        )
        self.assertEqual(result["supported_intervals"], [(0, 4), (6, 10)])
        self.assertEqual(result["alignment_deletions"], [(4, 6)])


if __name__ == "__main__":
    unittest.main()
