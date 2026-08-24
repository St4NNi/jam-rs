#!/usr/bin/env python3
"""Focused self-tests for the archive-seed experiment protocol."""

from __future__ import annotations

import unittest

from common import ExperimentError, ROOT, _measurement_envelope, load_manifest, run_manifest, workspace_path


class ExperimentProtocolTests(unittest.TestCase):
    def test_dry_run_is_bounded_and_has_measurement_fields(self) -> None:
        manifest = {
            "schema_version": "archive-seeds-experiment-manifest-v1",
            "name": "unit-plan",
            "kind": "seed_screen",
            "cases": [
                {
                    "name": "one",
                    "command": ["jam", "trace"],
                    "metadata": {"output_paths": ["tools/out/unit-plan.json"]},
                }
            ],
        }
        result = run_manifest(manifest, execute=False)
        self.assertEqual(result["mode"], "plan")
        measurement = result["cases"][0]["measurement"]
        self.assertEqual(measurement["timing"]["wall_micros"], 0)
        self.assertIn("metrics", measurement)
        self.assertIn("stages", measurement)
        self.assertIn("requests", measurement)
        self.assertEqual(result["cases"][0]["stage_metrics"], [])
        self.assertEqual(result["cases"][0]["rss_max_kib"], 0)
        self.assertEqual(result["tmp_policy"], "workspace-only")

    def test_workspace_rejects_external_paths(self) -> None:
        with self.assertRaises(ExperimentError):
            workspace_path("/tmp/archive-seeds-out.json", field="output")

    def test_stage_and_request_counters_are_normalized(self) -> None:
        envelope = _measurement_envelope(
            {
                "stages": [{"stage": 2, "name": "rescue", "bytes_read": 13}],
                "requests": [{"offset": 4, "length": 8, "bytes_returned": 8}],
                "resource_metrics": {"decoded_bytes": 21},
            },
            wall_micros=17,
            cpu_micros=11,
            rss_max_kib=99,
        )
        self.assertEqual(envelope["stages"][0]["bytes_read"], 13)
        self.assertEqual(envelope["requests"][0]["bytes_returned"], 8)
        self.assertEqual(envelope["metrics"]["decoded_bytes"], 21)
        self.assertEqual(envelope["rss_max_kib"], 99)

    def test_checked_in_manifest_has_expected_kind(self) -> None:
        manifest = load_manifest(
            ROOT / "evaluation/archive-seeds/manifests/seed-screen.example.json",
            "seed_screen",
        )
        self.assertGreaterEqual(len(manifest["cases"]), 1)


if __name__ == "__main__":
    unittest.main()
