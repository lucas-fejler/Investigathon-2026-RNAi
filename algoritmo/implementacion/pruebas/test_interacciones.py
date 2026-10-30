"""Coordinate, completeness, cache and fallback controls."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from evaluar_interacciones import pair_key, position_state, priority_tasks, run
from archivos import digest, write_json, write_rows
import rnahybrid as native


class NativePolicyTests(unittest.TestCase):
    """Synthetic controls; no RNAhybrid invocation or candidate measurement."""
    def test_coordinate_and_mask_validation(self):
            """Retained mapping handles antiparallel query and masked dangles.

            Args:
                None.

            Returns:
                None.
            """
            source = {"target_sequence": "ACCCCT", "query_sequence": "GGGG",
                      "query_id": "q", "query_sha256": "qhash",
                      "reference_id": "r", "reference_sha256": "rhash"}
            block = {"native_position": 1, "alignment_rows": [
                "A    U", " CCCC ", " GGGG ", "      "]}
            mapped = native.map_block(block, source, set())
            self.assertEqual(mapped["core_target_interval"], [1, 5])
            self.assertEqual(mapped["paired_edges"][0], [3, 1, "WC"])
            block["alignment_rows"][0] = "N    U"
            with self.assertRaises(ValueError):
                native.map_block(block, source, set())
            self.assertEqual(native.map_block(block, source, {0})[
                "masked_dangle_characters"], 1)
            block["native_position"] = 0
            with self.assertRaises(ValueError):
                native.map_block(block, source, set())

    def test_incomplete_fallback_and_failures_stay_unavailable(self):
            """An empty retrieval, prefix cap or failed dependency cannot give P=0.

            Args:
                None.

            Returns:
                None.
            """
            negative = {"state": "negative_reported_pair", "flag": False}
            positive = {"state": "positive_witness", "flag": True}
            self.assertEqual(position_state([], 4),
                             ("untested_or_incomplete", None))
            self.assertEqual(position_state([negative], 4)[1], None)
            self.assertEqual(position_state([negative] * 4, 4),
                             ("negative_complete", False))
            self.assertEqual(position_state([positive], 4),
                             ("positive_complete", True))
            for state in ("failed", "unavailable_cap", "unavailable_interrupted",
                          "unavailable_ambiguity", "omitted"):
                self.assertEqual(position_state([positive, {"state": state}], 4),
                                 ("unavailable_required_dependency", None))

    def test_evidence_reuse_requires_all_identities(self):
            """Query, reference, evaluator and classifier identities bind reuse.

            Args:
                None.

            Returns:
                None.
            """
            query = {"query_id": "q", "sequence_sha256": "qsha"}
            reference = {"reference_id": "r", "sequence_sha256": "rsha"}
            binding = {"evaluator_id": "e", "classification_id": "c"}
            key = pair_key(query, reference, binding)
            self.assertEqual(key, pair_key(dict(query), dict(reference),
                                           dict(binding)))
            for changed in ({**binding, "evaluator_id": "new"},
                            {**binding, "classification_id": "new"}):
                self.assertNotEqual(key, pair_key(query, reference, changed))
            self.assertNotEqual(key, pair_key({**query, "sequence_sha256": "new"},
                                             reference, binding))
            self.assertNotEqual(key, pair_key(
                query, {**reference, "sequence_sha256": "new"}, binding))

    def test_strict_native_truncation_and_missing_audit_end(self):
            """Malformed output and incomplete native work cannot certify no hit.

            Args:
                None.

            Returns:
                None.
            """
            with self.assertRaises(ValueError):
                native.parse_compact(b"truncated", {}, {})
            source = {"query_length": 1, "target_length": 1,
                      "query_sequence": "A", "target_sequence": "T"}
            result = native.parse_audit(b"BEGIN\tQ\tT\t1\t1\tA\tT\n",
                                        {("Q", "T"): source})
            self.assertTrue(result["errors"])

    def test_empty_priority_actually_executes_fallback(self):
            """The runner invokes remaining scope when retrieval offers no pair.

            Args:
                None.

            Returns:
                None.
            """
            manifest = {"panel_id": "synthetic_execution_control"}
            binding = {"evaluator_id": "synthetic"}
            query = {"query_id": "q", "sequence": "A" * 21}
            window = {"window_id": "w", "parent_id": "fixture@hash",
                      "occurrences": [{"query_id": "q", "orientation": "+"}]}
            reference = {"reference_id": "r", "length": 21}
            result = {"query_id": "q", "reference_id": "r",
                      "evidence_key": "synthetic", "state": "positive_witness",
                      "flag": True}
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                write_json(output / "execution_binding.json", binding)
                write_rows(output / "selected_queries.jsonl", [query])
                write_rows(output / "selected_windows.jsonl", [window])
                (output / "priority.tsv").write_text("")
                write_json(output / "priority_identity.json", {
                    "priority_sha256": digest(output / "priority.tsv")})
                with (
                    patch("evaluar_interacciones.validate_manifest",
                          return_value=manifest),
                    patch("evaluar_interacciones.catalog",
                          return_value={"r": reference}),
                    patch("evaluar_interacciones.build_identity", return_value=binding),
                    patch("evaluar_interacciones.references",
                          return_value=[("r", "T" * 21)]),
                    patch("evaluar_interacciones.assess_pair",
                          return_value=result) as native_call,
                    patch("evaluar_interacciones.priority_tasks", return_value=[]),
                ):
                    summary = run(output, 10, 2)
                self.assertEqual(native_call.call_count, 1)
                self.assertEqual(summary["positions"], {"positive_complete": 1})

    def test_priority_retains_intended_records_and_both_orientations(self):
            """Nonmember priority changes ordering, never the required task set.

            Args:
                None.

            Returns:
                None.
            """
            window = {"parent_id": "fixture@hash", "occurrences": [
                {"query_id": "plus", "orientation": "+"},
                {"query_id": "minus", "orientation": "-"}]}
            choices = {"plus": [("own", 21), ("other", 18)],
                       "minus": [("own", 20), ("other", 19)]}
            records = {"own": {"length": 100}, "other": {"length": 200}}
            config = {"candidates": {"fixture": {"curation": {
                "intended_targets": [{"reference_id": "own"}]}}}}
            tasks = priority_tasks(window, choices, records, config)
            self.assertEqual(tasks, [("minus", "other", 19),
                                     ("plus", "other", 18),
                                     ("plus", "own", 21),
                                     ("minus", "own", 20)])
