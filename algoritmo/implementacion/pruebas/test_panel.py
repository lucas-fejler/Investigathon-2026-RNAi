"""Release semantics with explicitly synthetic, in-memory test measurements."""

from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import ranking
import puntuacion
from cobertura import decision, position_coverage, search_summary
from consolidar_resultados import independently_score


def synthetic_fixture(flag=True):
    """Build three tied synthetic regions without writing candidate outputs.

    Args:
        flag: Synthetic P position flag.

    Returns:
        Frozen-shaped manifest, A/P measurements and producer M8 outputs.
    """
    pid = "fixture@" + "0" * 64
    regions, activity, positions, scored = [], {}, {}, []
    for start in range(3):
        rid = f"region:{pid}:{start}:{start + 200}"
        wids = [f"window:{pid}:{s}:{s + 21}"
                for s in range(start, start + 180)]
        row = {"region_id": rid, "parent_id": pid, "candidate_id": "fixture",
               "interval_1based_inclusive": [start + 1, start + 200],
               "window_ids": wids, "a_numerator": 90,
               "a_denominator": 180, "a_exact": "90/180"}
        regions.append(row)
        for offset, wid in zip(range(start, start + 180), wids):
            activity[wid] = {"status": "measured", "classification":
                             "favorable" if offset % 2 == 0 else "unfavorable"}
            positions[wid] = {"status": "positive_complete" if flag else
                              "negative_complete", "flagged": flag}
        scored.append({
            **row, "status": "scored", "favorable": 90,
            "flagged": 180 if flag else 0, "a_exact": "1/2",
            "p_exact": "1" if flag else "0",
            "scores": puntuacion.score_block(90, 180 if flag else 0)})
    manifest = {"regions": regions, "genes": [{
        "candidate_id": "fixture", "admitted_count": 3,
        "outside_panel_count": 0,
        "metadata": {"curation": {"organism": "synthetic fixture"}}}]}
    rankings = [r for key in ("0.5", "1", "2")
                for r in ranking.rank_gene(scored, key)]
    return (manifest, activity, positions,
            {r["region_id"]: r for r in scored}, rankings)


class PanelReleaseTests(unittest.TestCase):
    """Verify independent arithmetic, withholding and species scope."""

    def test_saturated_P_and_exact_ties_are_verified(self):
        """All admitted P=1 is explicit and preserves exact M8 tie groups.

        Args:
            None.

        Returns:
            None.
        """
        rows, genes, ranks, ties = independently_score(*synthetic_fixture())
        self.assertTrue(genes[0]["ranking_eligible"])
        self.assertEqual(genes[0]["constant_p_exact"], "1")
        self.assertEqual(genes[0]["completed_p_numerator_distribution"],
                         {180: 3})
        self.assertEqual(len(ranks), 9)
        self.assertTrue(all(len(t["region_ids"]) == 3 for t in ties))
        for row in rows:
            self.assertEqual(row["p_numerator"], 180)
            self.assertEqual(Fraction(row["scores"]["0.5"]), 0)
            self.assertEqual(Fraction(row["scores"]["1"]), Fraction(-1, 2))

    def test_constant_zero_P_is_distinct_from_missing_P(self):
        """A synthetic complete negative exports zero rather than null.

        Args:
            None.

        Returns:
            None.
        """
        rows, genes, _, _ = independently_score(*synthetic_fixture(False))
        self.assertEqual(genes[0]["constant_p_exact"], "0")
        self.assertTrue(all(r["p_complete"] and r["p_numerator"] == 0
                            for r in rows))

    def test_missing_P_withholds_whole_gene_and_constant_claim(self):
        """A and P remain separate and no partial panel releases a ranking.

        Args:
            None.

        Returns:
            None.
        """
        manifest, a, p, producer, ranks = synthetic_fixture()
        del p[manifest["regions"][0]["window_ids"][0]]
        rows, genes, released, ties = independently_score(
            manifest, a, p, producer, ranks)
        self.assertEqual(genes[0]["a_complete_regions"], 3)
        self.assertEqual(genes[0]["p_complete_regions"], 2)
        self.assertIsNone(genes[0]["p_constant_across_complete_panel"])
        self.assertFalse(genes[0]["ranking_eligible"])
        self.assertEqual((released, ties), ([], []))
        self.assertIsNone(rows[0]["p_numerator"])
        self.assertIsNone(rows[0]["scores"])
        self.assertTrue(rows[0]["a_complete"])

    def test_changed_flags_A_scores_or_ties_prevent_release(self):
        """Independent checks catch mismatched flags, M6 counts and M8 ties.

        Args:
            None.

        Returns:
            None.
        """
        for field in ("flag", "A", "score", "tie"):
            manifest, a, p, producer, ranks = synthetic_fixture()
            wid = manifest["regions"][0]["window_ids"][0]
            rid = manifest["regions"][0]["region_id"]
            if field == "flag":
                p[wid]["flagged"] = False
            elif field == "A":
                a[wid]["classification"] = "unfavorable"
            elif field == "score":
                producer[rid]["scores"]["1"]["exact"] = "999/180"
            else:
                ranks[0]["tie_group_size"] = 1
            _, genes, released, _ = independently_score(
                manifest, a, p, producer, ranks)
            self.assertFalse(genes[0]["ranking_eligible"], field)
            self.assertEqual(released, [])

    def test_complete_negative_requires_entire_greedy_policy_scope(self):
        """Partial or capped searches and repeated observations are not zero.

        Args:
            None.

        Returns:
            None.
        """
        pair = {"query_id": "q1", "reference_id": "s|r",
                "state": "negative_reported_pair", "flag": False,
                "greedy_complete": True}
        other = {**pair, "query_id": "q2"}
        self.assertEqual(decision([pair], 2), ("unresolved", None))
        self.assertEqual(decision([pair, pair], 2), ("unresolved", None))
        self.assertEqual(decision([pair, other], 2),
                         ("negative_complete", False))
        other["greedy_complete"] = False
        self.assertEqual(decision([pair, other], 2), ("unresolved", None))
        other.update(state="unavailable_timeout", flag=None)
        self.assertEqual(decision([pair, other], 2), ("failed", None))

    def test_union_positive_leaves_other_species_untested(self):
        """Positive stopping never certifies a negative in unsearched species.

        Args:
            None.

        Returns:
            None.
        """
        windows = [{"window_id": "w", "parent_id": "fixture@hash",
                    "interval_1based_inclusive": [1, 21],
                    "occurrences": [{"query_id": "q1"}, {"query_id": "q2"}]}]
        datasets = [{"id": s, "organism": s, "statistics": {"records": 1}}
                    for s in ("one", "two")]
        pair = {"evidence_key": "synthetic", "query_id": "q1",
                "reference_id": "one|r", "dataset_id": "one",
                "state": "positive_witness", "flag": True,
                "greedy_complete": False}
        pos, species, membership = position_coverage(
            windows, [pair], datasets, {"w": {"status": "measured",
                                             "flagged": True}})
        self.assertEqual(pos[0]["status"], "positive_complete")
        self.assertEqual(species[1]["status"], "untested")
        index = {(r["window_id"], r["dataset_id"]): r for r in species}
        summary = search_summary({"w"}, {"w": pos[0]}, index, membership,
                                  [pair], datasets, "panel", "synthetic")
        self.assertEqual(summary[1]["evaluated_pairs"], 0)
        self.assertEqual(summary[1]["not_required_for_union_P_pairs"], 2)
        self.assertFalse(summary[1]["all_required_pairs_evaluated"])
        self.assertEqual(summary[1]["negative_complete_positions"], 0)
        negative = {"evidence_key": "synthetic-negative", "query_id": "q2",
                    "reference_id": "two|r", "dataset_id": "two",
                    "state": "negative_reported_pair", "flag": False,
                    "greedy_complete": True}
        pos, species, _ = position_coverage(
            windows, [pair, negative], datasets,
            {"w": {"status": "measured", "flagged": True}})
        self.assertEqual(pos[0]["status"], "positive_complete")
        self.assertEqual(species[1]["status"], "partial_unresolved")
        self.assertEqual(species[1]["evaluated_pairs"], 1)
        self.assertEqual(species[1]["unattempted_pairs"], 1)


if __name__ == "__main__":
    unittest.main()
