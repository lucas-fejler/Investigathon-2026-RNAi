"""Exact ranking examples and multi-site aggregation invariants."""
import os
import sys
import unittest
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ranking
import puntuacion as scoring
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURES = os.path.join(os.path.dirname(__file__), "datos")
LAMBDAS = (Fraction(1, 2), Fraction(1), Fraction(2))

def synthetic_region(candidate_id, start, favorable, flagged,
                     status="scored"):
    """Build a labelled synthetic scored-region record for ranking tests.

    These exercise the ranking functions with component counts that a
    window-level fixture cannot reach, because consecutive regions share 179 of
    180 positions and so differ by at most one.

    Args:
        candidate_id: Synthetic candidate identifier.
        start: One-based region start coordinate.
        favorable: Favorable position count.
        flagged: Flagged position count.
        status: Region status, "scored" or "withheld".

    Returns:
        Region record shaped like `puntuacion.score_region` output.
    """
    parent = f"{candidate_id}@{'0' * 64}"
    return {
        "region_id": f"region:{parent}:{start - 1}:{start - 1 + 200}",
        "candidate_id": candidate_id,
        "parent_id": parent,
        "interval_1based_inclusive": [start, start + 199],
        "position_count": 180,
        "denominator": 180,
        "favorable": favorable,
        "flagged": flagged,
        "a_unavailable_positions": 0,
        "p_unavailable_positions": 0,
        "a_exact": scoring.exact_string(scoring.activity(favorable)),
        "p_exact": scoring.exact_string(scoring.penalty(flagged)),
        "scores": scoring.score_block(favorable, flagged),
        "status": status,
        "withheld_causes": None if status == "scored" else ["a_unavailable"],
        "fixture": "synthetic",
    }


class ExactScoringTest(unittest.TestCase):
    """Worked examples and exact-arithmetic properties (W01-W06)."""

    def test_w01_worked_examples(self):
        """W01: design examples E1-E5 at all three lambda values."""
        cases = [
            # favorable, flagged, S_0.5, S_1, S_2
            (90, 30, Fraction(75, 180), Fraction(60, 180), Fraction(30, 180)),
            (100, 0, Fraction(100, 180), Fraction(100, 180),
             Fraction(100, 180)),
            (60, 40, Fraction(40, 180), Fraction(20, 180), Fraction(-20, 180)),
            (0, 180, Fraction(-90, 180), Fraction(-180, 180),
             Fraction(-360, 180)),
            (45, 15, Fraction(75, 360), Fraction(30, 180), Fraction(15, 180)),
        ]
        for favorable, flagged, half, one, two in cases:
            self.assertEqual(scoring.score(favorable, flagged, LAMBDAS[0]),
                             half, (favorable, flagged, "0.5"))
            self.assertEqual(scoring.score(favorable, flagged, LAMBDAS[1]),
                             one, (favorable, flagged, "1"))
            self.assertEqual(scoring.score(favorable, flagged, LAMBDAS[2]),
                             two, (favorable, flagged, "2"))

    def test_w01_components(self):
        """W01: A and P use the fixed denominator and reduce exactly."""
        self.assertEqual(scoring.activity(90), Fraction(1, 2))
        self.assertEqual(scoring.penalty(180), Fraction(1))
        self.assertEqual(scoring.penalty(0), Fraction(0))
        self.assertEqual(
            scoring.exact_string(scoring.activity(79)), "79/180")

    def test_w02_negative_scores_order_below_zero(self):
        """W02: negative scores sort correctly against each other and zero."""
        records = [synthetic_region("SYN", 1, 60, 40),
                   synthetic_region("SYN", 301, 0, 180),
                   synthetic_region("SYN", 601, 90, 0)]
        rows = ranking.rank_gene(records, "2")
        self.assertEqual([r["score_exact"] for r in rows],
                         ["1/2", "-1/9", "-2/1"])
        self.assertEqual([r["rank"] for r in rows], [1, 2, 3])
        self.assertLess(Fraction(-1, 9), Fraction(0))
        self.assertLess(Fraction(-2), Fraction(-1, 9))

    def test_w03_exact_tie_from_different_components(self):
        """W03: (90,30) and (80,20) tie at lambda=1 and stay tied."""
        first = scoring.score(90, 30, Fraction(1))
        second = scoring.score(80, 20, Fraction(1))
        self.assertEqual(first, second)
        rows = ranking.rank_gene(
            [synthetic_region("SYN", 1, 90, 30),
             synthetic_region("SYN", 301, 80, 20)], "1")
        self.assertEqual([r["rank"] for r in rows], [1, 1])
        self.assertEqual([r["tie_group_size"] for r in rows], [2, 2])
        self.assertTrue(all(r["is_co_best"] for r in rows))

    def test_w04_tie_separates_and_reverses_across_lambda(self):
        """W04: the W03 pair separates at 0.5 and reverses at 2."""
        records = [synthetic_region("SYN", 1, 90, 30),
                   synthetic_region("SYN", 301, 80, 20)]
        half = ranking.rank_gene(records, "0.5")
        self.assertEqual(half[0]["favorable"], 90)
        self.assertEqual([r["score_exact"] for r in half],
                         ["5/12", "7/18"])
        self.assertEqual([r["rank"] for r in half], [1, 2])

        two = ranking.rank_gene(records, "2")
        self.assertEqual(two[0]["favorable"], 80)
        self.assertEqual([r["score_exact"] for r in two], ["2/9", "1/6"])
        self.assertEqual([r["rank"] for r in two], [1, 2])

    def test_w05_float_would_miss_this_tie(self):
        """W05: binary floats do not see the W03 tie; Fractions do."""
        float_first = 90 / 180 - 1 * (30 / 180)
        float_second = 80 / 180 - 1 * (20 / 180)
        self.assertNotEqual(float_first, float_second)
        self.assertEqual(scoring.score(90, 30, Fraction(1)),
                         scoring.score(80, 20, Fraction(1)))

    def test_w06_half_lambda_yields_denominator_360(self):
        """W06: lambda=1/2 legitimately produces denominator 360."""
        value = scoring.score(45, 15, Fraction(1, 2))
        self.assertEqual(value, Fraction(75, 360))
        self.assertEqual(value.denominator, 24)
        block = scoring.score_block(45, 15)
        self.assertEqual(block["0.5"]["exact"], "5/24")
        self.assertEqual(Fraction(block["0.5"]["numerator"],
                                  block["0.5"]["denominator"]),
                         Fraction(75, 360))

    def test_display_is_never_used_for_comparison(self):
        """Display decimals round without floats and stay display-only."""
        self.assertEqual(scoring.display_string(Fraction(1, 3)), "0.333333")
        self.assertEqual(scoring.display_string(Fraction(-1, 3)), "-0.333333")
        self.assertEqual(scoring.display_string(Fraction(0)), "0.000000")
        # Two distinct scores can share a display string; the exact strings
        # must still differ, which is why ties are decided on exact values.
        first = Fraction(1, 3)
        second = Fraction(1, 3) + Fraction(1, 10 ** 9)
        self.assertEqual(scoring.display_string(first),
                         scoring.display_string(second))
        self.assertNotEqual(scoring.exact_string(first),
                            scoring.exact_string(second))


class AggregationTest(unittest.TestCase):
    """Per-position aggregation and denominator invariants (W07, W08)."""

    @classmethod
    def setUpClass(cls):
        """Load the synthetic fixture inputs once.

        Args:
            None.

        Returns:
            None.
        """
        cls.findings = scoring.Findings()
        cls.a_inputs = scoring.load_activity_inputs(
            FIXTURES, cls.findings, "synthetic_a")
        cls.p_inputs = scoring.load_penalty_inputs(
            FIXTURES, cls.findings, "synthetic_p")

    def test_w07_multi_site_position_counts_once(self):
        """W07: several sites across organisms contribute one unit to P."""
        positions = self.p_inputs["positions_by_window"].values()
        flagged_rows = [r for r in positions if r.get("flagged")]
        self.assertTrue(flagged_rows)
        rich = [r for r in flagged_rows
                if r["qualifying_site_count"] > 1 and len(r["organisms"]) > 1]
        self.assertTrue(rich, "fixture should contain a multi-site position")

        region = self.a_inputs["regions_by_id"][
            "region:" + rich[0]["parent_id"] + ":0:200"]
        findings = scoring.Findings()
        record = scoring.score_region(region, self.a_inputs, self.p_inputs,
                                      findings)
        # SYN_GENE_A region 0 has 30 flagged windows, each carrying 7 sites
        # across two organisms: P must be 30/180, never 210/180.
        self.assertEqual(record["flagged"], 30)
        self.assertEqual(record["p_exact"], "1/6")
        self.assertEqual(scoring.penalty(record["flagged"]),
                         Fraction(30, 180))
        self.assertEqual(record["penalty_evidence"]["qualifying_site_total"],
                         210)
        self.assertEqual(record["penalty_evidence"]["organisms"],
                         ["Apis mellifera", "Homo sapiens"])

    def test_w07_regional_counts_match_closed_form(self):
        """W07: fixture regional counts equal the hand-derived formula."""
        findings = scoring.Findings()
        parent = None
        for record in self.a_inputs["regions_by_id"].values():
            if record["candidate_id"] == "SYN_GENE_A":
                parent = record["parent_id"]
                break
        for start in range(4):
            region = self.a_inputs["regions_by_id"][
                f"region:{parent}:{start}:{start + 200}"]
            record = scoring.score_region(region, self.a_inputs,
                                          self.p_inputs, findings)
            # favorable(s) = 90 - s and flagged(s) = 30 - s by construction.
            self.assertEqual(record["favorable"], 90 - start)
            self.assertEqual(record["flagged"], 30 - start)
            self.assertEqual(record["scores"]["1"]["exact"], "1/3")
        self.assertEqual(findings.items, [])

    def test_w08_denominator_fixed_under_failure(self):
        """W08: the denominator stays 180 for scored and withheld regions."""
        findings = scoring.Findings()
        for region in self.a_inputs["regions_by_id"].values():
            record = scoring.score_region(region, self.a_inputs,
                                          self.p_inputs, findings)
            self.assertEqual(record["denominator"], 180)
            self.assertEqual(record["position_count"], 180)
