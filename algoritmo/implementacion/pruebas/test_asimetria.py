"""Energy parameters, author oracle and aggregation invariants."""
import hashlib
import json
import os
import random
import sqlite3
import sys
import tempfile
import unittest
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import asimetria as m6
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENERGY_EXPORT = os.path.join(REPO, "parametros_energia", "dsrip_terminal_energies.tsv")
ENERGY_SOURCE = os.path.join(REPO, "parametros_energia", "dsrip_lookup.db")
PUBLISHED_INITIATION = 410
PUBLISHED_AU_END_ON_AU = 22
PUBLISHED_AU_END_ON_CG = 44
PUBLISHED_SYMMETRY = 43
PUBLISHED_STACKS = {"AA": -94, "AC": -225, "AG": -201, "AU": -109,
                    "CA": -207, "CC": -328, "CG": -233, "GA": -242,
                    "GC": -346, "UA": -129}

def _sha256_of(path):
    """Hash a file's bytes.

    Args:
        path: Path to an existing file.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _read_json(path):
    """Read a JSON document.

    Args:
        path: Path to a JSON file.

    Returns:
        Parsed document.
    """
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _collapse(dinucleotide):
    """Collapse a stack to its reverse-complement-equivalent representative.

    Args:
        dinucleotide: Two-nucleotide RNA string.

    Returns:
        Lexicographically smaller of the stack and its reverse complement.
    """
    return min(dinucleotide, m6.reverse_complement_rna(dinucleotide))


def zuber_energy(key):
    """Reconstruct a 4-mer duplex energy from published parameters.

    This is the independent derivation used by V01: it never consults the
    vendored table.

    Args:
        key: RNA 4-mer.

    Returns:
        Signed integer hundredths of a kcal/mol.
    """
    total = PUBLISHED_INITIATION
    for index in range(3):
        total += PUBLISHED_STACKS[_collapse(key[index:index + 2])]
    for terminal, penultimate in ((key[0], key[1]), (key[3], key[2])):
        if terminal in "AU":
            total += (PUBLISHED_AU_END_ON_AU if penultimate in "AU"
                      else PUBLISHED_AU_END_ON_CG)
    if key == m6.reverse_complement_rna(key):
        total += PUBLISHED_SYMMETRY
    return total


def dsrip_reference_asym(target_sequence, table):
    """Reimplement dsRIP's own duplex construction as an independent oracle.

    Mirrors `generate_siRNA`, `kmer_pairs_np` and `thermoAsym` from
    `dsRNA_efficiency.py` (figshare 10.6084/m9.figshare.28156787 v3), including
    the literal "XX" padding of the passenger 3' overhang.

    Args:
        target_sequence: Sense target sequence as RNA.
        table: Mapping from 4-mer key to integer hundredths.

    Returns:
        List of (guide, score_hundredths, passenger_3_5) per window, ordered by
        dsRIP's own window index.
    """
    complement = {"A": "U", "U": "A", "G": "C", "C": "G", "X": "X"}

    def comp(seq):
        """Complement an RNA sequence without reversing it.

        Args:
            seq: RNA sequence, possibly containing dsRIP's "X" padding.

        Returns:
            Base-by-base complement in the same order.
        """
        return "".join(complement[c] for c in seq)

    def revcomp(seq):
        """Reverse-complement an RNA sequence.

        Args:
            seq: RNA sequence, possibly containing dsRIP's "X" padding.

        Returns:
            Reverse complement.
        """
        return "".join(complement[c] for c in reversed(seq))

    whole_antisense = revcomp(target_sequence)
    padded_antisense = "XX" + whole_antisense
    padded_sense = "XX" + comp(whole_antisense)
    width = m6.WINDOW_LENGTH + 2
    results = []
    for index in range(len(target_sequence) - m6.WINDOW_LENGTH + 1):
        antisense = padded_antisense[index:index + width][2:]
        sense = padded_sense[index:index + width][:-2]
        key_a = antisense[:4]
        key_b = (sense[-4:])[::-1]
        results.append((antisense, table[key_a] - table[key_b], sense))
    return results


class EnergyTableTest(unittest.TestCase):
    """Energy-table provenance and structural invariants (V01-V02, V23)."""

    @classmethod
    def setUpClass(cls):
        """Load the vendored table once for the class.

        Args:
            None.

        Returns:
            None.
        """
        cls.table, cls.identity = m6.load_energy_table(
            ENERGY_EXPORT, ENERGY_SOURCE)

    def test_v01_table_reconstructs_from_published_parameters(self):
        """V01: every entry equals the published Zuber 2022 additive model."""
        for key, value in self.table.items():
            self.assertEqual(
                value, zuber_energy(key),
                f"{key} differs from the published-parameter reconstruction")

    def test_v02_structural_invariants(self):
        """V02: completeness, value range and reverse-complement symmetry."""
        self.assertEqual(len(self.table), m6.EXPECTED_TABLE_ENTRIES)
        self.assertEqual(min(self.table.values()), -592)
        self.assertEqual(max(self.table.values()), 200)
        self.assertEqual(sum(1 for v in self.table.values() if v > 0), 32)
        for key, value in self.table.items():
            self.assertEqual(self.table[m6.reverse_complement_rna(key)], value)

    def test_v23_export_matches_original_and_checksums(self):
        """V23: the text export reproduces the author's original table."""
        connection = sqlite3.connect(
            f"file:{ENERGY_SOURCE}?mode=ro", uri=True)
        rows = connection.execute("SELECT key,value FROM lookup").fetchall()
        self.assertEqual(len(rows), m6.EXPECTED_TABLE_ENTRIES)
        for key, value in rows:
            self.assertEqual(self.table[key], round(value * 100))
        self.assertEqual(self.identity["energy_table_source_md5"],
                         "e7e42098f352e74566850036b61f6433")
        with open(ENERGY_EXPORT, "rb") as handle:
            expected = hashlib.sha256(handle.read()).hexdigest()
        self.assertEqual(self.identity["energy_table_export_sha256"], expected)

    def test_table_rejects_incomplete_and_asymmetric_input(self):
        """A truncated or tampered table must fail loudly, never fall back."""
        with self.assertRaises(m6.EnergyTableError):
            m6.parse_energy_text("key\tdelta_g_37_kcal_per_mol\nAAAA\t1.7\n")
        with self.assertRaises(m6.EnergyTableError):
            m6.validate_energy_table({"AAAA": 172})
        tampered = dict(self.table)
        tampered["AAAA"] = tampered["AAAA"] + 1
        with self.assertRaises(m6.EnergyTableError):
            m6.validate_energy_table(tampered)


class CalculationTest(unittest.TestCase):
    """Checked numerical examples and the sign convention (V03-V09)."""

    @classmethod
    def setUpClass(cls):
        """Load the vendored table and its export checksum once.

        Args:
            None.

        Returns:
            None.
        """
        cls.table, cls.identity = m6.load_energy_table(
            ENERGY_EXPORT, ENERGY_SOURCE)
        cls.sha = cls.identity["energy_table_export_sha256"]

    def measure(self, window):
        """Measure a sense window with the loaded table.

        Args:
            window: 21 nt sense window.

        Returns:
            Measurement record.
        """
        return m6.measure_window(window, self.table, self.sha)

    def test_v03_favorable_hand_calculation(self):
        """V03: DM_ETH window 1-21, CAAA minus GGCU equals +4.58."""
        record = self.measure("ATGGCTTGCTACACTCTTTTG")
        self.assertEqual(record["guide_sequence"], "CAAAAGAGTGTAGCAAGCCAT")
        self.assertEqual(record["end_a_key"], "CAAA")
        self.assertEqual(record["end_b_key"], "GGCU")
        self.assertEqual(record["end_a_energy_kcal_mol"], "0.37")
        self.assertEqual(record["end_b_energy_kcal_mol"], "-4.21")
        self.assertEqual(record["asym_score_kcal_mol"], "4.58")
        self.assertEqual(record["classification"], "favorable")

    def test_v04_unfavorable_near_tie_is_exact(self):
        """V04: -2.15 minus -2.09 persists as exactly -0.06."""
        record = self.measure("TGGCTTGCTACACTCTTTTGC")
        self.assertEqual(record["end_a_energy_kcal_mol"], "-2.15")
        self.assertEqual(record["end_b_energy_kcal_mol"], "-2.09")
        self.assertEqual(record["asym_score_kcal_mol"], "-0.06")
        self.assertEqual(record["classification"], "unfavorable")
        # The exact string must not carry binary-float drift.
        self.assertNotIn("0000", record["asym_score_kcal_mol"])
        self.assertEqual(repr(-2.15 - (-2.09)), "-0.06000000000000005")
        # The hazard is a wide band, not a tiny epsilon: 1e-9 leaves the sign
        # intact, while any band of 0.06 or wider would reclassify this case.
        self.assertLess(record["asym_score_hundredths"], 0)
        self.assertGreaterEqual(abs(record["asym_score_hundredths"]), 6)

    def test_v05_neutral_from_real_sequence(self):
        """V05: DM_ETH window 34-54 scores exactly zero."""
        record = self.measure("GTGGTCTGCAGTCTTCTGACC")
        self.assertEqual(record["end_a_key"], "GGUC")
        self.assertEqual(record["end_b_key"], "GGUC")
        self.assertEqual(record["asym_score_kcal_mol"], "0.00")
        self.assertEqual(record["asym_score_hundredths"], 0)
        self.assertEqual(record["classification"], "neutral")

    def test_v06_synthetic_extremes(self):
        """V06: labelled synthetic fixtures at the attainable extremes."""
        self.assertEqual(self.table["AAAA"] - self.table["CCCC"], 746)
        self.assertEqual(self.table["GGGG"] - self.table["UUUU"], -746)
        self.assertEqual(
            max(self.table.values()) - min(self.table.values()), 792)

    def test_v07a_algebraic_subtraction_swap(self):
        """V07a: exchanging the operands negates the score."""
        for key_a, key_b in (("CAAA", "GGCU"), ("GCAA", "GCUU"),
                             ("GGUC", "GGUC")):
            forward = self.table[key_a] - self.table[key_b]
            swapped = self.table[key_b] - self.table[key_a]
            self.assertEqual(swapped, -forward)
            if forward == 0:
                self.assertEqual(m6.classify(swapped), "neutral")
            else:
                self.assertNotEqual(m6.classify(swapped), m6.classify(forward))

    def test_v07b_physical_duplex_strand_swap(self):
        """V07b: relabelling which strand of the same duplex is the guide.

        This is a statement about the modeled duplex rather than operand order:
        designating the passenger as guide exchanges helix END A with END B.
        """
        guide = m6.guide_for_window("ATGGCTTGCTACACTCTTTTG")
        end_a, end_b = m6.terminal_keys(guide)
        original = self.table[end_a] - self.table[end_b]
        # The same 19 bp duplex, now read with the passenger as intended guide:
        # its 5' terminal 4 bp is END B and the opposite end is END A.
        swapped = self.table[end_b] - self.table[end_a]
        self.assertEqual(swapped, -original)
        self.assertEqual(m6.classify(original), "favorable")
        self.assertEqual(m6.classify(swapped), "unfavorable")

    def test_v07c_opposite_oriented_query_is_not_a_strand_swap(self):
        """V07c: the plus-orientation query builds a different duplex."""
        window = "ATGGCTTGCTACACTCTTTTG"
        minus = self.measure(window)["asym_score_hundredths"]
        # Treating the sense window itself as a guide is a different duplex,
        # because the strands are staggered by 2 nt.
        plus = m6.measure_guide(
            window, self.table, self.sha)["asym_score_hundredths"]
        self.assertNotEqual(plus, -minus)
        differing = 0
        random.seed(11)
        for _ in range(300):
            sample = "".join(
                random.choice("ACGT") for _ in range(m6.WINDOW_LENGTH))
            a = self.measure(sample)["asym_score_hundredths"]
            b = m6.measure_guide(
                sample, self.table, self.sha)["asym_score_hundredths"]
            if b != -a:
                differing += 1
        self.assertGreater(differing, 200)

    def test_v08_coordinate_forms_agree(self):
        """V08: guide-coordinate and sense-window forms are identical."""
        random.seed(5)
        sequence = "".join(random.choice("ACGT") for _ in range(400))
        for start in range(len(sequence) - m6.WINDOW_LENGTH + 1):
            window = sequence[start:start + m6.WINDOW_LENGTH]
            guide = m6.guide_for_window(window)
            end_a, end_b = m6.terminal_keys(guide)
            # Sense-window form: END A is revcomp(w[18..21]), END B is w[3..6].
            self.assertEqual(
                end_a,
                m6.reverse_complement_rna(window[17:21].replace("T", "U")))
            self.assertEqual(end_b, window[2:6].replace("T", "U"))

    def test_v09_guide_tail_positions_are_unused(self):
        """V09: guide positions 20-21 cannot change the score."""
        base = "CAAAAGAGTGTAGCAAGCC"
        first = m6.measure_guide(base + "AT", self.table, self.sha)
        second = m6.measure_guide(base + "GC", self.table, self.sha)
        self.assertEqual(first["asym_score_hundredths"],
                         second["asym_score_hundredths"])
        self.assertNotEqual(first["measurement_id"], second["measurement_id"])


class AuthorOracleTest(unittest.TestCase):
    """Agreement with the author implementation (V10, V22)."""

    @classmethod
    def setUpClass(cls):
        """Load the vendored table and its export checksum once.

        Args:
            None.

        Returns:
            None.
        """
        cls.table, identity = m6.load_energy_table(
            ENERGY_EXPORT, ENERGY_SOURCE)
        cls.sha = identity["energy_table_export_sha256"]

    def test_v10_matches_dsrip_arithmetic(self):
        """V10: scores equal dsRIP's own padded-index construction."""
        random.seed(17)
        sequence = "".join(random.choice("ACGT") for _ in range(300))
        expected = dsrip_reference_asym(sequence.replace("T", "U"), self.table)
        # dsRIP iterates its window index along the antisense strand, which
        # runs opposite to the sense-window order used here.
        windows = len(sequence) - m6.WINDOW_LENGTH + 1
        for index, (guide, score, _) in enumerate(expected):
            start = windows - 1 - index
            record = m6.measure_window(
                sequence[start:start + m6.WINDOW_LENGTH], self.table, self.sha)
            self.assertEqual(record["guide_sequence"].replace("T", "U"), guide)
            self.assertEqual(record["asym_score_hundredths"], score)

    def test_v22_staggered_passenger_geometry(self):
        """V22: modeled passenger interval, core and overhang availability."""
        random.seed(23)
        sequence = "".join(random.choice("ACGT") for _ in range(120))
        expected = dsrip_reference_asym(sequence.replace("T", "U"), self.table)
        windows = len(sequence) - m6.WINDOW_LENGTH + 1
        padded_windows = []
        for index, (_, _, passenger_3_5) in enumerate(expected):
            start = windows - 1 - index
            geometry = m6.modeled_passenger_interval(start, len(sequence))
            low, high = geometry["modeled_passenger_interval_1based_inclusive"]
            self.assertEqual(low, start + 3)
            self.assertEqual(high, start + 23)
            span = "".join(
                sequence[i] if 0 <= i < len(sequence) else "X"
                for i in range(start + 2, start + 23))
            self.assertEqual(passenger_3_5[::-1], span.replace("T", "U"))
            if "X" in passenger_3_5:
                padded_windows.append(start)
                self.assertFalse(
                    geometry["modeled_passenger_overhang_available"])
            else:
                self.assertTrue(
                    geometry["modeled_passenger_overhang_available"])
        # Only the final two window positions lack a real passenger overhang.
        self.assertEqual(sorted(padded_windows), [windows - 2, windows - 1])


class AggregationAndFailureTest(unittest.TestCase):
    """Aggregation, withholding and the failure taxonomy (V16-V18, V20)."""

    @classmethod
    def setUpClass(cls):
        """Load the vendored table and its export checksum once.

        Args:
            None.

        Returns:
            None.
        """
        cls.table, identity = m6.load_energy_table(
            ENERGY_EXPORT, ENERGY_SOURCE)
        cls.sha = identity["energy_table_export_sha256"]

    def test_v16_denominator_is_fixed(self):
        """V16: the denominator stays 180 under duplicates and failures."""
        record = m6.aggregate_region(["favorable"] * 70 + ["unfavorable"] * 107
                                     + ["neutral"] * 3)
        self.assertEqual(record["a_denominator"], 180)
        self.assertEqual(record["a_exact"], "70/180")
        self.assertEqual(record["completion_status"], "complete")
        with_failure = m6.aggregate_region(["favorable"] * 179 + [None])
        self.assertEqual(with_failure["a_denominator"], 180)
        self.assertEqual(with_failure["unavailable"], 1)
        self.assertIsNone(with_failure["a_numerator"])
        with self.assertRaises(ValueError):
            m6.aggregate_region(["favorable"] * 179)

    def test_v17_withholding_and_propagation(self):
        """V17: unavailable and unresolved-orientation both withhold A."""
        unavailable = m6.aggregate_region(["favorable"] * 179 + [None])
        self.assertEqual(unavailable["completion_status"],
                         "withheld_unavailable_measurement")
        self.assertIsNone(unavailable["a_exact"])
        # Neutral and unfavorable are successful measurements, not failures.
        complete = m6.aggregate_region(["neutral"] * 180)
        self.assertEqual(complete["completion_status"], "complete")
        self.assertEqual(complete["a_exact"], "0/180")
        self.assertEqual(complete["unavailable"], 0)
        orientation = m6.aggregate_region(
            ["favorable"] * 180, guide_orientation_resolved=False)
        self.assertEqual(orientation["completion_status"],
                         "withheld_guide_orientation")
        self.assertIsNone(orientation["a_numerator"])
        self.assertEqual(orientation["favorable"], 180)

    def test_v18_failure_taxonomy(self):
        """V18: each cause produces its typed, non-zero failure record."""
        invalid = m6.measure_window(
            "ATGGCNTGCTACACTCTTTTG", self.table, self.sha)
        self.assertEqual(invalid["status"], "invalid_symbol")
        self.assertIsNone(invalid["classification"])
        self.assertNotIn("asym_score_hundredths", invalid)
        short = m6.measure_window("ATGGC", self.table, self.sha)
        self.assertEqual(short["status"], "missing_context")
        absent = m6.measure_guide(
            "CAAAAGAGTGTAGCAAGCCAT",
            {k: v for k, v in self.table.items() if k != "CAAA"}, self.sha)
        self.assertEqual(absent["status"], "unsupported_table_entry")
        self.assertIn("CAAA", absent["cause"])
        with self.assertRaises(m6.EnergyTableError):
            m6.load_energy_table(os.path.join(REPO, "parametros_energia", "absent.tsv"))

    def test_v20_exact_rational_scores(self):
        """V20: A is an exact rational, never a rounded float."""
        record = m6.aggregate_region(
            ["favorable"] * 60 + ["unfavorable"] * 120)
        self.assertEqual(record["a_exact"], "60/180")
        self.assertEqual(Fraction(60, 180), Fraction(1, 3))
        # Two distinct numerators must not collapse under display rounding.
        first = m6.aggregate_region(["favorable"] * 60 + ["unfavorable"] * 120)
        second = m6.aggregate_region(
            ["favorable"] * 61 + ["unfavorable"] * 119)
        self.assertNotEqual(first["a_exact"], second["a_exact"])

    def test_format_hundredths_round_trip(self):
        """Exact decimal formatting covers signs and sub-unit magnitudes."""
        self.assertEqual(m6.format_hundredths(0), "0.00")
        self.assertEqual(m6.format_hundredths(-6), "-0.06")
        self.assertEqual(m6.format_hundredths(458), "4.58")
        self.assertEqual(m6.format_hundredths(-592), "-5.92")
