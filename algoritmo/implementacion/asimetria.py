"""Model terminal duplex-energy asymmetry and aggregate regional A."""

import hashlib
import itertools
import json
from fractions import Fraction

REGION_LENGTH = 200
WINDOW_LENGTH = 21
POSITIONS_PER_REGION = REGION_LENGTH - WINDOW_LENGTH + 1
DUPLEX_CORE_LENGTH = 19
OVERHANG_LENGTH = 2
TERMINAL_SEGMENT_LENGTH = 4
METHOD_ID = "dsrip-thermoasym-v1"
EXPECTED_TABLE_ENTRIES = 256

DNA_ALPHABET = "ACGT"
RNA_ALPHABET = "ACGU"
_DNA_COMPLEMENT = str.maketrans(DNA_ALPHABET, "TGCA")
_RNA_COMPLEMENT = str.maketrans(RNA_ALPHABET, "UGCA")

CLASSIFICATIONS = ("favorable", "unfavorable", "neutral")
FAILURE_STATUSES = (
    "invalid_symbol",
    "missing_context",
    "unsupported_table_entry",
    "calculation_error",
)

SETTINGS = {
    "method_id": METHOD_ID,
    "region_length": REGION_LENGTH,
    "window_length": WINDOW_LENGTH,
    "positions_per_region": POSITIONS_PER_REGION,
    "duplex_core_length": DUPLEX_CORE_LENGTH,
    "overhang_length": OVERHANG_LENGTH,
    "terminal_segment_length": TERMINAL_SEGMENT_LENGTH,
    "guide_orientation": "-",
    "guide_definition": "reverse complement of the supplied sense window",
    "energy_units": "kcal/mol",
    "energy_model_state": "dG37, 1 M NaCl optical-melting basis",
    "arithmetic": "exact integer hundredths of kcal/mol",
    "classification_threshold": "favorable > 0, unfavorable < 0, neutral == 0",
    "neutral_tolerance_band": None,
    "a_denominator": POSITIONS_PER_REGION,
}


class EnergyTableError(Exception):
    """Raised when the vendored energy table fails its load-time validation."""


def _sha256_file(path):
    """Hash a file's bytes.

    Args:
        path: Path to an existing file.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def _md5_file(path):
    """Hash a file's bytes with MD5 to compare against the publisher manifest.

    Args:
        path: Path to an existing file.

    Returns:
        Lowercase hexadecimal MD5 digest.
    """
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement_dna(sequence):
    """Reverse-complement a normalized A/C/G/T sequence.

    Args:
        sequence: Uppercase A/C/G/T sequence.

    Returns:
        Reverse complement, leaving the input unchanged.

    Raises:
        ValueError: If unsupported symbols are present.
    """
    if not isinstance(sequence, str) or set(sequence) - set(DNA_ALPHABET):
        raise ValueError("Reverse complement requires normalized A/C/G/T")
    return sequence.translate(_DNA_COMPLEMENT)[::-1]


def reverse_complement_rna(sequence):
    """Reverse-complement a normalized A/C/G/U sequence.

    Args:
        sequence: Uppercase A/C/G/U sequence.

    Returns:
        Reverse complement, leaving the input unchanged.

    Raises:
        ValueError: If unsupported symbols are present.
    """
    if not isinstance(sequence, str) or set(sequence) - set(RNA_ALPHABET):
        raise ValueError("Reverse complement requires normalized A/C/G/U")
    return sequence.translate(_RNA_COMPLEMENT)[::-1]


def parse_energy_text(text):
    """Parse the vendored two-column energy export into integer hundredths.

    Values are kept as integer hundredths of a kcal/mol so that differences are
    exact and neutrality is decidable without a tolerance.

    Args:
        text: Contents of the TSV export, including its header line.

    Returns:
        Mapping from 4-mer key to signed integer hundredths.

    Raises:
        EnergyTableError: For a malformed header, a bad key, a value
            that is not exactly two decimal places, or a duplicate key.
    """
    lines = text.splitlines()
    if not lines or lines[0].split("\t") != ["key", "delta_g_37_kcal_per_mol"]:
        raise EnergyTableError(
            "Energy export header is not the expected two columns")
    table = {}
    for number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise EnergyTableError(
                f"Energy export line {number} is not two fields")
        key, raw = fields
        if len(key) != TERMINAL_SEGMENT_LENGTH or set(key) - set(RNA_ALPHABET):
            raise EnergyTableError(
                f"Energy export line {number} has a bad key: {key}")
        if key in table:
            raise EnergyTableError(f"Energy export repeats key {key}")
        sign = -1 if raw.startswith("-") else 1
        body = raw[1:] if raw[0] in "+-" else raw
        whole, _, fraction = body.partition(".")
        if not whole.isdigit() or len(fraction) != 2 or not fraction.isdigit():
            raise EnergyTableError(
                f"Energy export line {number} is not two decimals: {raw}")
        table[key] = sign * (int(whole) * 100 + int(fraction))
    return table


def validate_energy_table(table):
    """Check the structural invariants the specification requires.

    Args:
        table: Mapping from 4-mer key to integer hundredths.

    Returns:
        None.

    Raises:
        EnergyTableError: If the table is incomplete or not
            reverse-complement symmetric.
    """
    if len(table) != EXPECTED_TABLE_ENTRIES:
        raise EnergyTableError(
            f"Energy table has {len(table)} entries, expected "
            f"{EXPECTED_TABLE_ENTRIES}")
    for combination in itertools.product(
            RNA_ALPHABET, repeat=TERMINAL_SEGMENT_LENGTH):
        key = "".join(combination)
        if key not in table:
            raise EnergyTableError(f"Energy table is missing key {key}")
    # A duplex read from either strand is one physical object, so the
    # published values must agree; a violation would mean the file is not
    # the dsRIP table.
    for key, value in table.items():
        mirror = reverse_complement_rna(key)
        if table[mirror] != value:
            raise EnergyTableError(
                f"Energy table is not reverse-complement symmetric at {key}")


def load_energy_table(export_path, source_path=None):
    """Load and validate the vendored table, recording its identity.

    Args:
        export_path: Path to the TSV export that the calculation parses.
        source_path: Optional path to the publisher's original binary table,
            hashed for provenance when present.

    Returns:
        Tuple of the parsed table and a checksum/identity dictionary.

    Raises:
        EnergyTableError: If the file cannot be read or fails validation.
    """
    try:
        with open(export_path, "r", encoding="ascii") as handle:
            text = handle.read()
    except OSError as error:
        raise EnergyTableError(f"Energy export unreadable: {error}") from error
    table = parse_energy_text(text)
    validate_energy_table(table)
    identity = {
        "method_id": METHOD_ID,
        "energy_table_export_sha256": _sha256_file(export_path),
        "energy_table_entries": len(table),
    }
    if source_path is not None:
        identity["energy_table_source_md5"] = _md5_file(source_path)
        identity["energy_table_source_sha256"] = _sha256_file(source_path)
    return table, identity


def format_hundredths(value):
    """Render signed integer hundredths as an exact two-decimal string.

    Args:
        value: Signed integer hundredths of a kcal/mol.

    Returns:
        Decimal string such as "-0.06", never a binary-float repr.
    """
    sign = "-" if value < 0 else ""
    magnitude = abs(value)
    return f"{sign}{magnitude // 100}.{magnitude % 100:02d}"


def guide_for_window(window):
    """Build the intended guide strand for a sense window.

    Args:
        window: 21 nt sense window in supplied parent orientation, A/C/G/T.

    Returns:
        The guide strand 5'->3', which is the reverse complement of the window.

    Raises:
        ValueError: If the window is the wrong length or has bad symbols.
    """
    if len(window) != WINDOW_LENGTH:
        raise ValueError(f"Guide construction requires {WINDOW_LENGTH} nt")
    return reverse_complement_dna(window)


def terminal_keys(guide):
    """Return the two helix-end lookup keys for a guide strand.

    END A is the guide's own 5'-terminal 4 bp. END B is the opposite helix end,
    which is the passenger's 5'-terminal 4 bp; it is read from the guide as the
    reverse complement of guide positions 16-19. Both lie inside the 19 bp
    core, so neither 2 nt 3' overhang is consulted.

    Args:
        guide: 21 nt guide strand 5'->3', A/C/G/T or A/C/G/U.

    Returns:
        Tuple of the END A and END B keys as RNA 4-mers.
    """
    rna = guide.replace("T", "U")
    end_a = rna[0:TERMINAL_SEGMENT_LENGTH]
    core_offset = DUPLEX_CORE_LENGTH - TERMINAL_SEGMENT_LENGTH
    end_b = reverse_complement_rna(rna[core_offset:DUPLEX_CORE_LENGTH])
    return end_a, end_b


def modeled_passenger_interval(window_start_0based, parent_length):
    """Locate the modeled passenger strand in sense-parent coordinates.

    The duplex strands are staggered by 2 nt, so the passenger occupies
    s[P+2 .. P+22] for a window starting at P. Its 3' overhang falls outside
    the parent at the final two window positions, where dsRIP substitutes
    "XX"; the overhang is excluded from the calculation either way.

    Args:
        window_start_0based: Zero-based start of the 21 nt sense window.
        parent_length: Length of the parent sequence.

    Returns:
        Dictionary with the one-based inclusive interval, the paired-core
        interval, and whether the 3' overhang is available in the parent.
    """
    start = window_start_0based + OVERHANG_LENGTH
    core_end = start + DUPLEX_CORE_LENGTH
    overhang_end = start + WINDOW_LENGTH
    return {
        "modeled_passenger_interval_1based_inclusive":
            [start + 1, overhang_end],
        "modeled_passenger_core_1based_inclusive": [start + 1, core_end],
        "modeled_passenger_overhang_available": overhang_end <= parent_length,
    }


def classify(score_hundredths):
    """Classify an exact asymmetry score under the approved zero threshold.

    Args:
        score_hundredths: Signed integer hundredths of a kcal/mol.

    Returns:
        One of "favorable", "unfavorable" or "neutral".
    """
    if score_hundredths > 0:
        return "favorable"
    if score_hundredths < 0:
        return "unfavorable"
    return "neutral"


def measurement_identity(guide, table_sha256, method_id=METHOD_ID):
    """Derive the reusable measurement identifier.

    Every input that can change the result is hashed: the method, the energy
    table actually parsed, and the exact oriented guide sequence. Fields are
    separated so components cannot be ambiguously concatenated.

    Args:
        guide: Exact oriented 21 nt guide sequence.
        table_sha256: SHA-256 of the parsed energy export.
        method_id: Method identity string.

    Returns:
        Identifier of the form "asym:sha256:<digest>".
    """
    payload = "\x1f".join((method_id, table_sha256, guide)).encode("ascii")
    return "asym:sha256:" + hashlib.sha256(payload).hexdigest()


def measure_guide(guide, table, table_sha256, method_id=METHOD_ID):
    """Measure one guide strand, returning either a measurement or a failure.

    Args:
        guide: Candidate 21 nt guide strand 5'->3'.
        table: Mapping from 4-mer key to integer hundredths.
        table_sha256: SHA-256 of the parsed energy export.
        method_id: Method identity string.

    Returns:
        Record whose "status" is "measured" or one of the failure statuses.
        A measured record carries both keys, both energies, the exact score and
        the classification; a failure carries a cause and no score.
    """
    digest = None
    if isinstance(guide, str):
        digest = hashlib.sha256(guide.encode("ascii")).hexdigest()
    record = {
        "measurement_id": None,
        "method_id": method_id,
        "energy_table_sha256": table_sha256,
        "guide_sequence": guide,
        "guide_sequence_sha256": digest,
        "status": "measured",
        "classification": None,
    }
    if not isinstance(guide, str) or len(guide) != WINDOW_LENGTH:
        record["status"] = "missing_context"
        record["cause"] = f"guide is not {WINDOW_LENGTH} nt"
        return record
    unsupported = sorted(set(guide) - set(DNA_ALPHABET) - set(RNA_ALPHABET))
    if unsupported:
        record["status"] = "invalid_symbol"
        record["cause"] = "unsupported symbols: " + ",".join(unsupported)
        return record
    record["measurement_id"] = measurement_identity(
        guide, table_sha256, method_id)
    try:
        end_a, end_b = terminal_keys(guide)
    except ValueError as error:
        record["status"] = "calculation_error"
        record["cause"] = f"{type(error).__name__}: {error}"
        return record
    record["end_a_key"] = end_a
    record["end_b_key"] = end_b
    if end_a not in table or end_b not in table:
        record["status"] = "unsupported_table_entry"
        record["cause"] = "absent key(s): " + ",".join(
            k for k in (end_a, end_b) if k not in table)
        return record
    energy_a = table[end_a]
    energy_b = table[end_b]
    score = energy_a - energy_b
    record["end_a_energy_kcal_mol"] = format_hundredths(energy_a)
    record["end_b_energy_kcal_mol"] = format_hundredths(energy_b)
    record["asym_score_kcal_mol"] = format_hundredths(score)
    record["asym_score_hundredths"] = score
    record["classification"] = classify(score)
    return record


def measure_window(window, table, table_sha256, method_id=METHOD_ID):
    """Measure the intended guide derived from a sense window.

    Args:
        window: 21 nt sense window in supplied parent orientation.
        table: Mapping from 4-mer key to integer hundredths.
        table_sha256: SHA-256 of the parsed energy export.
        method_id: Method identity string.

    Returns:
        Measurement or failure record as returned by `measure_guide`.
    """
    if not isinstance(window, str) or len(window) != WINDOW_LENGTH:
        return {
            "measurement_id": None, "method_id": method_id,
            "energy_table_sha256": table_sha256, "guide_sequence": None,
            "guide_sequence_sha256": None, "status": "missing_context",
            "classification": None,
            "cause": f"window is not {WINDOW_LENGTH} nt",
        }
    unsupported = sorted(set(window) - set(DNA_ALPHABET))
    if unsupported:
        return {
            "measurement_id": None, "method_id": method_id,
            "energy_table_sha256": table_sha256, "guide_sequence": None,
            "guide_sequence_sha256": None, "status": "invalid_symbol",
            "classification": None,
            "cause": "unsupported symbols: " + ",".join(unsupported),
        }
    return measure_guide(
        guide_for_window(window), table, table_sha256, method_id)


def aggregate_region(classifications, guide_orientation_resolved=True):
    """Aggregate per-position classifications into a regional A record.

    The denominator is fixed at 180 and is never reduced. A is withheld when
    any required position is unavailable, and separately when the candidate's
    intended-guide mapping is unresolved, even if every measurement
    succeeded.

    Args:
        classifications: Sequence of 180 values, each a classification name or
            None for an unavailable measurement.
        guide_orientation_resolved: Whether curation evidence established the
            parent as the intended target's sense sequence.

    Returns:
        Dictionary of counts, the exact rational A, and the completion status.

    Raises:
        ValueError: If the position count is not exactly 180.
    """
    if len(classifications) != POSITIONS_PER_REGION:
        raise ValueError(
            "Region aggregation requires exactly "
            f"{POSITIONS_PER_REGION} positions")
    counts = {name: 0 for name in CLASSIFICATIONS}
    unavailable = 0
    for value in classifications:
        if value is None:
            unavailable += 1
        elif value in counts:
            counts[value] += 1
        else:
            raise ValueError(f"Unknown classification: {value}")
    record = {
        "position_count": POSITIONS_PER_REGION,
        "favorable": counts["favorable"],
        "unfavorable": counts["unfavorable"],
        "neutral": counts["neutral"],
        "unavailable": unavailable,
        "a_denominator": POSITIONS_PER_REGION,
    }
    if unavailable:
        record["completion_status"] = "withheld_unavailable_measurement"
    elif not guide_orientation_resolved:
        record["completion_status"] = "withheld_guide_orientation"
    else:
        record["completion_status"] = "complete"
    if record["completion_status"] == "complete":
        exact = Fraction(counts["favorable"], POSITIONS_PER_REGION)
        record["a_numerator"] = counts["favorable"]
        record["a_exact"] = f"{counts['favorable']}/{POSITIONS_PER_REGION}"
        record["a_decimal_approximation"] = f"{float(exact):.6f}"
    else:
        record["a_numerator"] = None
        record["a_exact"] = None
        record["a_decimal_approximation"] = None
    return record


def write_jsonl(path, rows):
    """Write records as deterministic newline-delimited JSON.

    Args:
        path: Destination path.
        rows: Iterable of JSON-serializable dictionaries, already ordered.

    Returns:
        Number of records written.
    """
    count = 0
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(
                json.dumps(row, sort_keys=True, separators=(",", ":")))
            handle.write("\n")
            count += 1
    return count


def write_json(path, payload):
    """Write a JSON document with deterministic key ordering.

    Args:
        path: Destination path.
        payload: JSON-serializable document.

    Returns:
        None.
    """
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, sort_keys=True, indent=2)
        handle.write("\n")
