"""Calculate exact A, P and weighted regional scores."""

import hashlib
import json
import os
import re
from fractions import Fraction

REGION_LENGTH = 200
WINDOW_LENGTH = 21
POSITIONS_PER_REGION = REGION_LENGTH - WINDOW_LENGTH + 1
DENOMINATOR = POSITIONS_PER_REGION

BASELINE_LAMBDA = Fraction(1)
LAMBDA_VALUES = (Fraction(1, 2), Fraction(1), Fraction(2))
P_INTERFACE_VERSION = "p-interface-v0-provisional"
STAGE_ID = "m8-ranking"
IMPLEMENTATION_VERSION = "m8-v1"

FAVORABLE = "favorable"
A_CLASSIFICATIONS = (FAVORABLE, "unfavorable", "neutral")

REGION_ID = re.compile(r"^region:(?P<parent>[^:]+@[0-9a-f]{64}):"
                       r"(?P<start>\d+):(?P<end>\d+)$")

WITHHELD_CAUSES = (
    "a_unavailable",
    "p_unavailable",
    "a_withheld_upstream",
    "p_withheld_upstream",
)

SETTINGS = {
    "stage": STAGE_ID,
    "implementation_version": IMPLEMENTATION_VERSION,
    "region_length": REGION_LENGTH,
    "window_length": WINDOW_LENGTH,
    "positions_per_region": POSITIONS_PER_REGION,
    "denominator": DENOMINATOR,
    "lambda_values": ["1/2", "1", "2"],
    "baseline_lambda": "1",
    "arithmetic": "exact rational (fractions.Fraction)",
    "rounding_policy": "display only; never read back for comparison",
    "p_interface_version": P_INTERFACE_VERSION,
    "position_contribution": "each position contributes at most once to A "
                             "and at most once to P",
    "m7_descriptor_weight": 0,
}


class M8InputError(Exception):
    """Raised when an input artifact cannot be read or parsed at all."""


def lambda_key(value):
    """Render a lambda weight as its canonical output key.

    Args:
        value: Lambda weight as a Fraction.

    Returns:
        "0.5", "1" or "2" for the supported weights, else the exact string.
    """
    if value == Fraction(1, 2):
        return "0.5"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def exact_string(value):
    """Render a Fraction as a canonical, sign-preserving exact string.

    Args:
        value: Fraction to render.

    Returns:
        "numerator/denominator", with denominator 1 kept explicit so the
        string is never mistaken for a rounded decimal.
    """
    return f"{value.numerator}/{value.denominator}"


def display_string(value, places=6):
    """Render a Fraction as a fixed-width decimal for display only.

    Args:
        value: Fraction to render.
        places: Decimal places to emit.

    Returns:
        Decimal string. Never used for comparison, ordering or ties.
    """
    scale = 10 ** places
    scaled = value * scale
    whole = int(scaled)
    # Round half away from zero without going through binary floating point.
    remainder = scaled - whole
    if remainder >= Fraction(1, 2):
        whole += 1
    elif remainder <= Fraction(-1, 2):
        whole -= 1
    sign = "-" if whole < 0 else ""
    whole = abs(whole)
    return f"{sign}{whole // scale}.{whole % scale:0{places}d}"


def activity(favorable):
    """Compute the exact activity component A.

    Args:
        favorable: Count of favorable positions in the region.

    Returns:
        A as a Fraction over the fixed denominator.
    """
    return Fraction(favorable, DENOMINATOR)


def penalty(flagged):
    """Compute the exact penalty component P.

    Args:
        flagged: Count of flagged positions in the region.

    Returns:
        P as a Fraction over the fixed denominator.
    """
    return Fraction(flagged, DENOMINATOR)


def score(favorable, flagged, weight):
    """Compute the exact combined score S_lambda.

    Args:
        favorable: Count of favorable positions.
        flagged: Count of flagged positions.
        weight: Lambda weight as a Fraction.

    Returns:
        S_lambda as an exact Fraction, which may be negative.
    """
    return activity(favorable) - weight * penalty(flagged)


def score_block(favorable, flagged, weights=LAMBDA_VALUES):
    """Build the per-weight score mapping for one region.

    Args:
        favorable: Count of favorable positions.
        flagged: Count of flagged positions.
        weights: Iterable of lambda weights as Fractions.

    Returns:
        Mapping from lambda key to exact components and a display decimal.
    """
    block = {}
    for weight in weights:
        value = score(favorable, flagged, weight)
        block[lambda_key(weight)] = {
            "exact": exact_string(value),
            "numerator": value.numerator,
            "denominator": value.denominator,
            "display": display_string(value),
        }
    return block


def sha256_file(path):
    """Fingerprint a consumed artifact.

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


def read_jsonl(path):
    """Read newline-delimited JSON records.

    Args:
        path: Path to the JSONL file.

    Returns:
        List of parsed records.

    Raises:
        M8InputError: If the file cannot be read or a line cannot be parsed.
    """
    rows = []
    try:
        with open(path, "r", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except ValueError as error:
                    raise M8InputError(
                        f"{path} line {number} is not JSON: "
                        f"{error}") from error
    except OSError as error:
        raise M8InputError(f"{path} unreadable: {error}") from error
    return rows


def read_json(path):
    """Read a JSON document.

    Args:
        path: Path to the JSON file.

    Returns:
        Parsed document.

    Raises:
        M8InputError: If the file cannot be read or parsed.
    """
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError) as error:
        raise M8InputError(f"{path} unreadable: {error}") from error


def parse_region_id(region_id):
    """Split a region identifier into its parent and zero-based interval.

    Args:
        region_id: Region identifier string.

    Returns:
        Tuple of parent id, start and end, or None when malformed.
    """
    match = REGION_ID.match(region_id or "")
    if match is None:
        return None
    return (match.group("parent"), int(match.group("start")),
            int(match.group("end")))


def window_starts_for_region(interval_1based):
    """List the one-based window starts a region contains.

    A region at 1-based [s, e] contains the 180 window positions whose starts
    run s..s+179, the last window ending exactly at e.

    Args:
        interval_1based: Two-element region interval, one-based inclusive.

    Returns:
        List of one-based window start coordinates.
    """
    start = interval_1based[0]
    return list(range(start, start + POSITIONS_PER_REGION))


class Findings:
    """Collects typed validation findings without aborting on the first one."""

    def __init__(self):
        """Initialise an empty finding list.

        Args:
            None.

        Returns:
            None.
        """
        self.items = []

    def add(self, kind, detail, **fields):
        """Record one typed finding.

        Args:
            kind: Failure class from the design's validation table.
            detail: Human-readable explanation.
            **fields: Additional identifying fields to retain.

        Returns:
            None.
        """
        entry = {"kind": kind, "detail": detail}
        entry.update(fields)
        self.items.append(entry)

    def kinds(self):
        """Summarise findings by class.

        Args:
            None.

        Returns:
            Mapping from failure class to count.
        """
        counts = {}
        for item in self.items:
            counts[item["kind"]] = counts.get(item["kind"], 0) + 1
        return counts

    def ordered(self):
        """Return findings in a deterministic order.

        Args:
            None.

        Returns:
            Findings sorted by class then serialised detail.
        """
        return sorted(self.items,
                      key=lambda item: (item["kind"],
                                        json.dumps(item, sort_keys=True)))


def load_activity_inputs(directory, findings, prefix="asymmetry"):
    """Load and validate the activity artifacts.

    Args:
        directory: Directory holding the activity stage outputs.
        findings: Findings collector.
        prefix: Basename prefix, "asymmetry" for real M6 output or a fixture
            prefix such as "synthetic_a".

    Returns:
        Dictionary with positions indexed by window id, region records indexed
        by region id, the stage summary, and input fingerprints.

    Raises:
        M8InputError: If a required artifact is missing or unparsable.
    """
    positions_path = os.path.join(directory, f"{prefix}_positions.jsonl")
    regions_path = os.path.join(directory, f"{prefix}_regions.jsonl")
    summary_path = os.path.join(directory, f"{prefix}.json")
    positions = read_jsonl(positions_path)
    regions = read_jsonl(regions_path)
    summary = read_json(summary_path)

    by_window = {}
    for row in positions:
        window_id = row.get("window_id")
        if window_id in by_window:
            findings.add("duplicate_position",
                         "window id repeated in activity positions",
                         window_id=window_id)
            continue
        classification = row.get("classification")
        if classification is not None and \
                classification not in A_CLASSIFICATIONS:
            findings.add("unknown_classification",
                         "classification outside the known vocabulary",
                         window_id=window_id,
                         classification=row["classification"])
        by_window[window_id] = row

    by_region = {}
    for row in regions:
        region_id = row.get("region_id")
        if region_id in by_region:
            findings.add("duplicate_region",
                         "region id repeated in activity regions",
                         region_id=region_id)
            continue
        if row.get("completion_status") is None:
            findings.add("completion_state_absent",
                         "activity region has no completion status",
                         region_id=region_id)
        by_region[region_id] = row

    return {
        "positions_by_window": by_window,
        "regions_by_id": by_region,
        "summary": summary,
        "fingerprints": {
            os.path.basename(positions_path): sha256_file(positions_path),
            os.path.basename(regions_path): sha256_file(regions_path),
            os.path.basename(summary_path): sha256_file(summary_path),
        },
        "method_identity": {
            "method_id": summary.get("settings", {}).get("method_id"),
            "energy_table": summary.get("energy_table_identity"),
            "completion_status": summary.get("completion_status"),
        },
    }


def load_penalty_inputs(directory, findings, prefix="p"):
    """Load and validate the provisional penalty artifacts.

    Args:
        directory: Directory holding the penalty inputs.
        findings: Findings collector.
        prefix: Basename prefix, "p" for production shape or a fixture prefix
            such as "synthetic_p".

    Returns:
        Dictionary with positions indexed by window id, the stage summary,
        fingerprints and the declared method/policy identity.

    Raises:
        M8InputError: If a required artifact is missing or unparsable.
    """
    positions_path = os.path.join(directory, f"{prefix}_positions.jsonl")
    summary_path = os.path.join(directory, f"{prefix}.json")
    positions = read_jsonl(positions_path)
    summary = read_json(summary_path)

    declared = summary.get("p_interface_version")
    if declared != P_INTERFACE_VERSION:
        findings.add("p_interface_version_mismatch",
                     "penalty input declares an unsupported interface version",
                     declared=declared, supported=P_INTERFACE_VERSION)

    by_window = {}
    for row in positions:
        window_id = row.get("window_id")
        if window_id in by_window:
            findings.add("duplicate_position",
                         "window id repeated in penalty positions",
                         window_id=window_id)
            continue
        status = row.get("status")
        if status is None:
            findings.add("completion_state_absent",
                         "penalty position has no status",
                         window_id=window_id)
        # A flag can only come from a completed measurement; the contrary
        # combination would silently inject a penalty with no evidence.
        if row.get("flagged") and status != "measured":
            findings.add("p_flag_contradiction",
                         "flagged position is not a measured result",
                         window_id=window_id, status=status)
        by_window[window_id] = row

    return {
        "positions_by_window": by_window,
        "summary": summary,
        "fingerprints": {
            os.path.basename(positions_path): sha256_file(positions_path),
            os.path.basename(summary_path): sha256_file(summary_path),
        },
        "method_identity": {
            "method_id": summary.get("method_id"),
            "policy_id": summary.get("policy_id"),
            "p_interface_version": declared,
            "evidence_class": summary.get("evidence_class"),
            "completion_status": summary.get("completion_status"),
        },
    }


def _upstream_region_status(record, stage):
    """Classify a producing stage's regional completion state.

    Args:
        record: Region record from the producing stage.
        stage: Either "a" or "p", used to name the cause.

    Returns:
        None when complete, otherwise the withheld cause string.
    """
    status = (record or {}).get("completion_status")
    if status == "complete":
        return None
    return f"{stage}_withheld_upstream"


def score_region(region_record, a_inputs, p_inputs, findings,
                 weights=LAMBDA_VALUES):
    """Score one region from per-position activity and penalty inputs.

    Both components are aggregated here so the at-most-once rule and the fixed
    denominator are enforced in one place, and the activity aggregation is
    cross-checked against the producing stage's own regional result.

    Args:
        region_record: Activity region record for this region.
        a_inputs: Loaded activity inputs.
        p_inputs: Loaded penalty inputs.
        findings: Findings collector.
        weights: Lambda weights as Fractions.

    Returns:
        Region scoring record, with status "scored" or "withheld".
    """
    region_id = region_record["region_id"]
    parent_id = region_record["parent_id"]
    interval = region_record["interval_1based_inclusive"]
    favorable = 0
    flagged = 0
    a_unavailable = 0
    p_unavailable = 0
    site_total = 0
    organisms = set()
    orientations = set()

    for start in window_starts_for_region(interval):
        window_id = (f"window:{parent_id}:{start - 1}:"
                     f"{start - 1 + WINDOW_LENGTH}")
        a_row = a_inputs["positions_by_window"].get(window_id)
        p_row = p_inputs["positions_by_window"].get(window_id)
        if a_row is None:
            findings.add("missing_position",
                         "activity position absent for a region window",
                         region_id=region_id, window_id=window_id)
            a_unavailable += 1
        elif a_row.get("status") != "measured" or \
                a_row.get("classification") is None:
            a_unavailable += 1
        elif a_row["classification"] == FAVORABLE:
            favorable += 1

        if p_row is None:
            findings.add("missing_position",
                         "penalty position absent for a region window",
                         region_id=region_id, window_id=window_id)
            p_unavailable += 1
        elif p_row.get("status") != "measured":
            p_unavailable += 1
        elif p_row.get("flagged"):
            # One flagged position contributes exactly one unit regardless of
            # how many sites, organisms or records produced the flag.
            flagged += 1
            site_total += int(p_row.get("qualifying_site_count") or 0)
            organisms.update(p_row.get("organisms") or [])
            orientations.update(p_row.get("orientations") or [])
        elif p_row.get("qualifying_site_count"):
            site_total += int(p_row.get("qualifying_site_count") or 0)

    # The activity adapter must reproduce M6's own aggregation exactly.
    reported = region_record.get("a_numerator")
    if not a_unavailable and reported is not None and reported != favorable:
        findings.add("a_adapter_disagreement",
                     "re-aggregated activity count differs from the "
                     "producing stage's regional value",
                     region_id=region_id, reaggregated=favorable,
                     reported=reported)

    causes = []
    upstream_a = _upstream_region_status(region_record, "a")
    if upstream_a:
        causes.append(upstream_a)
    if a_unavailable:
        causes.append("a_unavailable")
    if p_unavailable:
        causes.append("p_unavailable")

    record = {
        "region_id": region_id,
        "candidate_id": region_record.get("candidate_id"),
        "parent_id": parent_id,
        "interval_1based_inclusive": interval,
        "position_count": POSITIONS_PER_REGION,
        "denominator": DENOMINATOR,
        "favorable": favorable,
        "flagged": flagged,
        "a_unavailable_positions": a_unavailable,
        "p_unavailable_positions": p_unavailable,
        "a_exact": exact_string(activity(favorable)),
        "p_exact": exact_string(penalty(flagged)),
        "scores": score_block(favorable, flagged, weights),
        "penalty_evidence": {
            "qualifying_site_total": site_total,
            "organisms": sorted(organisms),
            "orientations": sorted(orientations),
            "note": "evidence only; never summed into P",
        },
        "status": "withheld" if causes else "scored",
        "withheld_causes": sorted(set(causes)) or None,
    }
    return record


def validate_region_population(region_records, parent_lengths, findings):
    """Check each parent's region population against its sequence length.

    Args:
        region_records: Iterable of activity region records.
        parent_lengths: Mapping from parent id to sequence length.
        findings: Findings collector.

    Returns:
        None.
    """
    seen = {}
    for record in region_records:
        seen.setdefault(record["parent_id"], 0)
        seen[record["parent_id"]] += 1
    for parent_id, count in sorted(seen.items()):
        length = parent_lengths.get(parent_id)
        if length is None:
            findings.add("unexpected_region",
                         "regions present for an unknown parent",
                         parent_id=parent_id)
            continue
        expected = max(0, length - REGION_LENGTH + 1)
        if count != expected:
            findings.add("region_population_mismatch",
                         "region count differs from the enumerated population",
                         parent_id=parent_id, observed=count,
                         expected=expected)


def validate_region_identity(record, findings):
    """Check one activity region record's internal identifier consistency.

    Args:
        record: Activity region record.
        findings: Findings collector.

    Returns:
        None.
    """
    region_id = record.get("region_id")
    parsed = parse_region_id(region_id)
    if parsed is None:
        findings.add("identifier_inconsistent",
                     "region id does not match the expected form",
                     region_id=region_id)
        return
    parent, start, end = parsed
    if parent != record.get("parent_id"):
        findings.add("identifier_inconsistent",
                     "region id parent differs from the parent field",
                     region_id=region_id, parent_id=record.get("parent_id"))
    interval = record.get("interval_1based_inclusive") or [None, None]
    if interval != [start + 1, end]:
        findings.add("identifier_inconsistent",
                     "region id interval differs from the interval field",
                     region_id=region_id, interval=interval)
    if end - start != REGION_LENGTH:
        findings.add("identifier_inconsistent",
                     "region id span is not the configured region length",
                     region_id=region_id, span=end - start)
    candidate = parent.split("@", 1)[0]
    if candidate != record.get("candidate_id"):
        findings.add("identifier_inconsistent",
                     "parent id candidate differs from the candidate field",
                     region_id=region_id,
                     candidate_id=record.get("candidate_id"))
    if record.get("position_count") not in (None, POSITIONS_PER_REGION):
        findings.add("position_count_mismatch",
                     "region declares a position count other than the fixed "
                     "value",
                     region_id=region_id,
                     position_count=record.get("position_count"))


def cross_check_regions(a_regions, p_positions_by_window, findings):
    """Detect penalty positions that belong to no enumerated region.

    Args:
        a_regions: Mapping of activity region id to record.
        p_positions_by_window: Mapping of penalty window id to record.
        findings: Findings collector.

    Returns:
        None.
    """
    expected = set()
    for record in a_regions.values():
        parent_id = record["parent_id"]
        for start in window_starts_for_region(
                record["interval_1based_inclusive"]):
            expected.add(f"window:{parent_id}:{start - 1}:"
                         f"{start - 1 + WINDOW_LENGTH}")
    for window_id in sorted(set(p_positions_by_window) - expected):
        findings.add("unexpected_position",
                     "penalty position is outside every enumerated region",
                     window_id=window_id)


def gene_eligibility(region_records):
    """Decide per gene whether a final ranking may be produced.

    A single withheld region withholds the whole gene; the gene's regions are
    still reported so the decision is auditable.

    Args:
        region_records: Iterable of scored region records.

    Returns:
        Mapping from candidate id to an eligibility summary.
    """
    genes = {}
    for record in region_records:
        candidate_id = record["candidate_id"]
        entry = genes.setdefault(candidate_id, {
            "candidate_id": candidate_id,
            "regions": 0,
            "regions_scored": 0,
            "regions_withheld": 0,
            "withheld_causes": set(),
            "ranking_withheld": False,
        })
        entry["regions"] += 1
        if record["status"] == "scored":
            entry["regions_scored"] += 1
        else:
            entry["regions_withheld"] += 1
            entry["withheld_causes"].update(record["withheld_causes"] or [])
    for entry in genes.values():
        entry["ranking_withheld"] = entry["regions_withheld"] > 0
        entry["withheld_causes"] = sorted(entry["withheld_causes"])
        entry["ranking_eligible"] = not entry["ranking_withheld"]
    return genes


def propagate_shared_causes(genes, shared_causes):
    """Mark every gene withheld by a cause that is not region specific.

    Args:
        genes: Eligibility mapping from `gene_eligibility`.
        shared_causes: Iterable of run-level cause strings.

    Returns:
        Sorted list of candidate ids withheld by a shared cause.
    """
    shared = sorted(set(shared_causes))
    if not shared:
        return []
    for entry in genes.values():
        entry["ranking_withheld"] = True
        entry["ranking_eligible"] = False
        entry["withheld_causes"] = sorted(
            set(entry["withheld_causes"]) | set(shared))
        entry["shared_cause"] = True
    return sorted(genes)


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
