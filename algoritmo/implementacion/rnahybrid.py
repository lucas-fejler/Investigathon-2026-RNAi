"""Parse RNAhybrid evidence and classify intended and unintended interactions."""


from collections import Counter
import hashlib
import math

from archivos import identity

WC = {"AT", "TA", "CG", "GC"}
GU = {"GT", "TG"}


def fasta_bytes(rows):
    """Encode complete native FASTA inputs at safe line widths.

    Args:
        rows: Native alias and preserved sequence pairs.

    Returns:
        ASCII FASTA bytes with no cropping or sequence substitutions.
    """
    return "".join(
        ">" + name + "\n" + "".join(sequence[i:i + 80] + "\n"
                                     for i in range(0, len(sequence), 80))
        for name, sequence in rows).encode("ascii")


def parse_compact(raw, expected, profile):
    """Strictly parse native compact blocks with explicit filter decisions.

    Args:
        raw: Exact native stdout bytes.
        expected: Native alias-pair to exact source/length mapping.
        profile: Effective comparison filter/profile.

    Returns:
        All raw blocks, parsed fields, identities and filter decisions.
    """
    text = raw.decode("ascii")
    if raw and not raw.endswith(b"\n"):
        raise ValueError("Truncated native block framing")
    blocks = []
    for number, line in enumerate(text.splitlines(), 1):
        fields = line.split(":")
        if len(fields) != 11:
            raise ValueError("Malformed native block or diagnostic on stdout")
        target, m, query, n, energy, pvalue, position, *alignments = fields
        source = expected.get((query, target))
        if source is None:
            raise ValueError("Unexpected native target/query identifier")
        try:
            m, n, position = int(m), int(n), int(position)
            energy, pvalue = float(energy), float(pvalue)
        except ValueError as exc:
            raise ValueError("Malformed native numeric field") from exc
        if (m != source["target_length"] or n != source["query_length"]
                or hashlib.sha256(source["query_sequence"].encode()
                                  ).hexdigest() != source["query_sha256"]
                or hashlib.sha256(source["target_sequence"].encode()
                                  ).hexdigest() != source["reference_sha256"]
                or not math.isfinite(energy) or not math.isfinite(pvalue)
                or not 0 <= pvalue <= 1 or not 0 <= position <= m
                or len({len(row) for row in alignments}) != 1
                or any(set(row) - set("ACGUTNacgutn ") for row in alignments)):
            raise ValueError(
                "Native alignment/length/position is inconsistent")
        aligned = sum(len(row.replace(" ", ""))
                      for row in alignments[2:])
        if aligned > n:
            raise ValueError("Native query alignment exceeds query length")
        printed = "".join(a if a != " " else b for a, b in zip(
            alignments[2], alignments[3])).replace(" ", "").upper()
        if (any(a != " " and b != " " for a, b in zip(
                alignments[2], alignments[3]))
                or printed != source["query_sequence"][::-1].replace(
                    "T", "U")):
            raise ValueError("Native printed query disagrees with exact input")
        qualifies = (aligned >= profile["min_align_length"]
                     and energy <= profile["energy_cutoff"]
                     and pvalue <= profile["pvalue_cutoff"])
        site = {"query_id": source["query_id"],
                "query_sha256": source["query_sha256"],
                "reference_id": source["reference_id"],
                "reference_sha256": source["reference_sha256"],
                "native_position": position, "alignment_rows": alignments}
        blocks.append({"block_number": number, "raw_block": line,
                       **site, "mfe": energy, "pvalue": pvalue,
                       "aligned_query_nucleotides": aligned,
                       "qualifies": qualifies,
                       "filter_reason": "passed" if qualifies else
                       "fixture_filter_rejected", "site_id":
                       "native-site:sha256:" + identity(site)})
    return blocks


def parse_audit(raw, expected):
    """Verify actual parsed bytes and entry/return for all scheduled pairs.

    Args:
        raw: Exact audit file bytes.
        expected: Scheduled alias-pair identity and sequence ledger.

    Returns:
        Entered/finished ledgers and any precise accounting errors.
    """
    entered, finished, errors = [], [], []
    active = None
    for line in raw.decode("ascii").splitlines():
        fields = line.split("\t")
        if len(fields) == 7 and fields[0] == "BEGIN":
            _, query, target, n, m, qseq, tseq = fields
            pair = query, target
            source = expected.get(pair)
            if source is None or active is not None or pair in entered:
                errors.append("Unexpected/duplicate native evaluation entry")
                continue
            if (int(n) != source["query_length"]
                    or int(m) != source["target_length"]
                    or qseq != source["query_sequence"]
                    or tseq != source["target_sequence"]):
                errors.append("Native parsed bytes/lengths differ from input")
            entered.append(pair)
            active = pair
        elif len(fields) == 3 and fields[0] == "END":
            pair = tuple(fields[1:])
            if pair != active or pair in finished:
                errors.append("Unmatched native evaluation return")
            finished.append(pair)
            active = None
        else:
            errors.append("Malformed native audit event")
    if (set(entered) != set(expected) or set(finished) != set(expected)
            or active is not None):
        errors.append("Scheduled targets skipped or evaluation incomplete")
    return {"entered": entered, "finished": finished, "errors": errors}


def alignment_metrics(rows):
    """Count printed paired columns without interpreting native coordinates.

    Args:
        rows: Four native alignment rows, with spaces preserved.

    Returns:
        Watson-Crick/wobble counts, unpaired counts and exact-column run.
    """
    counts = Counter()
    run = longest = 0
    for target, query in zip(rows[1].upper().replace("U", "T"),
                             rows[2].upper().replace("U", "T")):
        pair = target + query
        if target != " " and query != " ":
            assert pair in WC | GU, (rows, pair)
            counts["wc_pairs" if pair in WC else "gu_pairs"] += 1
        run = run + 1 if pair in WC else 0
        longest = max(longest, run)
    return {"wc_pairs": counts["wc_pairs"], "gu_pairs": counts["gu_pairs"],
            "longest_exact_paired_tract": longest,
            "unpaired_target_characters": len(rows[0].replace(" ", "")),
            "unpaired_query_characters": len(rows[3].replace(" ", ""))}


def normalize(sequence):
    """Reproduce the pinned native spelling without altering source bytes.

    Args:
        sequence: Uppercase source DNA/IUPAC string.

    Returns:
        RNA spelling with non-ACGT symbols represented as N.
    """
    return "".join(b if b in "ACGU" else "N" for b in
                   sequence.replace("T", "U"))


def map_block(block, source, masked):
    """Validate one printed interaction against source and earlier masks.

    Args:
        block: Parsed native block with preserved four alignment rows.
        source: Exact expected query and target input ledger row.
        masked: Zero-based target positions masked by prior reported blocks.

    Returns:
        Verified source intervals, paired edges and masked-dangle count.
    """
    rows = block["alignment_rows"]
    if len(rows) != 4 or len({len(r) for r in rows}) != 1:
        raise ValueError("Unequal native alignment rows")
    target = normalize(source["target_sequence"])
    query = normalize(source["query_sequence"])
    start = block["native_position"] - 1
    if not 0 <= start < len(target):
        raise ValueError("Native position is outside the target")
    position, offset = start, len(query) - 1
    edges, masked_dangles = [], 0
    for a, b, c, d in zip(*rows):
        if (a != " " and b != " ") or (c != " " and d != " "):
            raise ValueError("Overlapping native row characters")
        t, q = (a if a != " " else b), (c if c != " " else d)
        if t != " ":
            if position >= len(target):
                raise ValueError("Printed target exceeds source")
            if position in masked:
                if a != "N" or b != " ":
                    raise ValueError("Masked target is not an unpaired N")
                masked_dangles += 1
            elif t != target[position]:
                raise ValueError("Printed target disagrees with source")
        if q != " " and (offset < 0 or q != query[offset]):
            raise ValueError("Printed query disagrees with source")
        if b != " " or c != " ":
            if b == " " or c == " ":
                raise ValueError("Unmatched paired-row base")
            pair = (b + c).replace("U", "T")
            if pair not in WC | GU:
                raise ValueError("Unsupported paired column")
            edges.append([offset, position, "WC" if pair in WC
                          else "GU"])
        position += t != " "
        offset -= q != " "
    if offset != -1 or not edges:
        raise ValueError("Incomplete query or empty duplex")
    core = [min(e[1] for e in edges), max(e[1] for e in edges) + 1]
    if len(edges) != sum(len(r.replace(" ", "")) for r in rows[2:3]):
        raise ValueError("Paired-row accounting differs")
    duplex = {"query_id": source["query_id"],
              "query_sha256": source["query_sha256"],
              "reference_id": source["reference_id"],
              "reference_sha256": source["reference_sha256"],
              "reference_orientation": "+", "paired_edges": edges}
    return {"printed_target_interval": [start, position],
            "core_target_interval": core, "paired_edges": edges,
            "masked_dangle_characters": masked_dangles,
            "duplex_id": "duplex:sha256:" + identity(duplex)}


def qualifies(block):
    """Apply the draft predicate to a retained energy-gated native block.

    Args:
        block: Strictly parsed block from an unchanged native -e -10 run.

    Returns:
        Whether the draft paired-base predicate qualifies this reported block.
    """
    metrics = alignment_metrics(block["alignment_rows"])
    return block["mfe"] <= -10 and (
        metrics["wc_pairs"] + metrics["gu_pairs"] >= 18
        or metrics["longest_exact_paired_tract"] >= 16)


def classify(candidate, reference_id, core, config):
    """Classify a verified footprint against the accepted candidate intervals.

    Args:
        candidate: Original candidate identifier, not a gene-label inference.
        reference_id: Exact frozen dataset/accession/version identifier.
        core: Half-open target span between first and last paired bases.
        config: Current accepted candidate curation registry.

    Returns:
        Intended, outside-scope, or nonmember reference classification.
    """
    targets = config["candidates"][candidate]["curation"]["intended_targets"]
    members = [r for r in targets if r["reference_id"] == reference_id]
    if not members:
        return "unintended_reference"
    if any(r["interval_0based_half_open"][0] <= core[0]
           and core[1] <= r["interval_0based_half_open"][1] for r in members):
        return "intended_interval"
    return "unintended_outside_or_crossing_interval"


def mapped_blocks(outcome, source):
    """Replay native masking order for a completed single-pair invocation.

    Args:
        outcome: Verified native outcome with ordered reported blocks.
        source: Exact expected input ledger row.

    Returns:
        Ordered original blocks enriched with mappings and draft decisions.
    """
    masked, result = set(), []
    for block in outcome["blocks"]:
        mapping = map_block(block, source, masked)
        result.append({**block, **mapping, "draft_qualifies": qualifies(block),
                       **alignment_metrics(block["alignment_rows"])})
        masked.update(range(*mapping["core_target_interval"]))
    return result


def report_ceiling(length):
    """Bound report count without using an empirical truncation threshold.

    Args:
        length: Positive target length within the verified native ceiling.

    Returns:
        L+1, exceeding L reports if each consumes a fresh target position.
    """
    if type(length) is not int or not 0 < length <= 109224:
        raise ValueError("Invalid target length for the closure bound")
    return length + 1


def check_progress(blocks, length):
    """Check the fresh-mask premise against every mapped report in order.

    Args:
        blocks: Strictly source-verified mapped native blocks.
        length: Exact whole-target length.

    Returns:
        Counts proving finite progress for this observed greedy trace.
    """
    masked, fresh_counts = set(), []
    for block in blocks:
        start, end = block["core_target_interval"]
        paired = {edge[1] for edge in block["paired_edges"]}
        if (not 0 <= start < end <= length or not paired
                or not paired <= set(range(start, end))
                or paired & masked):
            raise ValueError("Empty, invalid or previously masked paired core")
        fresh = set(range(start, end)) - masked
        if not fresh:
            raise ValueError("Native report made no masking progress")
        fresh_counts.append(len(fresh))
        masked.update(range(start, end))
    if len(blocks) > length:
        raise ValueError("Report count contradicts the fresh-mask bound")
    return {"reports": len(blocks), "masked_positions": len(masked),
            "minimum_fresh_positions_per_report": min(fresh_counts)
            if fresh_counts else None, "progress_verified": True}
