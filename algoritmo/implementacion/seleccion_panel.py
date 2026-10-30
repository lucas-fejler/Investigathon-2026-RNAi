"""Freeze A-based admission before interaction assessment."""

from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import ranking
import puntuacion as scoring
from archivos import (digest, fasta, identity, read_json, read_rows,
                          reverse_complement, sequence_sha, write_json,
                          write_rows)


ROOT = Path(__file__).resolve().parent
POLICY = {
    "qualifies": "MFE <= -10 AND (WC+GU >=18 OR consecutive WC >=16)",
    "native_options": "-c -e -10 -p 1 -d 0,1 -u 15 -v 15 -m L -n 21",
    "cap": "L+1; verified fresh paired-core masking progress",
    "query_orientations": ["+", "-"],
    "reference_orientation": "whole_record_stored_forward",
    "negative_scope": "complete_native_greedy_reports_not_all_structures",
    "fallback": "all_remaining_records_for_undecided_positions",
    "positive_stop": "one qualifying unintended witness settles position",
    "required_failure": "unavailable_even_with_another_positive",
    "intended": "candidate_specific_paired_core_containment",
    "denominator": 180,
}


def admission(rows, limit):
    """Select a deterministic administrative sample and retain whole A ties.

    Args:
        rows: Complete M6 region rows for one gene.
        limit: Maximum admitted regions, a positive integer.

    Returns:
        Selected rows and tie groups intersecting the admitted sample.
    """
    if type(limit) is not int or limit < 1 or not rows:
        raise ValueError("Admission needs a positive budget and regions")
    ordered = sorted(rows, key=lambda r: (
        -Fraction(r["a_numerator"], r["a_denominator"]),
        r["interval_1based_inclusive"], r["region_id"]))
    selected = ordered[:limit]
    admitted = {r["region_id"] for r in selected}
    ties = []
    for value in sorted({Fraction(r["a_exact"]) for r in selected},
                        reverse=True):
        members = [r["region_id"] for r in ordered
                   if Fraction(r["a_exact"]) == value]
        ties.append({"a_exact": str(value), "all_members": members,
                     "admitted": [r for r in members if r in admitted],
                     "excluded": [r for r in members if r not in admitted]})
    return selected, ties


def load_inputs():
    """Bind synchronized parents, M3 mappings and completed M6 measurements.

    Args:
        None.

    Returns:
        Validated original inputs and their artifact identities.
    """
    config = read_json(ROOT / "configuracion.json")
    m3 = ROOT / "datos_preparados/enumeracion"
    m6 = ROOT / "datos_preparados/asimetria"
    parents = {r["parent_id"]: r for r in read_json(m3 / "parents.json")}
    sequences = dict(fasta(m3 / "normalized.fasta"))
    if sequences != dict(fasta(ROOT / "sequences/genes_9.fasta")):
        raise ValueError("Supplied and normalized parents disagree")
    if set(sequences) != set(config["candidates"]) or len(sequences) != 9:
        raise ValueError("Expected the original nine candidates")
    for pid, row in parents.items():
        gene = row["candidate_id"]
        sha = sequence_sha(sequences[gene])
        if (sha != row["normalized_sequence_sha256"]
                or pid != gene + "@" + sha
                or len(sequences[gene]) != row["sequence_length"]
                or sha != config["candidates"][gene][
                    "normalized_sequence_sha256"]):
            raise ValueError("Parent identity changed")
    queries = {r["query_id"]: r for r in read_rows(m3 / "queries.jsonl")}
    windows = {r["window_id"]: dict(r, occurrences=[]) for r in
               read_rows(m3 / "windows.jsonl")}
    for row in read_rows(m3 / "occurrences.jsonl"):
        windows[row["window_id"]]["occurrences"].append(row)
    for wid, row in windows.items():
        pid = row["parent_id"]
        gene = parents[pid]["candidate_id"]
        start, end = row["interval_0based_half_open"]
        seq = sequences[gene][start:end]
        occs = row["occurrences"]
        if (wid != f"window:{pid}:{start}:{end}" or end - start != 21
                or row["interval_1based_inclusive"] != [start + 1, end]
                or len(occs) != 2
                or {o["orientation"] for o in occs} != {"+", "-"}):
            raise ValueError("Invalid original window mapping")
        for occ in occs:
            orient = occ["orientation"]
            seq_o = seq if orient == "+" else reverse_complement(seq)
            query = queries[occ["query_id"]]
            sha = sequence_sha(seq_o)
            if (query["sequence"] != seq_o
                    or query["sequence_sha256"] != sha
                    or occ["query_id"] != "query:sha256:" + sha
                    or occ["occurrence_id"] !=
                    f"occurrence:{pid}:{start}:{end}:{orient}"):
                raise ValueError("Oriented query differs from source")
    positions = {r["window_id"]: r for r in
                 read_rows(m6 / "asymmetry_positions.jsonl")}
    if set(windows) != set(positions):
        raise ValueError("M3/M6 window population disagreement")
    summary = read_json(m6 / "asymmetry.json")
    if summary["completion_status"] != "completed":
        raise ValueError("M6 is not complete")
    regions = list(read_rows(m6 / "asymmetry_regions.jsonl"))
    by_id = {r["region_id"]: r for r in
             read_rows(m3 / "regions.jsonl")}
    findings = scoring.Findings()
    scoring.validate_region_population(
        regions, {k: v["sequence_length"] for k, v in parents.items()},
        findings)
    for row in regions:
        scoring.validate_region_identity(row, findings)
        src = by_id[row["region_id"]]
        if (row["completion_status"] != "complete"
                or row["a_denominator"] != 180
                or row["position_count"] != 180
                or src["interval_1based_inclusive"] !=
                row["interval_1based_inclusive"]):
            raise ValueError("Invalid complete M6 region")
    if findings.items or len(by_id) != len(regions):
        raise ValueError("Original exhaustive population invalid")
    paths = [m3 / name for name in (
        "manifest.json", "parents.json", "normalized.fasta", "regions.jsonl",
        "windows.jsonl", "occurrences.jsonl", "queries.jsonl")]
    paths += [m6 / name for name in (
        "asymmetry.json", "asymmetry_positions.jsonl",
        "asymmetry_regions.jsonl", "orientation_evidence.json")]
    paths += [ROOT / "configuracion.json", ROOT / "sequences/genes_9.fasta"]
    return {"config": config, "parents": parents, "regions": regions,
            "m3_regions": by_id, "windows": windows, "queries": queries,
            "positions": positions, "sequences": sequences,
            "inputs_sha256": {str(p.relative_to(ROOT)): digest(p)
                              for p in paths}}


def freeze(output, limit=3):
    """Create an immutable panel; refuse to overwrite an existing freeze.

    Args:
        output: New run directory.
        limit: Admission budget per gene.

    Returns:
        Sealed frozen manifest.
    """
    output = Path(output)
    if (output / "panel_manifest.json").exists():
        raise ValueError("Panel already frozen; use run/export to resume")
    data = load_inputs()
    config = data["config"]
    scope_path = ROOT / config["reference_manifest"]["path"]
    if digest(scope_path) != config["reference_manifest"]["sha256"]:
        raise ValueError("Reference manifest identity changed")
    scope = read_json(scope_path)
    groups = defaultdict(list)
    for row in data["regions"]:
        groups[row["candidate_id"]].append(row)
    genes, regions, selected_windows = [], [], {}
    for gene, rows in sorted(groups.items()):
        selected, ties = admission(rows, limit)
        genes.append({"candidate_id": gene,
                      "metadata": config["candidates"][gene],
                      "original_region_count": len(rows),
                      "admitted_count": len(selected),
                      "outside_panel_count": len(rows) - len(selected),
                      "a_ties": ties,
                      "overlaps": ranking.overlap_evidence(selected)})
        for ordinal, row in enumerate(selected, 1):
            pid = row["parent_id"]
            start = row["interval_1based_inclusive"][0] - 1
            wids = [f"window:{pid}:{s}:{s + 21}"
                    for s in range(start, start + 180)]
            pos = [data["positions"][wid] for wid in wids]
            favorable = sum(p["classification"] == "favorable" for p in pos)
            if (favorable != row["a_numerator"]
                    or Fraction(row["a_exact"]) != Fraction(favorable, 180)
                    or any(p["status"] != "measured" for p in pos)):
                raise ValueError("Frozen A disagrees with M6 positions")
            regions.append({**row, "m3_region": data["m3_regions"][
                row["region_id"]], "window_ids": wids,
                "admission_order": ordinal,
                "admission_reason": "exact_A_desc_then_coordinates_then_ID",
                "expected_positions": 180})
            selected_windows.update({wid: data["windows"][wid]
                                     for wid in wids})
    qids = {o["query_id"] for w in selected_windows.values()
            for o in w["occurrences"]}
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "selected_windows.jsonl", [selected_windows[k]
               for k in sorted(selected_windows)])
    write_rows(output / "selected_queries.jsonl", [data["queries"][k]
               for k in sorted(qids)])
    manifest = {
        "kind": "real_A_selected_preliminary_panel",
        "frozen_before_panel_P": True,
        "base_commit": "06a58c407f3d8c8acda3ba321b4824744a564efd",
        "preexisting_changes": ["untracked references/", "untracked tools/"],
        "admission": {"maximum_regions_per_gene": limit,
                      "ordering": "exact A descending, coordinates, ID",
                      "tie_sampling": (
                          "administrative; no biological superiority"),
                      "overlap_filter": None},
        "assessment_policy": POLICY,
        "classification_id": identity({k: v["curation"]["intended_targets"]
                                       for k, v in
                                       config["candidates"].items()}),
        "reference_scope": {"manifest": config["reference_manifest"],
                            "datasets": scope["datasets"],
                            "required_records": sum(d["statistics"]["records"]
                                                    for d in
                                                    scope["datasets"])},
        "inputs_sha256": data["inputs_sha256"],
        "native_policy_sha256": digest(ROOT / "rnahybrid.py"),
        "genes": genes, "regions": regions,
        "expected_coverage": {"genes": len(genes), "regions": len(regions),
                              "regional_positions": len(regions) * 180,
                              "unique_windows": len(selected_windows),
                              "unique_oriented_queries": len(qids)},
        "panel_artifacts_sha256": {name: digest(output / name) for name in (
            "selected_windows.jsonl", "selected_queries.jsonl")},
    }
    manifest["panel_id"] = identity(manifest)
    write_json(output / "panel_manifest.json", manifest)
    return manifest


def validate_manifest(output):
    """Reject changed admission, policy, mappings or upstream input identity.

    Args:
        output: Existing frozen run directory.

    Returns:
        Verified manifest.
    """
    output = Path(output)
    manifest = read_json(output / "panel_manifest.json")
    if manifest["panel_id"] != identity({k: v for k, v in manifest.items()
                                         if k != "panel_id"}):
        raise ValueError("Frozen manifest edited")
    for name, sha in manifest["inputs_sha256"].items():
        if digest(ROOT / name) != sha:
            raise ValueError("Changed upstream input: " + name)
    for name, sha in manifest["panel_artifacts_sha256"].items():
        if digest(output / name) != sha:
            raise ValueError("Changed frozen panel mapping: " + name)
    if (manifest["assessment_policy"] != POLICY
            or manifest["native_policy_sha256"] != digest(
                ROOT / "rnahybrid.py")):
        raise ValueError("Frozen evaluator/classification policy changed")
    ref = manifest["reference_scope"]["manifest"]
    if digest(ROOT / ref["path"]) != ref["sha256"]:
        raise ValueError("Changed declared reference scope")
    return manifest
