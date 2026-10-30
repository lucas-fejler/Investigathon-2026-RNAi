"""Aggregate reported positions into regional results and rankings."""

from collections import Counter, defaultdict
import csv
from pathlib import Path

import ranking
import puntuacion as scoring
from evaluar_interacciones import (assess_pair, build_identity, catalog,
                              position_state)
from archivos import (digest, fasta, read_json, read_rows, write_json,
                          write_rows)
from seleccion_panel import ROOT, validate_manifest


def validate_population(manifest, regions):
    """Require exactly the frozen population without relaxing exhaustive M8.

    Args:
        manifest: Frozen admission manifest.
        regions: Rows offered for panel scoring.

    Returns:
        None; missing, duplicate, extra or altered identities raise ValueError.
    """
    expected = {r["region_id"]: r for r in manifest["regions"]}
    rows = list(regions)
    if len(rows) != len(expected) or {r["region_id"] for r in rows} != set(
            expected):
        raise ValueError("Panel population differs from frozen admission")
    findings = scoring.Findings()
    for row in rows:
        scoring.validate_region_identity(row, findings)
        frozen = expected[row["region_id"]]
        for field in ("candidate_id", "parent_id", "interval_1based_inclusive",
                      "position_count", "a_numerator", "a_denominator"):
            if row[field] != frozen[field]:
                raise ValueError("Changed frozen regional identity: " + field)
    if findings.items:
        raise ValueError("Invalid panel region identity")


def score_panel(manifest, activity_positions, penalty_positions):
    """Apply M8 aggregation/arithmetic and hold ranks at gene completeness.

    Args:
        manifest: Explicit frozen panel.
        activity_positions: Original M6 positions keyed by window ID.
        penalty_positions: Verified nullable P positions keyed by window ID.

    Returns:
        Region measurements, gene coverage and exact within-panel rankings.
    """
    validate_population(manifest, manifest["regions"])
    findings = scoring.Findings()
    rows = []
    for frozen in manifest["regions"]:
        row = scoring.score_region(
            frozen, {"positions_by_window": activity_positions},
            {"positions_by_window": penalty_positions}, findings)
        row["population_status"] = "selected_complete" if row[
            "status"] == "scored" else "selected_incomplete"
        row["p_withholding_details"] = dict(Counter(
            cause for wid in frozen["window_ids"]
            for cause in penalty_positions.get(wid, {}).get(
                "withholding_reasons", ["missing_required_position"])
            if penalty_positions.get(wid, {}).get("status") != "measured"))
        if row["p_unavailable_positions"]:
            row.update({"flagged": None, "p_exact": None, "scores": None})
        if row["a_unavailable_positions"]:
            row.update({"favorable": None, "a_exact": None, "scores": None})
        row["organism"] = next(g["metadata"]["curation"]["organism"]
                               for g in manifest["genes"]
                               if g["candidate_id"] == row["candidate_id"])
        rows.append(row)
    serious = [r for r in findings.items if r["kind"] != "missing_position"]
    if serious:
        raise ValueError("M6/M8 adapter inconsistency: " + str(serious))
    genes, ranking = [], []
    for gene in manifest["genes"]:
        selected = [r for r in rows
                    if r["candidate_id"] == gene["candidate_id"]]
        complete = sum(r["status"] == "scored" for r in selected)
        eligible = complete == gene["admitted_count"]
        causes = sorted({cause for r in selected
                         for cause in (r["withheld_causes"] or [])})
        coverage = {"candidate_id": gene["candidate_id"],
                    "organism": gene["metadata"]["curation"]["organism"],
                    "admitted": gene["admitted_count"],
                    "complete": complete,
                    "incomplete": len(selected) - complete,
                    "outside_panel": gene["outside_panel_count"],
                    "ranking_released": eligible,
                    "withholding_reasons": causes}
        genes.append(coverage)
        for row in selected:
            row["ranking_withholding_reasons"] = [] if eligible else causes
        if eligible:
            for weight in ("0.5", "1", "2"):
                ranking.extend(ranking.rank_gene(selected, weight))
    return rows, genes, ranking


def write_table(path, rows, columns):
    """Export a compact TSV with blank numeric nulls.

    Args:
        path: Destination.
        rows: Flat export rows.
        columns: Ordered column names.

    Returns:
        None.
    """
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def export(output):
    """Reverify native checkpoints and write real-data panel results/report.

    Args:
        output: Existing frozen run directory.

    Returns:
        Summary of completed and incomplete regional coverage.
    """
    output = Path(output)
    manifest = validate_manifest(output)
    binding = read_json(output / "execution_binding.json")
    if build_identity()["evaluator_id"] != binding["evaluator_id"]:
        raise ValueError("Native evaluator changed since assessment")
    records = catalog(manifest)
    config = read_json(ROOT / "configuracion.json")
    queries = {r["query_id"]: r for r in
               read_rows(output / "selected_queries.jsonl")}
    windows = list(read_rows(output / "selected_windows.jsonl"))
    by_query = defaultdict(list)
    pair_rows, call_counts, call_states = [], Counter(), Counter()
    for directory in sorted((output / "evidence").glob("*")):
        saved = read_json(directory / "binding.json")
        qid = saved["query"]["query_id"]
        rid = saved["reference"]["reference_id"]
        if (saved["binding"] != binding or saved["query"] != queries[qid]
                or saved["reference"] != records[rid]):
            raise ValueError("Cached query/reference/classification changed")
        source_file = next(iter(sorted(directory.glob("*/source.json"))), None)
        if source_file is None:
            result = {"state": "unavailable_interrupted", "flag": None,
                      "evidence_key": directory.name, "query_id": qid,
                      "reference_id": rid, "blocks": []}
        else:
            sequence = read_json(source_file)["target_sequence"]
            result = assess_pair(output, queries[qid], records[rid], sequence,
                                 saved["candidate"], binding, config,
                                 execute=False)
        by_query[qid].append(result)
        for call in directory.glob("*/call.json"):
            call_counts[read_json(call)["status"]] += 1
        for assessment in directory.glob("*/assessment.json"):
            call_states[read_json(assessment)["state"]] += 1
        pair_rows.append({k: result.get(k) for k in (
            "evidence_key", "query_id", "reference_id", "state", "flag",
            "organism", "witnesses", "error")})
    penalties = {}
    for window in windows:
        qids = {o["query_id"] for o in window["occurrences"]}
        results = [r for qid in qids for r in by_query[qid]]
        state, flag = position_state(results, len(records) * len(qids))
        witness_pairs = [r for r in results
                         if r["state"] == "positive_witness"]
        orientations = [o["orientation"] for o in window["occurrences"]
                        if any(r["state"] == "positive_witness"
                               for r in by_query[o["query_id"]])]
        penalties[window["window_id"]] = {
            "window_id": window["window_id"], "parent_id": window["parent_id"],
            "interval_1based_inclusive": window["interval_1based_inclusive"],
            "status": "measured" if flag is not None else "unavailable",
            "coverage_state": state, "flagged": flag,
            "qualifying_site_count": sum(len(r.get("witnesses", []))
                                         for r in witness_pairs),
            "organisms": sorted({r["organism"] for r in witness_pairs}),
            "orientations": orientations,
            "evidence_keys": [r["evidence_key"] for r in results],
            "required_pairs": len(records) * len(qids),
            "attempted_pairs": len(results),
            "unattempted_pairs": len(records) * len(qids) - len(results),
            "unattempted_state": "not_required_after_positive" if flag is True
            else "required_untested" if flag is None else "none",
            "withholding_reasons": sorted({r["state"] for r in results
                                           if r["flag"] is None}) or (
                ["missing_required_coverage"] if flag is None else []),
        }
    a_path = ROOT / "datos_preparados/asimetria/asymmetry_positions.jsonl"
    a_positions = {r["window_id"]: r for r in read_rows(a_path)}
    rows, genes, ranking = score_panel(manifest, a_positions, penalties)
    rank_by = {(r["region_id"], r["lambda"]): r for r in ranking}
    flat = []
    for row in rows:
        entry = {"candidate_id": row["candidate_id"],
                 "region_id": row["region_id"], "organism": row["organism"],
                 "start_1based": row["interval_1based_inclusive"][0],
                 "end_1based": row["interval_1based_inclusive"][1],
                 "A": row["a_exact"], "P": row["p_exact"],
                 "favorable": row["favorable"], "flagged": row["flagged"],
                 "status": row["population_status"],
                 "missing_P_positions": row["p_unavailable_positions"],
                 "withholding_reasons": ";".join(row["withheld_causes"] or []),
                 "ranking_withholding_reasons": ";".join(row[
                     "ranking_withholding_reasons"])}
        for weight in ("0.5", "1", "2"):
            rank = rank_by.get((row["region_id"], weight), {})
            entry["S_" + weight] = row["scores"][weight]["exact"] if row[
                "scores"] else None
            entry["rank_" + weight] = rank.get("rank")
            entry["tie_" + weight] = rank.get("tie_group_id")
            entry["tie_size_" + weight] = rank.get("tie_group_size")
        flat.append(entry)
    write_rows(output / "p_positions.jsonl", penalties.values())
    write_rows(output / "pair_status.jsonl", pair_rows)
    write_rows(output / "region_results.jsonl", rows)
    write_rows(output / "rankings.jsonl", ranking)
    write_json(output / "gene_coverage.json", genes)
    write_table(output / "region_results.tsv", flat, list(flat[0]))
    write_table(output / "gene_coverage.tsv", genes, list(genes[0]))
    selected_states = {r["region_id"]: r["population_status"] for r in rows}
    write_rows(output / "population_status.jsonl", (
        {"region_id": r["region_id"], "candidate_id": r["candidate_id"],
         "interval_1based_inclusive": r["interval_1based_inclusive"],
         "status": selected_states.get(r["region_id"],
                                       "outside_declared_panel")} for r in
        read_rows(ROOT / "datos_preparados/asimetria/asymmetry_regions.jsonl")))
    complete = sum(r["status"] == "scored" for r in rows)
    recovery_states = Counter(read_json(p)["status"] for p in
                              (output / "evidence").glob("*/recovery.json"))
    summary = {"kind": "preliminary_combined_assessment_of_A_selected_panel",
               "panel_id": manifest["panel_id"], "binding": binding,
               "complete_regions": complete,
               "incomplete_regions": len(rows) - complete,
               "genes_ranked": sum(g["ranking_released"] for g in genes),
               "position_states": dict(Counter(p["coverage_state"]
                                               for p in penalties.values())),
               "native_call_status": dict(call_counts),
               "native_assessment_states": dict(call_states),
               "recovery_states": dict(recovery_states),
               "unavailable_attempted_pairs": sum(r["state"] not in (
                   "positive_witness", "negative_reported_pair")
                   for r in pair_rows),
               "pair_states": dict(Counter(r["state"] for r in pair_rows)),
               "unique_pairs_attempted": len(pair_rows),
               "scope_reference_count": len(records),
               "fixture_P_used": False,
               "exhaustive_milestones_complete": False}
    write_json(output / "results_summary.json", summary)
    inventory = []
    for directory in sorted((output / "evidence").glob("*")):
        inventory.append({
            "evidence_key": directory.name,
            "binding_sha256": digest(directory / "binding.json"),
            "result_sha256": digest(directory / "result.json"),
            "recovery_sha256": digest(directory / "recovery.json")
            if (directory / "recovery.json").exists() else None,
            "calls_sha256": {str(p.relative_to(directory)): digest(p)
                             for p in sorted(directory.glob("*/call.json"))},
            "raw_dependency_hashes": "Each call.json binds raw stdout, "
            "stderr, audit, exact FASTA inputs and source ledger.",
        })
    write_rows(output / "native_evidence_inventory.jsonl", inventory)
    lines = ["# Preliminary combined assessment of an A-selected panel", "",
             "The 2026-10-29 deadline amendment freezes up to three regions "
             "per gene from finished M6 measurements before inspecting P. "
             "Admission uses exact A descending, then ascending coordinates "
             "and identifier solely for budget sampling. Full A ties and "
             "excluded tied regions remain in panel_manifest.json. "
             "Administrative sampling gives no biological superiority. "
             "Overlap is reported without filtering.", "",
             "A = favorable modeled-asymmetry positions / 180; "
             "P = positions with a qualifying unintended interaction / 180; "
             "Sλ = A − λP for λ = 0.5, 1, 2. The 200 bp regions each contain "
             "180 overlapping 21 nt positions. Either query orientation may "
             "flag a position once. M8 exact Fraction arithmetic and tied "
             "competition ranks are reused unchanged.", "",
             "RNAhybrid 2.1.2 uses whole stored-forward records, energy ≤ −10 "
             "kcal/mol and (at least 18 WC/GU or 16 consecutive WC pairs), "
             "with loops ≤15, -p 1 -d 0,1 and no forced helix. Strict audit "
             "BEGIN/END and source coordinate/masking checks are required. "
             "Qualifying paired cores inside candidate-specific intended "
             "features are excluded; other hits in the same record remain "
             "eligible. A five-report prefix extends to L+1 when undecided. "
             "Fresh-mask progress and prefix agreement are checked.", "",
             "The priority scanner uses exact 8 nt seeds and ungapped "
             "WC/GU similarity to order a sparse list of whole records. "
             "It cannot certify negatives. Every undecided position falls "
             "back to all remaining query/reference pairs. A verified "
             "unintended witness settles Boolean P; remaining unattempted "
             "pairs are labeled not required after positive, never negative. "
             "Failed, interrupted, ambiguous or unresolved capped work stays "
             "unavailable and blocks its position. Reuse binds query, "
             "reference, evaluator policy and classification identity.", "",
             "An interrupted attempt remains archived as unavailable. "
             "An explicit recovery can complete only that exact task under "
             "unchanged identities; archived bytes and replacement calls "
             "are bound together. Another positive cannot resolve failed "
             "required work. Recovery states: " + str(dict(recovery_states))
             + ". Capped five-hit prefixes in native_assessment_states may "
             "have been resolved by their complete continuation; final pair "
             "states are reported separately.", "",
             "Reference scope: original snapshot manifest SHA-256 `" +
             manifest["reference_scope"]["manifest"]["sha256"] + "`, " +
             str(len(records)) + " whole records:", "",
             "| Dataset | Organism | Records |", "| --- | --- | ---: |"]
    for ds in manifest["reference_scope"]["datasets"]:
        lines.append(f"| {ds['id']} | {ds['organism']} | "
                     f"{ds['statistics']['records']} |")
    lines += ["", "The original dataset inclusion policies, releases, "
              "taxonomy and limitations are retained in the frozen manifest. "
              "Seven supplied candidates are Dalbulus maidis; two are Maize "
              "rayado fino virus (current taxonomic assignment Marafivirus "
              "maydis). These are separate from reference species.",
              "", f"Actual results: **{complete}/{len(rows)} regions complete**."
              f"**{summary['genes_ranked']}/9 genes released for ranking**. "
              f"{len(pair_rows)} unique pairs attempted; native call states "
              f"{dict(call_counts)}. Position states: "
              f"{summary['position_states']}.", "",
              "| Candidate | Organism | Complete / admitted | Outside panel | "
              "Rank released | Withholding |",
              "| --- | --- | ---: | ---: | --- | --- |"]
    for gene in genes:
        lines.append(f"| {gene['candidate_id']} | {gene['organism']} | "
                     f"{gene['complete']}/{gene['admitted']} | "
                     f"{gene['outside_panel']} | {gene['ranking_released']} | "
                     f"{', '.join(gene['withholding_reasons'])} |")
    lines += ["", "| Candidate | Coordinates (1-based inclusive) | A | P | "
              "S0.5 | S1 | S2 | Ranks (0.5 / 1 / 2) |",
              "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in flat:
        ranks = " / ".join(str(row["rank_" + w] or "withheld")
                           for w in ("0.5", "1", "2"))
        scores = " | ".join(str(row["S_" + w]) if row["S_" + w] is not None
                            else "null" for w in ("0.5", "1", "2"))
        lines.append(f"| {row['candidate_id']} | {row['start_1based']}–"
                     f"{row['end_1based']} | {row['A']} | {row['P']} | "
                     f"{scores} | {ranks} |")
    lines += ["", "Exact tie group IDs and sizes for each λ are in "
              "rankings.jsonl and region_results.tsv; coordinates arrange "
              "tied rows for display only. Gene rankings are withheld unless "
              "all frozen admissions have complete A and P. Missing P and "
              "dependent scores remain null. Original IDs are preserved. "
              "The complete original population is labeled in "
              "population_status.jsonl.", "",
              "Historical M4/M5 summaries have been retained as provenance, "
              "but their referenced raw invocation directories are absent. "
              "No completed native result was reused from those summaries. "
              "Synthetic controls and retained-case replays are under "
              "runs/deadline-verification, separate from panel evidence. "
              "No synthetic-P ranking selected or scored a candidate.", "",
              "Focused validation and exact commands are recorded in "
              "../REPRODUCIBILIDAD.md. Timing, initial workload, "
              "environment and resumable progress are in workload.json, "
              "environment.json, priority_identity.json and execution.json. "
              "Every native invocation is retained under evidence/.", "",
              "This is a preliminary combined assessment of an A-selected "
              "panel. Excluded lower-A regions could outperform it under the "
              "combined score. A is a modeled descriptor. Computational "
              "interaction evidence establishes neither biological efficacy "
              "nor safety. Negatives, if available, concern the pinned native "
              "greedy reporting policy, not all duplex structures. The "
              "reference snapshots do not cover all isoforms, populations, "
              "isolates or biological contexts. Exhaustive production "
              "coverage and exhaustive M4/M5/M9 acceptance are deferred.",
              ""]
    (output / "REPORT.md").write_text("\n".join(lines))
    exports = [p for p in output.iterdir() if p.is_file()
               and p.name not in ("export_manifest.json",)
               and p.suffix in (".json", ".jsonl", ".tsv", ".md")]
    write_json(output / "export_manifest.json", {
        "panel_id": manifest["panel_id"],
        "artifacts_sha256": {p.name: digest(p) for p in exports}})
    return summary
