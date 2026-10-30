"""Consolidate checked regional scores, ties and species coverage."""

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path

from evaluar_interacciones import (build_identity, catalog, check_recovery,
                              native_options, pair_key)
from archivos import (digest, read_json, read_rows, write_json, write_rows)
from cobertura import (COMPLETE, VALID_PAIRS, position_coverage,
                            search_summary)
from seleccion_panel import ROOT, validate_manifest
import rnahybrid as native
from exportar_resultados import write_table


WEIGHTS = {"0.5": Fraction(1, 2), "1": Fraction(1), "2": Fraction(2)}
POSITION_MEANINGS = {
    "a_complete": "All 180 original M6 positions are measured and the exact "
    "favorable count agrees with frozen A and M6.",
    "p_complete": "All 180 union flags are independently decided and agree "
    "with retained P; species searches may remain incomplete.",
    "ranking_eligible": "Every frozen admission for this gene has complete "
    "A/P, verified exact M8 scores and verified within-panel ranks/ties.",
    "selected_complete": "Frozen admission with complete A/P measurements; "
    "score verification and gene ranking eligibility are separate fields.",
    "selected_incomplete": "Frozen admission with unavailable A and/or P; "
    "missing measurements and dependent scores are null.",
    "positive_complete": "A verified unintended witness settles union P; "
    "unattempted species/records are not thereby assessed.",
    "negative_complete": "Both orientations completed all required records "
    "without qualifying unintended evidence, under the complete native "
    "greedy-report policy; not all possible RNA structures.",
    "unresolved": "Missing, partial, ambiguous or capped required work; "
    "P and dependent scores remain null.",
    "failed": "An unresolved execution, recovery, identity or verification "
    "failure blocks the position even with another positive.",
    "species_partial_unresolved": "Some native pairs assessed in this "
    "species, without a witness or complete-policy negative.",
    "species_untested": "No native pair assessed in this species/position.",
    "not_required_for_union_P_after_positive": "Unattempted work unnecessary "
    "for settled union P; species search coverage remains incomplete.",
    "outside_panel": "Outside frozen admission; no combined panel rank.",
    "negative_reported_pair": "One complete native greedy-report trace has "
    "no qualifying unintended witness; other records/orientations may remain.",
    "historical_attempt_failure": "Retained failed attempt, never a zero. "
    "An exact-task verified retry can resolve the required task; unrelated "
    "positives cannot do so.",
}


def verify_cached_pairs(output, manifest):
    """Bind retained parsed results to raw evidence and unchanged identities.

    Args:
        output: Completed or bounded partial run directory.
        manifest: Verified frozen admission.

    Returns:
        Verified task rows, native attempt counts and verification errors.
        A bad task remains failed; no evaluator invocation is possible here.
    """
    binding = read_json(output / "execution_binding.json")
    if (binding["panel_id"] != manifest["panel_id"]
            or build_identity()["evaluator_id"] != binding["evaluator_id"]):
        raise ValueError("Changed panel/evaluator identity")
    records = catalog(manifest)
    queries = {q["query_id"]: q for q in
               read_rows(output / "selected_queries.jsonl")}
    inventory = {r["evidence_key"]: r for r in
                 read_rows(output / "native_evidence_inventory.jsonl")}
    config = read_json(ROOT / "configuracion.json")
    rows, errors, attempts, seen = [], [], Counter(), set()
    for number, stub in enumerate(read_rows(output / "pair_status.jsonl"), 1):
        key = stub["evidence_key"]
        directory = output / "evidence" / key
        qid, rid = stub["query_id"], stub["reference_id"]
        pair = {**stub, "dataset_id": rid.split("|", 1)[0],
                "greedy_complete": False}
        if (qid, rid) in seen:
            raise ValueError("Duplicate final query/reference task")
        seen.add((qid, rid))
        try:
            item = inventory[key]
            for name, field in (("binding.json", "binding_sha256"),
                                ("result.json", "result_sha256")):
                if digest(directory / name) != item[field]:
                    raise ValueError("Changed retained " + name)
            saved = read_json(directory / "binding.json")
            expected = {"query": queries[qid], "reference": records[rid],
                        "candidate": saved["candidate"], "binding": binding}
            if (saved != expected
                    or saved["candidate"] not in queries[qid]["candidate_ids"]
                    or key != pair_key(queries[qid], records[rid], binding)):
                raise ValueError("Query/reference/classification mismatch")
            recovery = directory / "recovery.json"
            if (digest(recovery) if recovery.exists() else None) != item.get(
                    "recovery_sha256"):
                raise ValueError("Changed recovery evidence")
            if check_recovery(directory, expected):
                raise ValueError("Unresolved required recovery")
            for name, sha in item["calls_sha256"].items():
                path = directory / name
                if digest(path) != sha:
                    raise ValueError("Changed native call metadata")
                call = read_json(path)
                attempts[call["status"]] += 1
                if call["status"] != "completed":
                    if path.parent.name.isdigit():
                        raise ValueError("Unresolved native execution")
                    continue
                if (call["exit_code"] != 0 or call["options"] !=
                        native_options({"target_length": records[rid][
                            "length"]}, call["cap"])):
                    raise ValueError("Native exit/options mismatch")
                for artifact, checksum in call["artifacts_sha256"].items():
                    if digest(path.parent / artifact) != checksum:
                        raise ValueError("Changed raw native dependency")
            result = read_json(directory / "result.json")
            if result["evidence_binding"] != expected:
                raise ValueError("Cached result policy/classifier mismatch")
            native.check_progress(result["blocks"], records[rid]["length"])
            witnesses = set()
            for block in result["blocks"]:
                qualifies = native.qualifies(block)
                classification = native.classify(
                    saved["candidate"], rid,
                    block["core_target_interval"], config)
                if (qualifies != block["draft_qualifies"]
                        or classification != block["classification"]):
                    raise ValueError("Qualification/classification changed")
                if qualifies and classification.startswith("unintended"):
                    witnesses.add(block["duplex_id"])
            if (set(result.get("witnesses", [])) != witnesses
                    or any(result.get(k) != stub.get(k) for k in (
                        "state", "flag", "query_id", "reference_id"))):
                raise ValueError("Final task differs from retained evidence")
            if (result["state"] == "positive_witness") != bool(witnesses):
                raise ValueError("Unintended witness/state disagreement")
            pair["greedy_complete"] = result.get("greedy_complete", False)
            if (result["state"] == "negative_reported_pair"
                    and pair["greedy_complete"] is not True):
                raise ValueError("Negative task has incomplete greedy search")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            pair.update(state="unavailable_verification", flag=None)
            errors.append({"evidence_key": key, "error": str(exc)})
        rows.append(pair)
        if number % 500 == 0:
            print(f"Verified {number} retained tasks; native launches: 0",
                  flush=True)
    return rows, dict(attempts), errors


def independently_score(manifest, activity, positions, producer_regions,
                        producer_rankings):
    """Reaggregate 180 flags and verify scores/ranks without M8 arithmetic.

    Args:
        manifest: Frozen panel with original region identities.
        activity: Original M6 per-window measurements.
        positions: Independently checked union-P positions keyed by ID.
        producer_regions: Existing M8 region records keyed by original ID.
        producer_rankings: Existing exact M8 within-gene ranking rows.

    Returns:
        Separate A/P coverage, gene eligibility, verified rankings and ties.
        Incomplete genes have no released ranking rows.
    """
    rows, genes, released, ties = [], [], [], []
    for frozen in manifest["regions"]:
        rid, pid = frozen["region_id"], frozen["parent_id"]
        start, end = frozen["interval_1based_inclusive"]
        wids = [f"window:{pid}:{s}:{s + 21}"
                for s in range(start - 1, end - 20)]
        if len(wids) != 180 or wids != frozen["window_ids"]:
            raise ValueError("Invalid frozen 180-position identity")
        a = [activity.get(wid, {}) for wid in wids]
        p = [positions.get(wid, {"status": "unresolved"}) for wid in wids]
        counts = Counter(row["status"] for row in p)
        a_count = sum(r.get("status") == "measured" and r.get(
            "classification") in ("favorable", "neutral", "unfavorable")
                      for r in a)
        favorable = sum(r.get("classification") == "favorable" for r in a)
        a_complete = a_count == 180
        reasons = []
        if a_complete and (favorable != frozen["a_numerator"]
                           or Fraction(frozen["a_exact"]) !=
                           Fraction(favorable, 180)):
            a_complete = False
            reasons.append("a_disagrees_with_M6")
        if not a_complete:
            reasons.append("a_incomplete")
        p_complete = sum(counts[s] for s in COMPLETE) == 180 and all(
            r.get("flagged") is (r["status"] == "positive_complete")
            for r in p)
        if not p_complete:
            reasons.append("p_incomplete")
        flagged = (sum(r.get("flagged") is True for r in p)
                   if p_complete else None)
        expected = ({key: str(Fraction(favorable, 180)
                              - weight * Fraction(flagged, 180))
                     for key, weight in WEIGHTS.items()}
                    if a_complete and p_complete else None)
        scores_verified = expected is not None
        producer = producer_regions.get(rid, {})
        if scores_verified:
            try:
                scores_verified = (
                    producer["candidate_id"] == frozen["candidate_id"]
                    and producer["parent_id"] == pid
                    and producer["interval_1based_inclusive"] == [start, end]
                    and producer["status"] == "scored"
                    and producer["favorable"] == favorable
                    and producer["flagged"] == flagged
                    and Fraction(producer["a_exact"]) == Fraction(
                        favorable, 180)
                    and Fraction(producer["p_exact"]) == Fraction(flagged, 180)
                    and all(Fraction(producer["scores"][key]["exact"])
                            == Fraction(value) == Fraction(
                                producer["scores"][key]["numerator"],
                                producer["scores"][key]["denominator"])
                            for key, value in expected.items()))
            except (KeyError, ValueError, TypeError, ZeroDivisionError):
                scores_verified = False
            if not scores_verified:
                reasons.append("m8_score_disagreement")
        row = {"region_id": rid, "candidate_id": frozen["candidate_id"],
               "parent_id": pid, "start_1based": start, "end_1based": end,
               "a_complete": a_complete, "a_measured_positions": a_count,
               "a_numerator": favorable if a_complete else None,
               "a_denominator": 180,
               "a_exact": str(Fraction(favorable, 180))
               if a_complete else None,
               "p_complete": p_complete, "p_numerator": flagged,
               "p_denominator": 180,
               "p_exact": str(Fraction(flagged, 180)) if p_complete else None,
               "positive_complete_positions": counts["positive_complete"],
               "negative_complete_positions": counts["negative_complete"],
               "unresolved_positions": counts["unresolved"],
               "failed_positions": counts["failed"],
               "scores_verified": scores_verified,
               "scores": expected if scores_verified else None,
               "ranking_eligible": False, "withholding_reasons": reasons}
        row["population_status"] = ("selected_complete"
                                    if a_complete and p_complete
                                    else "selected_incomplete")
        rows.append(row)
    for frozen_gene in manifest["genes"]:
        gene = frozen_gene["candidate_id"]
        selected = [r for r in rows if r["candidate_id"] == gene]
        complete_p = [r["p_numerator"] for r in selected if r["p_complete"]]
        eligible = (len(selected) == frozen_gene["admitted_count"] and all(
            r["a_complete"] and r["p_complete"] and r["scores_verified"]
            for r in selected))
        gene_ranks, gene_ties = [], []
        if eligible:
            for key in WEIGHTS:
                values = {r["region_id"]: Fraction(r["scores"][key])
                          for r in selected}
                distinct = sorted(set(values.values()), reverse=True)
                actual = [r for r in producer_rankings
                          if r["candidate_id"] == gene and r["lambda"] == key]
                valid = (len(actual) == len(selected)
                         and {r["region_id"] for r in actual} == set(values))
                if not valid:
                    eligible = False
                    break
                for r in actual:
                    value = values[r["region_id"]]
                    member = next(row for row in selected
                                  if row["region_id"] == r["region_id"])
                    group = sorted(rid for rid, score in values.items()
                                   if score == value)
                    valid &= (r["rank"] == 1 + sum(v > value for v in
                                                    values.values())
                              and Fraction(r["score_exact"]) == value
                              and r["interval_1based_inclusive"] == [
                                  member["start_1based"], member["end_1based"]]
                              and r["favorable"] == member["a_numerator"]
                              and r["flagged"] == member["p_numerator"]
                              and r["tie_group_size"] == len(group)
                              and r["tie_group_id"]
                              == distinct.index(value) + 1)
                if not valid:
                    eligible = False
                    break
                gene_ranks.extend(actual)
                for value in distinct:
                    gene_ties.append({"candidate_id": gene, "lambda": key,
                                      "score_exact": str(value),
                                      "region_ids": sorted(rid for rid, v in
                                                           values.items()
                                                           if v == value)})
        causes = sorted({cause for r in selected
                         for cause in r["withholding_reasons"]})
        if not eligible and not causes:
            causes.append("ranking_disagreement")
        for r in selected:
            r["ranking_eligible"] = eligible
            r["ranking_withholding_reasons"] = [] if eligible else causes
        constant = len(set(complete_p)) == 1 if len(complete_p) == len(
            selected) else None
        genes.append({
            "candidate_id": gene,
            "organism": frozen_gene["metadata"]["curation"]["organism"],
            "admitted_regions": len(selected),
            "a_complete_regions": sum(r["a_complete"] for r in selected),
            "p_complete_regions": len(complete_p),
            "scores_verified_regions": sum(r["scores_verified"]
                                           for r in selected),
            "ranking_eligible": eligible,
            "withholding_reasons": [] if eligible else causes,
            "outside_panel_regions": frozen_gene["outside_panel_count"],
            "completed_p_numerator_distribution": dict(Counter(complete_p)),
            "p_constant_across_complete_panel": constant,
            "constant_p_exact": str(Fraction(complete_p[0], 180))
            if constant else None,
        })
        if eligible:
            released.extend(gene_ranks)
            ties.extend(gene_ties)
    return rows, genes, released, ties


def save_tables(output, stem, rows, jsonl=False):
    """Write explicit JSON and flat TSV exports, preserving numeric nulls.

    Args:
        output: Release directory.
        stem: File basename.
        rows: Structured table rows.
        jsonl: Use JSONL rather than a JSON list.

    Returns:
        Paths written.
    """
    path = output / (stem + (".jsonl" if jsonl else ".json"))
    if jsonl:
        write_rows(path, rows)
    else:
        write_json(path, rows)
    flat = [{k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list))
             else v for k, v in row.items()} for row in rows]
    tsv = output / (stem + ".tsv")
    if flat:
        write_table(tsv, flat, list(flat[0]))
    return [path, tsv] if flat else [path]


def release_report(output, summary, genes, species, region_rows):
    """Describe actual completion, P saturation and conditional search scope.

    Args:
        output: Existing frozen run directory.
        summary: Release verification and accounting.
        genes: Independently verified gene coverage.
        species: Panel-level species search rows.
        region_rows: Original-coordinate measurements and verified scores.

    Returns:
        Path of the preliminary release report.
    """
    position_counts = summary["position_status_counts"]
    counts = ", ".join(f"{key}: {value}" for key, value in
                       position_counts.items())
    stage = summary["run"]["stage"]
    stop_description = (
        "The existing bounded run finished before its configured limit. "
        if stage == "assessment_finished" else
        "The bounded segment stopped at its configured budget boundary. "
        if stage == "segment_budget_exhausted" else
        "This is an incomplete checkpoint of the existing bounded run. ")
    lines = [
        "# Preliminary combined assessment of the frozen A-selected panel",
        "",
        f"Frozen panel `{summary['panel_id']}` is unchanged. "
        f"**{summary['p_complete_regions']}/{len(region_rows)} regions have "
        f"complete P; {summary['ranking_eligible_genes']}/{len(genes)} genes "
        "are ranking-eligible.** "
        "This release concerns only the declared panel, not exhaustive "
        "production coverage.",
        "",
        stop_description +
        "No duplicate assessment or expanded panel was launched for release. "
        "The resumed segment used "
        f"{summary['run']['segment_wall_seconds']:.3f} / 3600 seconds and "
        f"{summary['run']['segment_pairs']:,} / 20,000 new pairs. "
        "Earlier completed evidence and its interruption/recovery history "
        "remain retained. Segment timing excludes preparation and export; "
        "UTC elapsed time includes the interruption.",
        "",
        f"Unique position states: **{counts}**. There are "
        f"{len(region_rows) * 180:,} regional position memberships and "
        f"{summary['unique_positions']:,} unique original windows. "
        "Overlapping regions reuse evidence without changing their fixed "
        "denominator of 180.",
        "",
        f"**{summary['evaluated_unique_pairs']:,} unique pairs evaluated**: "
        f"{summary['positive_pairs']:,} with an unintended witness and "
        f"{summary['negative_reported_pairs']:,} complete-policy pair "
        "negatives. A negative pair alone is not a negative position. "
        f"Of the {summary['initial_unique_pairs']:,} declared pairs, "
        f"**{summary['not_required_after_positive_pairs']:,} were not "
        "attempted because a positive already settled union P**. "
        f"Still-required untested pairs: {summary['still_required_pairs']:,}. "
        "Unattempted work is never counted as evaluated or negative.",
        "",
        "Native invocations, including continuations and recovery, are "
        f"counted separately: {summary['native_attempt_status']}. "
        "The historical interrupted attempt remains failed/unavailable as an "
        "attempt; its exact task was successfully repeated under unchanged "
        "identities. No unresolved required failure is hidden by another hit.",
        "",
        "## A, P and within-panel ranking",
        "",
        "A = favorable modeled-asymmetry positions / 180. "
        "P = positions with at least one qualifying unintended interaction "
        "/ 180. Sλ = A − λP for λ = 0.5, 1, 2. Each position contributes "
        "at most once across sites, orientations and organisms.",
        "",
        "Completed P numerator distribution: "
        f"**{summary['completed_p_numerator_distribution']}**, denominator "
        "180. Constant P is declared only when every admission in a gene "
        "has complete P. A partial set of equal measurements cannot establish "
        "constant P for a withheld panel.",
        "",
        "| Gene | Organism | A complete | P complete | P numerator counts | "
        "Constant P | Rank eligible |",
        "| --- | --- | ---: | ---: | --- | --- | --- |",
    ]
    for gene in genes:
        lines.append(
            f"| {gene['candidate_id']} | {gene['organism']} | "
            f"{gene['a_complete_regions']}/{gene['admitted_regions']} | "
            f"{gene['p_complete_regions']}/{gene['admitted_regions']} | "
            f"{gene['completed_p_numerator_distribution']} | "
            f"{gene['constant_p_exact']} | {gene['ranking_eligible']} |")
    if all(g["constant_p_exact"] == "1" for g in genes):
        lines += [
            "", "**Every one of the nine gene panels has constant P = 1 "
            "(180/180 in all three admitted regions).** Thus S0.5 = A − 0.5, "
            "S1 = A − 1 and S2 = A − 2. Within this selected panel all three "
            "rankings preserve the A ordering and its exact ties. The "
            "penalty does not distinguish these admissions; this is not a "
            "species-specific completeness or biological safety result."]
    lines += [
        "", "## Search coverage by reference species", "",
        "The five reference catalogs were verified during preparation. "
        "Native interaction search coverage is the sparse evaluated subset "
        "below. A union-positive position does not mean every species was "
        "assessed. A positive in a species is a witness, not exhaustive "
        "search of that species. Partial no-witness searches stay unresolved.",
        "",
        "| Reference organism | Records | Native pairs evaluated | Positions "
        "positive in species | Complete-policy negative positions | Partial "
        "no-witness positions | Untested positions |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in species:
        lines.append(
            f"| {row['organism']} | {row['reference_records']:,} | "
            f"{row['evaluated_pairs']:,} | "
            f"{row['positive_complete_positions']} | "
            f"{row['negative_complete_positions']} | "
            f"{row['partial_unresolved_positions']} | "
            f"{row['untested_positions']} |")
    lines += [
        "", "Position-, region- and gene-level species search tables "
        "retain these distinctions, raw pair counts, unattempted work, "
        "exhausted greedy traces and `all_required_pairs_evaluated`. Work "
        "unnecessary for union P remains unassessed for species coverage.",
        "", "## Verification and status meanings", "",
        "The release verifier checked recorded raw dependency hashes, "
        "query/reference/evaluator/classification identities, qualification, "
        "intended-feature classification and failed-attempt recovery links. "
        "It launched zero native evaluations. Independently generated "
        "180-window memberships reaggregate P from Boolean flags and A from "
        "the original M6 classifications. Exact Fraction arithmetic and "
        "competition ranks/tie groups independently verify the existing M8 "
        "outputs. Coordinates and original identifiers are preserved.",
        "",
        "A completion, P completion, score verification and gene ranking "
        "eligibility are separate fields. A gene releases ranks only when "
        "all frozen admissions pass every required check. Missing P and its "
        "dependent scores remain null. Selected incomplete regions and "
        "regions outside the declared panel are distinct statuses. "
        "`release_status_meanings.json` defines every status.",
        "", "## Policy, amendment and limitations", "",
        "Admission remains up to three regions/gene by exact A descending, "
        "then ascending coordinates/identifier solely for budget sampling. "
        "The complete A ties, 12 tied exclusions, overlaps and organism "
        "metadata remain in the original frozen manifest. No overlap filter "
        "was added; coordinate ordering gives no biological superiority.",
        "",
        "RNAhybrid 2.1.2 uses whole stored-forward records, both query "
        "orientations, MFE ≤ −10 kcal/mol and (at least 18 WC/GU pairs OR "
        "16 consecutive WC pairs), -p 1 -d 0,1, loop limits 15 and no forced "
        "helix. Candidate-specific paired-core containment excludes intended "
        "features. Intended first hits do not exempt later unintended hits "
        "in the same record. An undecided five-hit prefix continues to L+1 "
        "with strict parsing, audited completion and fresh-mask checks.",
        "",
        "**Negative means no qualifying unintended interaction in a complete "
        "search under this declared native greedy-report policy. It does "
        "not mean exhaustive enumeration of all possible RNA structures.** "
        "There are no negative-complete union positions in this actual panel. "
        "Empty retrieval, unfinished work and exhausted budgets cannot "
        "establish a negative.",
        "",
        "This is a preliminary combined assessment of an A-selected panel. "
        "The 10,966 excluded regions could outperform it under combined "
        "scoring. A is a modeled descriptor; computational interactions "
        "establish neither biological efficacy nor safety. Reference scope "
        "does not cover every isoform, population, isolate or biological "
        "context. Original exhaustive M4/M5/M9 milestones remain incomplete.",
        "", "## Release files", "",
        "- `release_gene_coverage.tsv` and `release_region_coverage.tsv`: "
        "separate completion, coordinates, exact measurements, scores and "
        "withholding reasons.",
        "- `release_position_status.jsonl`: positive, negative-complete, "
        "unresolved and failed positions.",
        "- `species_position_coverage.jsonl`, `species_region_coverage.tsv`, "
        "`species_gene_coverage.tsv`, `species_search_summary.tsv`: actual "
        "native search coverage without extrapolating union positives.",
        "- `release_rankings.tsv`, `release_exact_ties.json`, "
        "`p_distribution.json`: verified within-panel results.",
        "- `release_verification.json`, `release_manifest.json`: independent "
        "checks, budget accounting, input/output IDs and status meanings.",
        "",
        "The earlier `REPORT.md` and producer exports remain unchanged as "
        "assessment provenance. This release checkpoint supplies the explicit "
        "coverage and independently verified ranking eligibility.", "",
    ]
    path = output / "RELEASE_CHECKPOINT.md"
    path.write_text("\n".join(lines))
    return path


def release(output):
    """Produce a release checkpoint entirely from the existing bounded run.

    Args:
        output: Frozen run directory with retained producer exports.

    Returns:
        Independently checked completion and conditional search accounting.
    """
    output = Path(output).resolve()
    frozen_sha = digest(output / "panel_manifest.json")
    manifest = validate_manifest(output)
    producer_manifest = read_json(output / "export_manifest.json")
    for name, sha in producer_manifest["artifacts_sha256"].items():
        if digest(output / name) != sha:
            raise ValueError("Changed producer export: " + name)
    pairs, attempts, errors = verify_cached_pairs(output, manifest)
    datasets = manifest["reference_scope"]["datasets"]
    windows = list(read_rows(output / "selected_windows.jsonl"))
    producer_p = {r["window_id"]: r for r in
                  read_rows(output / "p_positions.jsonl")}
    pos, species_pos, membership = position_coverage(
        windows, pairs, datasets, producer_p)
    positions = {p["window_id"]: p for p in pos}
    species_index = {(p["window_id"], p["dataset_id"]): p
                     for p in species_pos}
    activity = {p["window_id"]: p for p in read_rows(
        ROOT / "datos_preparados/asimetria/asymmetry_positions.jsonl")}
    m6_regions = {r["region_id"]: r for r in read_rows(
        ROOT / "datos_preparados/asimetria/asymmetry_regions.jsonl")}
    for region in manifest["regions"]:
        if any(region[field] != m6_regions[region["region_id"]][field]
               for field in ("a_numerator", "a_denominator", "a_exact",
                             "interval_1based_inclusive", "parent_id")):
            raise ValueError("Frozen region differs from original M6")
    producer_regions = {r["region_id"]: r for r in
                        read_rows(output / "region_results.jsonl")}
    producer_rankings = list(read_rows(output / "rankings.jsonl"))
    regions, genes, rankings, ties = independently_score(
        manifest, activity, positions, producer_regions, producer_rankings)
    rank_index = {(r["region_id"], r["lambda"]): r for r in rankings}
    for region in regions:
        for weight in WEIGHTS:
            rank = rank_index.get((region["region_id"], weight), {})
            region["S_" + weight] = (region["scores"][weight]
                                     if region["scores"] else None)
            region["rank_" + weight] = rank.get("rank")
            region["tie_size_" + weight] = rank.get("tie_group_size")
    region_species, gene_species = [], []
    for region in manifest["regions"]:
        region_species.extend(search_summary(
            set(region["window_ids"]), positions, species_index, membership,
            pairs, datasets, "region", region["region_id"]))
    for gene in genes:
        wids = {w for r in manifest["regions"]
                if r["candidate_id"] == gene["candidate_id"]
                for w in r["window_ids"]}
        gene["unique_positions"] = len(wids)
        gene["regional_position_memberships"] = gene["admitted_regions"] * 180
        gene["position_status_counts"] = dict(Counter(
            positions[wid]["status"] for wid in wids))
        gene["a_complete"] = gene["a_complete_regions"] == gene[
            "admitted_regions"]
        gene["p_complete"] = gene["p_complete_regions"] == gene[
            "admitted_regions"]
        gene_species.extend(search_summary(
            wids, positions, species_index, membership, pairs, datasets,
            "gene", gene["candidate_id"]))
    species = search_summary(set(positions), positions, species_index,
                             membership, pairs, datasets, "panel",
                             manifest["panel_id"])
    initial_pairs = read_json(output / "workload.json")[
        "unique_query_reference_pairs"]
    if sum(s["required_pairs"] for s in species) != initial_pairs:
        raise ValueError("Search accounting differs from frozen scope")
    counts = Counter(p["status"] for p in pos)
    execution = read_json(output / "execution.json")
    if (execution["segment_wall_seconds"] > 3600
            or execution["segment_pairs"] > 20000):
        raise ValueError("Saved segment exceeds its configured bounds")
    summary = {
        "kind": "preliminary_A_selected_consolidar_resultados",
        "panel_id": manifest["panel_id"],
        "admitted_regions": len(regions), "unique_positions": len(pos),
        "position_status_counts": {s: counts[s] for s in (
            "positive_complete", "negative_complete", "unresolved", "failed")},
        "a_complete_regions": sum(r["a_complete"] for r in regions),
        "p_complete_regions": sum(r["p_complete"] for r in regions),
        "ranking_eligible_genes": sum(g["ranking_eligible"] for g in genes),
        "completed_p_numerator_distribution": dict(Counter(
            r["p_numerator"] for r in regions if r["p_complete"])),
        "initial_unique_pairs": read_json(output / "workload.json")[
            "unique_query_reference_pairs"],
        "attempted_unique_pairs": len(pairs),
        "evaluated_unique_pairs": sum(p["state"] in VALID_PAIRS
                                      for p in pairs),
        "positive_pairs": sum(p["state"] == "positive_witness" for p in pairs),
        "negative_reported_pairs": sum(p["state"] == "negative_reported_pair"
                                       for p in pairs),
        "not_required_after_positive_pairs": sum(
            s["not_required_for_union_P_pairs"] for s in species),
        "still_required_pairs": sum(
            s["still_required_untested_pairs"] + s["unavailable_pairs"]
            for s in species),
        "native_attempt_status": attempts, "cache_verification_errors": errors,
        "run": execution, "native_evaluations_launched_for_release": 0,
        "run_limits": {"wall_seconds": 3600, "new_pairs": 20000},
        "budget_boundary_reached": execution["stage"]
        == "segment_budget_exhausted",
        "m6_regions_checked": len(regions),
        "independently_checked_M8_ranking_rows": len(rankings),
        "exact_score_groups": len(ties),
        "nontrivial_exact_tie_groups": sum(len(t["region_ids"]) > 1
                                           for t in ties),
        "frozen_manifest_sha256": frozen_sha,
        "exhaustive_production_complete": False,
    }
    files = []
    for stem, rows, jsonl in (
            ("release_gene_coverage", genes, False),
            ("release_region_coverage", regions, True),
            ("release_position_status", pos, True),
            ("species_region_coverage", region_species, True),
            ("species_gene_coverage", gene_species, True),
            ("species_search_summary", species, False),
            ("release_rankings", rankings, True)):
        files.extend(save_tables(output, stem, rows, jsonl))
    write_rows(output / "species_position_coverage.jsonl", species_pos)
    write_json(output / "release_exact_ties.json", ties)
    write_json(output / "p_distribution.json", {
        "denominator": 180,
        "completed_region_numerators": summary[
            "completed_p_numerator_distribution"],
        "per_gene": [{k: g[k] for k in (
            "candidate_id", "completed_p_numerator_distribution",
            "p_constant_across_complete_panel", "constant_p_exact")}
                     for g in genes]})
    write_json(output / "release_status_meanings.json", POSITION_MEANINGS)
    write_json(output / "release_verification.json", summary)
    files += [output / name for name in (
        "species_position_coverage.jsonl", "release_exact_ties.json",
        "p_distribution.json", "release_status_meanings.json",
        "release_verification.json")]
    files.append(release_report(output, summary, genes, species, regions))
    if digest(output / "panel_manifest.json") != frozen_sha:
        raise ValueError("Frozen manifest changed during release")
    write_json(output / "release_manifest.json", {
        "panel_id": manifest["panel_id"],
        "frozen_manifest_unchanged": True,
        "frozen_manifest_sha256": frozen_sha,
        "producer_export_manifest_sha256": digest(
            output / "export_manifest.json"),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "release_code_sha256": {name: digest(ROOT / name) for name in (
            "consolidar_resultados.py", "cobertura.py")},
        "artifacts_sha256": {p.name: digest(p) for p in files}})
    return summary


def main():
    """Verify and release only existing evidence within the frozen scope.

    Args:
        None. Command-line arguments select the existing output directory.

    Returns:
        None.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/deadline-panel"))
    args = parser.parse_args()
    print(json.dumps(release(args.output_dir), indent=2))


if __name__ == "__main__":
    main()
