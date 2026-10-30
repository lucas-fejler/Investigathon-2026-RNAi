"""Recalculate A from FASTA and panel scores from reported interaction flags.

Python 3.11+, standard library only. This command never launches RNAhybrid.
"""

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "implementacion"))

import asimetria as asymmetry
import ranking
import puntuacion
from evaluar_interacciones import position_state
from archivos import fasta, read_json, read_rows, reverse_complement
from seleccion_panel import admission


def require(condition, message):
    """Stop at the first incompatible input or numerical disagreement.

    Args:
        condition: Required invariant.
        message: Explanation of failure.

    Returns:
        None.
    """
    if not condition:
        raise ValueError(message)


def digest(path):
    """Hash exact artifact bytes without newline conversion.

    Args:
        path: Existing input file.

    Returns:
        SHA-256 string.
    """
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def measure_population(sequences, table, table_sha):
    """Enumerate every 200 nt region and calculate its modeled asymmetry.

    Args:
        sequences: Canonical candidate sequences keyed by identifier.
        table: Original dsRIP terminal energy table.
        table_sha: Exact table identity.

    Returns:
        All regional A records and all per-position classifications.
    """
    regions, activity = [], {}
    for gene, sequence in sorted(sequences.items()):
        require(len(sequence) >= 200 and set(sequence) <= set("ACGT"),
                f"Invalid candidate sequence: {gene}")
        parent = gene + "@" + hashlib.sha256(sequence.encode("ascii")).hexdigest()
        classifications = []
        for start in range(len(sequence) - 20):
            row = asymmetry.measure_window(sequence[start:start + 21], table,
                                           table_sha)
            require(row["status"] == "measured", f"Unavailable A: {gene}:{start}")
            classifications.append(row["classification"])
            activity[f"window:{parent}:{start}:{start + 21}"] = row
        for start in range(len(sequence) - 199):
            count = Counter(classifications[start:start + 180])
            favorable = count["favorable"]
            regions.append({
                "candidate_id": gene, "parent_id": parent,
                "region_id": f"region:{parent}:{start}:{start + 200}",
                "interval_1based_inclusive": [start + 1, start + 200],
                "a_numerator": favorable, "a_denominator": 180,
                "a_exact": str(Fraction(favorable, 180)),
                "favorable": favorable, "neutral": count["neutral"],
                "unfavorable": count["unfavorable"],
            })
    return regions, activity


def verify_positions(panel, sequences, manifest):
    """Join frozen positions to FASTA, query identities and saved pair decisions.

    This verifies consistency of reported witnesses, not their raw alignments.

    Args:
        panel: Retained panel directory.
        sequences: Canonical sequences.
        manifest: Frozen admission and reference scope.

    Returns:
        Position flags keyed by original window identifier and pair count.
    """
    windows = list(read_rows(panel / "selected_windows.jsonl"))
    queries = {r["query_id"]: r for r in read_rows(panel / "selected_queries.jsonl")}
    positions = {r["window_id"]: r for r in
                 read_rows(panel / "release_position_status.jsonl")}
    required_ids = {wid for r in manifest["regions"] for wid in r["window_ids"]}
    require({w["window_id"] for w in windows} == set(positions) == required_ids,
            "Frozen position population changed")
    require(len(windows) == len(positions), "Duplicate frozen position")
    pairs = list(read_rows(panel / "pair_status.jsonl"))
    require(len({(p["query_id"], p["reference_id"]) for p in pairs}) == len(pairs),
            "Duplicate query/reference task")
    by_query = defaultdict(list)
    for pair in pairs:
        require(pair["query_id"] in queries, "Pair outside selected query population")
        if pair["state"] == "positive_witness":
            require(pair["flag"] is True and bool(pair["witnesses"]),
                    "Positive task without a reported witness")
        by_query[pair["query_id"]].append(pair)
    reference_count = manifest["reference_scope"]["required_records"]
    for window in windows:
        wid = window["window_id"]
        parent = window["parent_id"]
        gene = parent.split("@")[0]
        start, end = window["interval_0based_half_open"]
        sequence = sequences[gene][start:end]
        require(end - start == 21 and len(sequence) == 21 and
                window["interval_1based_inclusive"] == [start + 1, end] and
                wid == f"window:{parent}:{start}:{end}", "Position mapping changed")
        occurrences = window["occurrences"]
        require(len(occurrences) == 2 and
                {o["orientation"] for o in occurrences} == {"+", "-"},
                "Both query orientations are required")
        qids = set()
        for occurrence in occurrences:
            oriented = (sequence if occurrence["orientation"] == "+"
                        else reverse_complement(sequence))
            sha = hashlib.sha256(oriented.encode("ascii")).hexdigest()
            qid = occurrence["query_id"]
            require(qid == "query:sha256:" + sha and
                    queries[qid]["sequence"] == oriented, "Query identity changed")
            qids.add(qid)
        tasks = [p for qid in qids for p in by_query[qid]]
        state, flag = position_state(tasks, len(qids) * reference_count)
        saved = positions[wid]
        require(state in {"positive_complete", "negative_complete"} and
                saved["status"] == state and saved["flagged"] is flag,
                f"Unavailable or inconsistent reported P: {wid}")
    return positions, len(pairs)


def write_table(path, rows, fields):
    """Write a tabular recalculation to a new output directory.

    Args:
        path: Destination TSV.
        rows: Dictionaries with calculated values.
        fields: Ordered column names.

    Returns:
        None.
    """
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reproduce(output):
    """Rebuild admission, asymmetry, P aggregation and all panel rankings.

    Args:
        output: New or empty directory for derived results.

    Returns:
        Scientific scope and successful check counts.
    """
    require(not output.exists() or not any(output.iterdir()),
            "Use a new or empty output directory; existing results are preserved")
    integrity = read_json(ROOT / "resultados/integridad.json")
    for name, checksum in integrity["files_sha256"].items():
        require(digest(ROOT / name) == checksum, f"Artifact bytes changed: {name}")
    fasta_path = ROOT.parent / "seleccion_genes/secuencias/genes_9.fasta"
    require(digest(fasta_path) == integrity["input_fasta_sha256"], "FASTA identity changed")
    entries = list(fasta(fasta_path))
    sequences = dict(entries)
    require(len(entries) == len(sequences) == 9, "Expected nine unique FASTA entries")
    panel = ROOT / "resultados/panel"
    manifest = read_json(panel / "panel_manifest.json")
    require(digest(fasta_path) == manifest["inputs_sha256"]["sequences/genes_9.fasta"],
            "FASTA differs from the original assessment")
    table, identity = asymmetry.load_energy_table(
        ROOT / "implementacion/parametros_energia/dsrip_terminal_energies.tsv",
        ROOT / "implementacion/parametros_energia/dsrip_lookup.db")
    all_regions, activity = measure_population(
        sequences, table, identity["energy_table_export_sha256"])
    by_gene = defaultdict(list)
    for row in all_regions:
        by_gene[row["candidate_id"]].append(row)
    require(set(by_gene) == {g["candidate_id"] for g in manifest["genes"]},
            "Candidate population differs from frozen panel")
    selected = []
    for gene in manifest["genes"]:
        population = by_gene[gene["candidate_id"]]
        chosen, ties = admission(population, manifest["admission"]["maximum_regions_per_gene"])
        require(len(population) == gene["original_region_count"] and ties == gene["a_ties"],
                f"A admission or excluded ties changed: {gene['candidate_id']}")
        selected.extend(chosen)
    expected = {r["region_id"]: r for r in manifest["regions"]}
    require({r["region_id"] for r in selected} == set(expected), "Frozen admission changed")
    positions, pair_count = verify_positions(panel, sequences, manifest)
    saved_regions = {r["region_id"]: r for r in
                     read_rows(panel / "release_region_coverage.jsonl")}
    scored = defaultdict(list)
    for row in selected:
        original = expected[row["region_id"]]
        saved = saved_regions[row["region_id"]]
        require(all(row[k] == original[k] for k in
                    ("favorable", "neutral", "unfavorable")), "A disagrees with frozen measurement")
        flags = [positions[wid]["flagged"] for wid in original["window_ids"]]
        require(len(flags) == 180 and len(set(original["window_ids"])) == 180,
                "Each admitted region requires 180 distinct positions")
        flagged = sum(flags)
        require(saved["p_complete"] is True and saved["a_complete"] is True and
                saved["p_numerator"] == flagged and
                saved["a_numerator"] == row["favorable"], "Regional components disagree")
        scores = puntuacion.score_block(row["favorable"], flagged)
        require(all(Fraction(scores[w]["exact"]) == Fraction(saved["scores"][w])
                    for w in scores), "Combined scores disagree")
        scored[row["candidate_id"]].append({**row, "flagged": flagged,
                                           "status": "scored", "scores": scores})
    rankings = [r for gene in sorted(scored) for weight in ("0.5", "1", "2")
                for r in ranking.rank_gene(scored[gene], weight)]
    saved_ranks = list(read_rows(panel / "release_rankings.jsonl"))
    key = lambda row: (row["candidate_id"], row["lambda"], row["region_id"])
    require(sorted(rankings, key=key) == sorted(saved_ranks, key=key),
            "Exact scores, competition ranks or full tie groups disagree")
    output.mkdir(parents=True, exist_ok=True)
    write_table(output / "asimetria_regiones.tsv", [
        {**r, "start_1based": r["interval_1based_inclusive"][0],
         "end_1based": r["interval_1based_inclusive"][1]} for r in all_regions],
        ("candidate_id", "region_id", "start_1based", "end_1based", "favorable",
         "neutral", "unfavorable", "a_exact"))
    write_table(output / "puntajes_panel.tsv", [
        {**r, "start_1based": r["interval_1based_inclusive"][0],
         "end_1based": r["interval_1based_inclusive"][1]} for r in rankings],
        ("candidate_id", "region_id", "start_1based", "end_1based", "lambda",
         "favorable", "flagged", "score_exact", "rank", "tie_group_size", "is_co_best"))
    summary = {
        "candidates": len(sequences), "A_regions_recalculated": len(all_regions),
        "panel_regions": len(selected), "reported_positions_checked": len(positions),
        "reported_pair_decisions_checked": pair_count, "ranking_rows_matched": len(rankings),
        "P_numerator_distribution": dict(Counter(r["flagged"] for rows in scored.values() for r in rows)),
        "artifact_hashes_matched": len(integrity["files_sha256"]),
        "RNAhybrid_evaluations_launched": 0, "raw_native_alignments_revalidated": False,
        "scope": "A recomputed from FASTA; P reconstructed from reported decisions; scores and ties matched exactly.",
    }
    (output / "resumen.json").write_text(json.dumps(summary, indent=2) + "\n",
                                         encoding="utf-8", newline="\n")
    return summary


def main():
    """Run portable numerical reproduction with explicit failure status.

    Args:
        None. Arguments are parsed from the command line.

    Returns:
        Zero on success, one for an incompatible input or calculation.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--salida", type=Path, default=ROOT.parent / "salidas/algoritmo")
    args = parser.parse_args()
    try:
        print(json.dumps(reproduce(args.salida.resolve()), indent=2))
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
