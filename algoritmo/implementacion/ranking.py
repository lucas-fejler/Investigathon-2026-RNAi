"""Rank within genes using exact rational scores and preserve full ties."""

from fractions import Fraction

import puntuacion as scoring

DISPLAY_ORDER_BASIS = "interval start, interval end, region_id"
SELECTION_PENDING_STATUS = "pending_q11"
PENDING_SELECTION_DECISIONS = ("Q11-a", "Q11-b", "Q11-c")


def _score_fraction(record, weight_key):
    """Recover a region's exact score from its stored components.

    Args:
        record: Scored region record.
        weight_key: Lambda key such as "1".

    Returns:
        The score as an exact Fraction.
    """
    entry = record["scores"][weight_key]
    return Fraction(entry["numerator"], entry["denominator"])


def _display_key(record):
    """Build the deterministic display-ordering key for a region.

    Args:
        record: Scored region record.

    Returns:
        Tuple ordering by interval start, interval end, then region id.
    """
    interval = record["interval_1based_inclusive"]
    return (interval[0], interval[1], record["region_id"])


def rank_gene(region_records, weight_key):
    """Rank one gene's scored regions for a single lambda weight.

    Tied regions share the smallest rank in their group, and the next group
    starts at the previous rank plus the group size, so rank numbers reflect
    how many regions score strictly higher.

    Args:
        region_records: Scored region records for one candidate gene.
        weight_key: Lambda key such as "0.5", "1" or "2".

    Returns:
        List of ranking rows ordered by rank then display key.
    """
    scored = [r for r in region_records if r["status"] == "scored"]
    ordered = sorted(scored,
                     key=lambda r: (-_score_fraction(r, weight_key),
                                    _display_key(r)))

    rows = []
    rank = 0
    index = 0
    group_id = 0
    while index < len(ordered):
        value = _score_fraction(ordered[index], weight_key)
        group = [r for r in ordered[index:]
                 if _score_fraction(r, weight_key) == value]
        rank = index + 1
        group_id += 1
        for member in group:
            rows.append({
                "lambda": weight_key,
                "candidate_id": member["candidate_id"],
                "region_id": member["region_id"],
                "interval_1based_inclusive":
                    member["interval_1based_inclusive"],
                "favorable": member["favorable"],
                "flagged": member["flagged"],
                "score_exact": scoring.exact_string(value),
                "score_display": member["scores"][weight_key]["display"],
                "rank": rank,
                "tie_group_id": group_id,
                "tie_group_size": len(group),
                "is_co_best": group_id == 1,
                "display_order_basis": DISPLAY_ORDER_BASIS,
                "biological_claim": False,
            })
        index += len(group)
    return rows


def overlap_evidence(records):
    """Report pairwise interval overlap among regions as neutral evidence.

    Overlap is factual and M9 will need it, but no approved overlap rule exists
    (Q11-b), so nothing is filtered or preferred here.

    Args:
        records: Scored region records, typically a co-best set.

    Returns:
        List of pairwise overlap descriptions, ordered deterministically.
    """
    pairs = []
    ordered = sorted(records, key=_display_key)
    for i, first in enumerate(ordered):
        for second in ordered[i + 1:]:
            a_start, a_end = first["interval_1based_inclusive"]
            b_start, b_end = second["interval_1based_inclusive"]
            overlap = min(a_end, b_end) - max(a_start, b_start) + 1
            pairs.append({
                "region_a": first["region_id"],
                "region_b": second["region_id"],
                "overlap_bp": max(0, overlap),
                "intervals_intersect": overlap > 0,
            })
    return pairs


def co_best(region_records, weight_key):
    """Identify the full set of regions attaining the maximum exact score.

    Args:
        region_records: Scored region records for one candidate gene.
        weight_key: Lambda key.

    Returns:
        Tuple of the co-best records and their shared exact Fraction score, or
        an empty list and None when no region is scored.
    """
    scored = [r for r in region_records if r["status"] == "scored"]
    if not scored:
        return [], None
    best = max(_score_fraction(r, weight_key) for r in scored)
    members = [r for r in scored if _score_fraction(r, weight_key) == best]
    return sorted(members, key=_display_key), best


def selection_for_gene(region_records, weight_key, eligible):
    """Build the pending-selection record for one gene and weight.

    Args:
        region_records: Scored region records for one candidate gene.
        weight_key: Lambda key.
        eligible: Whether the gene's ranking is eligible at all.

    Returns:
        Selection record naming the pending decisions, with no winner chosen.
    """
    if not eligible:
        return {
            "lambda": weight_key,
            "selection_status": "withheld",
            "co_best_size": 0,
            "co_best": [],
            "pending_decisions": list(PENDING_SELECTION_DECISIONS),
            "note": "gene ranking withheld; no selection attempted",
        }
    members, best = co_best(region_records, weight_key)
    return {
        "lambda": weight_key,
        "selection_status": SELECTION_PENDING_STATUS,
        "co_best_score_exact": scoring.exact_string(best) if best is not None
        else None,
        "co_best_size": len(members),
        "co_best": [{
            "region_id": r["region_id"],
            "interval_1based_inclusive":
                r["interval_1based_inclusive"],
            "favorable": r["favorable"],
            "flagged": r["flagged"],
        } for r in members],
        "overlap_evidence": overlap_evidence(members),
        "pending_decisions": list(PENDING_SELECTION_DECISIONS),
        "alternatives_selected": None,
        "overlap_rule_applied": None,
        "winner_chosen": False,
        "note": "Full co-best set reported. The number of alternatives "
                "(Q11-a), the overlap rule (Q11-b) and any diversity "
                "requirement (Q11-c) are unresolved, so no narrowing was "
                "applied.",
    }


def sensitivity_for_gene(region_records, weight_keys, eligible):
    """Summarise how the top set and ranks change across lambda values.

    Args:
        region_records: Scored region records for one candidate gene.
        weight_keys: Ordered lambda keys to compare.
        eligible: Whether the gene's ranking is eligible.

    Returns:
        Sensitivity record for this gene.
    """
    if not eligible:
        return {"status": "withheld", "per_lambda": {}, "baseline_lambda": "1"}

    baseline_members, _ = co_best(region_records, "1")
    baseline_ids = {r["region_id"] for r in baseline_members}
    per_lambda = {}
    for key in weight_keys:
        members, best = co_best(region_records, key)
        member_ids = {r["region_id"] for r in members}
        scores = [_score_fraction(r, key) for r in region_records
                  if r["status"] == "scored"]
        groups = len({s for s in scores})
        per_lambda[key] = {
            "co_best_size": len(members),
            "co_best_score_exact": scoring.exact_string(best)
            if best is not None else None,
            "co_best_region_ids": sorted(member_ids),
            "identical_to_baseline": member_ids == baseline_ids,
            "intersection_with_baseline": len(member_ids & baseline_ids),
            "distinct_score_groups": groups,
            "score_min_exact": scoring.exact_string(min(scores))
            if scores else None,
            "score_max_exact": scoring.exact_string(max(scores))
            if scores else None,
        }

    # Where each baseline co-best region lands under the other weights.
    baseline_tracking = []
    rank_by_lambda = {}
    for key in weight_keys:
        rank_by_lambda[key] = {row["region_id"]: row["rank"]
                               for row in rank_gene(region_records, key)}
    for region_id in sorted(baseline_ids):
        baseline_tracking.append({
            "region_id": region_id,
            "rank_by_lambda": {key: rank_by_lambda[key].get(region_id)
                               for key in weight_keys},
        })

    return {
        "status": "compared",
        "baseline_lambda": "1",
        "per_lambda": per_lambda,
        "baseline_co_best_tracking": baseline_tracking,
        "top_set_stable_across_lambda": all(
            per_lambda[key]["identical_to_baseline"] for key in weight_keys),
    }


def rank_all(region_records, genes, organisms, weight_keys):
    """Rank every eligible gene and build selection and sensitivity records.

    Args:
        region_records: All scored region records.
        genes: Eligibility mapping from `puntuacion.gene_eligibility`.
        organisms: Mapping from candidate id to organism name.
        weight_keys: Ordered lambda keys.

    Returns:
        Tuple of ranking rows, the selection document and the sensitivity
        document.
    """
    by_gene = {}
    for record in region_records:
        by_gene.setdefault(record["candidate_id"], []).append(record)

    ranking_rows = []
    selection = {"by_candidate": {}, "organism_groups": {}}
    sensitivity = {"by_candidate": {}}

    for candidate_id in sorted(by_gene):
        records = by_gene[candidate_id]
        eligible = genes[candidate_id]["ranking_eligible"]
        organism = organisms.get(candidate_id)
        selection["organism_groups"].setdefault(
            organism or "unknown", []).append(candidate_id)

        per_lambda_selection = {}
        for key in weight_keys:
            if eligible:
                ranking_rows.extend(rank_gene(records, key))
            per_lambda_selection[key] = selection_for_gene(
                records, key, eligible)
        selection["by_candidate"][candidate_id] = {
            "candidate_id": candidate_id,
            "organism": organism,
            "ranking_eligible": eligible,
            "regions": len(records),
            "per_lambda": per_lambda_selection,
        }
        sensitivity["by_candidate"][candidate_id] = {
            "candidate_id": candidate_id,
            "organism": organism,
            **sensitivity_for_gene(records, weight_keys, eligible),
        }

    for organism in selection["organism_groups"]:
        selection["organism_groups"][organism].sort()

    ranking_rows.sort(key=lambda row: (
        row["lambda"], row["candidate_id"], row["rank"],
        row["interval_1based_inclusive"][0], row["region_id"]))
    selection["selection_status"] = SELECTION_PENDING_STATUS
    selection["pending_decisions"] = list(PENDING_SELECTION_DECISIONS)
    selection["rankings_are_separate_per_gene"] = True
    return ranking_rows, selection, sensitivity
