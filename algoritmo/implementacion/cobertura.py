"""Separate union-P completion from native search coverage by species."""

from collections import Counter, defaultdict


COMPLETE = {"positive_complete", "negative_complete"}
VALID_PAIRS = {"positive_witness", "negative_reported_pair"}
SOFT_UNAVAILABLE = {"unavailable_cap", "unavailable_ambiguity", "unexamined"}


def decision(pairs, required_pairs):
    """Decide a Boolean interaction flag without treating omissions as zero.

    Args:
        pairs: Unique final task results, not individual native attempts.
        required_pairs: Whole declared oriented query/reference universe.

    Returns:
        Positive-complete, negative-complete, unresolved or failed, and flag.
    """
    tasks = {(p["query_id"], p["reference_id"]): p for p in pairs}
    if len(tasks) > required_pairs:
        raise ValueError("Evaluated pairs exceed the declared scope")
    soft, positive = False, False
    for pair in tasks.values():
        state = pair["state"]
        if state not in VALID_PAIRS:
            if state in SOFT_UNAVAILABLE:
                soft = True
                continue
            return "failed", None
        expected = state == "positive_witness"
        if pair.get("flag") is not expected:
            return "failed", None
        if not expected and pair.get("greedy_complete") is not True:
            soft = True
        positive |= expected
    if soft:
        return "unresolved", None
    if positive:
        return "positive_complete", True
    if len(tasks) == required_pairs:
        return "negative_complete", False
    return "unresolved", None


def position_coverage(windows, pairs, datasets, producer_positions):
    """Build independently checked union flags and per-species observations.

    Args:
        windows: Original selected window mappings.
        pairs: Verified final task records, including greedy completeness.
        datasets: Frozen reference datasets and declared record counts.
        producer_positions: Original exported P positions keyed by window ID.

    Returns:
        Union position rows, species rows, and query membership by position.
    """
    by_query = defaultdict(list)
    for pair in pairs:
        by_query[pair["query_id"]].append(pair)
    positions, species, membership = [], [], {}
    total = sum(d["statistics"]["records"] for d in datasets)
    for window in windows:
        wid = window["window_id"]
        qids = {o["query_id"] for o in window["occurrences"]}
        membership[wid] = qids
        observed = {p["evidence_key"]: p for q in qids
                    for p in by_query[q]}
        tasks = list(observed.values())
        state, flag = decision(tasks, len(qids) * total)
        reported = producer_positions.get(wid)
        reason = None
        if reported is None:
            state, flag = "unresolved", None
            reason = "missing_producer_position"
        elif state in COMPLETE and (
                reported.get("status") != "measured"
                or reported.get("flagged") is not flag):
            state, flag, reason = "failed", None, "producer_flag_disagreement"
        positions.append({
            "window_id": wid, "parent_id": window["parent_id"],
            "candidate_id": window["parent_id"].split("@")[0],
            "interval_1based_inclusive": window["interval_1based_inclusive"],
            "status": state, "flagged": flag,
            "reported_flag": reported.get("flagged") if reported else None,
            "withholding_reason": reason,
            "required_pairs": len(qids) * total,
            "attempted_pairs": len(tasks),
            "evaluated_pairs": sum(p["state"] in VALID_PAIRS for p in tasks),
            "unattempted_pairs": len(qids) * total - len(tasks),
            "unattempted_disposition": (
                "not_required_for_union_P_after_positive"
                if state == "positive_complete" else "required_untested"),
        })
        for dataset in datasets:
            selected = [p for p in tasks
                        if p["dataset_id"] == dataset["id"]]
            required = len(qids) * dataset["statistics"]["records"]
            species_state, _ = decision(selected, required)
            if not selected:
                species_state = "untested"
            elif species_state == "unresolved":
                species_state = "partial_unresolved"
            species.append({
                "window_id": wid,
                "candidate_id": positions[-1]["candidate_id"],
                "dataset_id": dataset["id"], "organism": dataset["organism"],
                "status": species_state, "union_status": state,
                "required_pairs": required, "attempted_pairs": len(selected),
                "evaluated_pairs": sum(p["state"] in VALID_PAIRS
                                       for p in selected),
                "positive_pairs": sum(p["state"] == "positive_witness"
                                      for p in selected),
                "negative_reported_pairs": sum(
                    p["state"] == "negative_reported_pair" for p in selected),
                "greedy_traces_exhausted": sum(
                    p.get("greedy_complete") is True for p in selected),
                "unattempted_pairs": required - len(selected),
                "unattempted_disposition": positions[-1][
                    "unattempted_disposition"],
            })
    return positions, species, membership


def search_summary(window_ids, positions, species_positions, membership,
                   pairs, datasets, level, identifier):
    """Summarize a region/gene/panel with unique query/reference accounting.

    Args:
        window_ids: Set of original window IDs in this population.
        positions: Union position rows keyed by window ID.
        species_positions: Species position rows keyed by window/dataset.
        membership: Original query membership per window.
        pairs: Verified final pair records.
        datasets: Frozen required datasets.
        level: Region, gene or panel.
        identifier: Original region/gene ID or frozen panel ID.

    Returns:
        One coverage row for each species; unknown searches remain explicit.
    """
    qids = {q for wid in window_ids for q in membership[wid]}
    still_required_qids = {q for wid in window_ids
                           if positions[wid]["status"] != "positive_complete"
                           for q in membership[wid]}
    not_required_qids = qids - still_required_qids
    result = []
    for dataset in datasets:
        did = dataset["id"]
        tasks = {(p["query_id"], p["reference_id"]): p for p in pairs
                 if p["query_id"] in qids and p["dataset_id"] == did}
        count = dataset["statistics"]["records"]
        states = Counter(species_positions[wid, did]["status"]
                         for wid in window_ids)
        unattempted = len(qids) * count - len(tasks)
        unnecessary = len(not_required_qids) * count - sum(
            q in not_required_qids for q, _ in tasks)
        failed = sum(p["state"] not in VALID_PAIRS for p in tasks.values())
        result.append({
            "level": level, "population_id": identifier,
            "dataset_id": did, "organism": dataset["organism"],
            "reference_records": count, "unique_positions": len(window_ids),
            "unique_queries": len(qids), "required_pairs": len(qids) * count,
            "attempted_pairs": len(tasks),
            "evaluated_pairs": len(tasks) - failed,
            "unavailable_pairs": failed,
            "positive_pairs": sum(p["state"] == "positive_witness"
                                  for p in tasks.values()),
            "negative_reported_pairs": sum(
                p["state"] == "negative_reported_pair"
                for p in tasks.values()),
            "greedy_traces_exhausted": sum(
                p.get("greedy_complete") is True for p in tasks.values()),
            "unattempted_pairs": unattempted,
            "not_required_for_union_P_pairs": unnecessary,
            "still_required_untested_pairs": unattempted - unnecessary,
            "positive_complete_positions": states["positive_complete"],
            "negative_complete_positions": states["negative_complete"],
            "partial_unresolved_positions": states["partial_unresolved"],
            "untested_positions": states["untested"],
            "failed_positions": states["failed"],
            "all_required_pairs_evaluated": len(tasks) - failed
            == len(qids) * count,
            "every_position_boolean_decided_for_species": sum(
                states[s] for s in COMPLETE) == len(window_ids),
        })
    return result
