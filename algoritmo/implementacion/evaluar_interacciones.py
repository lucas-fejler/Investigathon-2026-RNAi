"""Evaluate interactions with bounded execution and explicit completeness."""

from collections import Counter
from datetime import datetime, timezone
from itertools import zip_longest
import os
from pathlib import Path
import signal
import subprocess
import time

import rnahybrid as native
from archivos import (digest, fasta, identity, read_json, read_rows,
                          sequence_sha, write_json)
from seleccion_panel import POLICY, ROOT, validate_manifest


def catalog(manifest):
    """Load exact reference identities without copying reference sequences.

    Args:
        manifest: Verified frozen panel.

    Returns:
        Reference bindings keyed by original identifier.
    """
    base = ROOT / Path(manifest["reference_scope"]["manifest"]["path"]).parent
    records = {}
    for dataset in manifest["reference_scope"]["datasets"]:
        path = base / dataset["records"]
        if digest(path) != dataset["records_sha256"]:
            raise ValueError("Changed reference catalog: " + dataset["id"])
        count = 0
        for row in read_rows(path):
            rid = row["reference_id"]
            if rid in records:
                raise ValueError("Duplicate reference identifier")
            records[rid] = {k: row[k] for k in (
                "reference_id", "sequence_sha256", "length")}
            records[rid]["organism"] = dataset["organism"]
            count += 1
        if count != dataset["statistics"]["records"]:
            raise ValueError("Incomplete required catalog")
    if len(records) != manifest["reference_scope"]["required_records"]:
        raise ValueError("Reference scope cannot be narrowed")
    return records


def references(manifest):
    """Stream and verify whole records against their curated identities.

    Args:
        manifest: Verified panel.

    Yields:
        Reference ID and complete stored-forward sequence.

    Returns:
        None.
    """
    base = ROOT / Path(manifest["reference_scope"]["manifest"]["path"]).parent
    for dataset in manifest["reference_scope"]["datasets"]:
        count = 0
        pairs = zip_longest(fasta(base / dataset["fasta"]),
                           read_rows(base / dataset["records"]))
        for record, meta in pairs:
            if record is None or meta is None:
                raise ValueError("Reference FASTA/catalog count mismatch")
            rid, seq = record
            if (rid != meta["reference_id"] or len(seq) != meta["length"]
                    or sequence_sha(seq) != meta["sequence_sha256"]):
                raise ValueError("Reference record identity mismatch: " + rid)
            count += 1
            yield rid, seq
        if count != dataset["statistics"]["records"]:
            raise ValueError("Reference count differs from declared scope")


def build_identity():
    """Bind the synchronized evaluator and its actual laptop libraries.

    Args:
        None.

    Returns:
        Portable evaluator identity and local runtime details.
    """
    root = ROOT / "tools/rnahybrid/2.1.2"
    build = read_json(root / "manifest.json")
    files = {}
    for name in ("audit-launcher", "audit/src/RNAhybrid", "audit.patch",
                 "stock-launcher", "stock/src/RNAhybrid"):
        sha = digest(root / name)
        if sha != build["files_sha256"][name]:
            raise ValueError("Synchronized evaluator changed: " + name)
        files[name] = sha
    env = {**os.environ, "LD_LIBRARY_PATH": str(root / "deps/usr/lib")}
    result = subprocess.run(["ldd", str(root / "audit/src/RNAhybrid")],
                            env=env, capture_output=True, text=True,
                            check=True)
    libraries = {}
    for line in result.stdout.splitlines():
        if "not found" in line:
            raise ValueError("Native dependency unavailable: " + line)
        for field in line.split():
            if field.startswith("/"):
                path = Path(field).resolve()
                libraries[path.name] = {"path": str(path),
                                        "sha256": digest(path)}
    portable = {"original_build_id": build["build_id"], "files": files,
                "runtime_libraries": {k: v["sha256"]
                                      for k, v in libraries.items()},
                "policy": POLICY, "parser_sha256": digest(
                    ROOT / "rnahybrid.py")}
    return {"evaluator_id": identity(portable), "binding": portable,
            "runtime": libraries, "python": os.sys.executable,
            "python_version": os.sys.version}


def prepare(output):
    """Record initial workload and prioritize evidence without excluding work.

    Args:
        output: Frozen run directory.

    Returns:
        Workload summary after valid evidence reuse.
    """
    output = Path(output)
    manifest = validate_manifest(output)
    records = catalog(manifest)
    build = build_identity()
    binding = {"panel_id": manifest["panel_id"],
               "evaluator_id": build["evaluator_id"],
               "classification_id": manifest["classification_id"]}
    if (output / "execution_binding.json").exists():
        if read_json(output / "execution_binding.json") != binding:
            raise ValueError("Cannot reuse a different evaluator or panel")
    else:
        write_json(output / "execution_binding.json", binding)
        write_json(output / "environment.json", build)
    queries = list(read_rows(output / "selected_queries.jsonl"))
    pairs = len(queries) * len(records)
    summary = {**binding, "unique_query_reference_pairs": pairs,
               "valid_prior_completed_pairs_reused": 0,
               "initial_outstanding_unique_pairs": pairs,
               "unique_queries": len(queries),
               "reference_records": len(records),
               "prior_evidence_disposition": "Native raw invocation files "
               "absent; retained original summaries are historical only.",
               "completion_time_estimate": None,
               "stopping": "positive Boolean witness may make other pairs "
               "not required; all undecided positions retain full fallback"}
    if not (output / "workload.json").exists():
        write_json(output / "workload.json", summary)
    priority = output / "priority.tsv"
    stamp = output / "priority_identity.json"
    required = {**binding,
                "scanner_sha256": digest(ROOT / "priorizar_referencias.cpp"),
                "query_sha256": digest(output / "selected_queries.jsonl")}
    if stamp.exists():
        saved = read_json(stamp)
        if (saved["binding"] != required
                or saved["priority_sha256"] != digest(priority)):
            raise ValueError("Stale heuristic priority evidence")
        return summary
    binary = output / "local/priority"
    binary.parent.mkdir(exist_ok=True)
    subprocess.run(["g++", "-O3", "-std=c++17",
                    str(ROOT / "priorizar_referencias.cpp"), "-o", str(binary)],
                   check=True)
    qfile = output / "local/queries.tsv"
    qfile.write_text("".join(q["query_id"] + "\t" + q["sequence"] + "\n"
                             for q in queries))
    started = time.monotonic()
    with priority.with_suffix(".tmp").open("w") as handle:
        process = subprocess.Popen([str(binary), str(qfile)],
                                   stdin=subprocess.PIPE, stdout=handle,
                                   text=True)
        try:
            for rid, seq in references(manifest):
                process.stdin.write(rid + "\t" + seq + "\n")
            process.stdin.close()
            if process.wait() != 0:
                raise ValueError("Priority scanner failed")
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    priority.with_suffix(".tmp").replace(priority)
    write_json(stamp, {"binding": required,
                       "priority_sha256": digest(priority),
                       "records_verified": len(records),
                       "wall_seconds": time.monotonic() - started,
                       "assurance": "heuristic_priority_only_not_exclusion"})
    return summary


def pair_key(query, reference, binding):
    """Bind a reusable task to query, whole reference and scientific policy.

    Args:
        query: Original M3 query row.
        reference: Whole-reference identity.
        binding: Frozen panel/evaluator/classification identity.

    Returns:
        Content-addressed task key.
    """
    return identity({"query": query, "reference": reference,
                     "binding": binding})


def native_options(source, cap):
    """Specify the pinned options shared by launch and evidence validation.

    Args:
        source: Exact query and whole-reference ledger.
        cap: Five-hit prefix or verified finite greedy-report ceiling.

    Returns:
        Ordered native command options.
    """
    length = source["target_length"]
    if cap not in (5, native.report_ceiling(length)):
        raise ValueError("Unsupported native cap")
    return ["-c", "-b", str(cap), "-e", "-10", "-p", "1", "-d", "0,1",
            "-u", "15", "-v", "15", "-m", str(length), "-n", "21",
            "-t", "targets.fasta", "-q", "queries.fasta"]


def parse_call(directory, source, config, candidate):
    """Reparse raw audited evidence and classify every qualifying block.

    Args:
        directory: One immutable invocation directory.
        source: Exact whole input ledger.
        config: Original curated config.
        candidate: Candidate-specific intended exclusions.

    Returns:
        Verified pair decision, blocks and completeness.
    """
    directory = Path(directory)
    call = read_json(directory / "call.json")
    if call["status"] != "completed" or call["exit_code"] != 0:
        return {"state": call["status"], "blocks": [], "flag": None}
    if (call["options"] != native_options(source, call["cap"])
            or read_json(directory / "source.json") != source):
        raise ValueError("Cached evaluator options/input ledger mismatch")
    expected = {("Q0000", "T0000"): source}
    for name, sha in call["artifacts_sha256"].items():
        if digest(directory / name) != sha:
            raise ValueError("Changed native evidence artifact: " + name)
    for name, value in (("queries.fasta", source["query_sequence"]),
                        ("targets.fasta", source["target_sequence"])):
        rows = list(fasta(directory / name))
        if len(rows) != 1 or rows[0][1] != value:
            raise ValueError("Native FASTA input identity mismatch")
    if (directory / "stderr.txt").read_bytes():
        raise ValueError("Native stderr is not empty")
    audit = native.parse_audit(
        (directory / "audit.tsv").read_bytes(), expected)
    if audit["errors"]:
        raise ValueError("; ".join(audit["errors"]))
    profile = {"min_align_length": 18, "energy_cutoff": -10,
               "pvalue_cutoff": 1}
    raw = native.parse_compact((directory / "stdout.txt").read_bytes(),
                               expected, profile)
    if len(raw) > call["cap"]:
        raise ValueError("Native cap exceeded")
    blocks = native.mapped_blocks({"blocks": raw}, source)
    native.check_progress(blocks, source["target_length"])
    for block in blocks:
        block["classification"] = native.classify(
            candidate, source["reference_id"],
            block["core_target_interval"], config)
    witnesses = [b for b in blocks if b["draft_qualifies"]
                 and b["classification"].startswith("unintended")]
    if witnesses:
        state, flag = "positive_witness", True
    elif set(source["target_sequence"]) - set("ACGT"):
        state, flag = "unavailable_ambiguity", None
    elif len(raw) == call["cap"]:
        state, flag = "unavailable_cap", None
    else:
        state, flag = "negative_reported_pair", False
    return {"state": state, "flag": flag, "blocks": blocks,
            "greedy_complete": len(raw) < call["cap"],
            "witnesses": [b["duplex_id"] for b in witnesses]}


def invoke(directory, source, cap, timeout):
    """Checkpoint before launch and preserve failure or completed raw output.

    Args:
        directory: New per-call directory.
        source: Exact query/whole-reference ledger.
        cap: Native report cap, never used as a completed negative.
        timeout: Bound for this invocation in seconds.

    Returns:
        None. Evidence is persisted even after failure.
    """
    directory = Path(directory).resolve()
    directory.mkdir(parents=True)
    (directory / "queries.fasta").write_bytes(native.fasta_bytes(
        [("Q0000", source["query_sequence"])]))
    (directory / "targets.fasta").write_bytes(native.fasta_bytes(
        [("T0000", source["target_sequence"])]))
    write_json(directory / "source.json", source)
    options = native_options(source, cap)
    command = [str(ROOT / "tools/rnahybrid/2.1.2/audit-launcher"), *options]
    call = {"status": "interrupted_or_running", "exit_code": None,
            "cap": cap, "options": options, "timeout_seconds": timeout,
            "started_utc": datetime.now(timezone.utc).isoformat()}
    write_json(directory / "call.json", call)
    env = {"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
           "RNAHYBRID_AUDIT": str(directory / "audit.tsv")}
    started = time.monotonic()
    (directory / "audit.tsv").touch()
    with (directory / "stdout.txt").open("wb") as stdout, (
            directory / "stderr.txt").open("wb") as stderr:
        process = subprocess.Popen(command, cwd=directory, env=env,
                                   stdout=stdout, stderr=stderr,
                                   start_new_session=True)
        try:
            call["exit_code"] = process.wait(timeout=timeout)
            call["status"] = ("completed" if process.returncode == 0
                              else "failed")
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            call["status"] = "unavailable_timeout"
        except BaseException:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            call["status"] = "unavailable_interrupted"
            write_json(directory / "call.json", call)
            raise
    call["wall_seconds"] = time.monotonic() - started
    call["artifacts_sha256"] = {name: digest(directory / name) for name in (
        "stdout.txt", "stderr.txt", "audit.tsv", "queries.fasta",
        "targets.fasta", "source.json")}
    write_json(directory / "call.json", call)


def assess_pair(output, query, reference, sequence, candidate, binding,
                config, timeout=120, execute=True, repairing=False):
    """Reuse verified evidence or run a cap-safe native pair with checkpoint.

    Args:
        output: Frozen run directory.
        query: Original query row.
        reference: Original whole-reference identity.
        sequence: Complete stored-forward reference sequence.
        candidate: Candidate-specific classifier.
        binding: Exact policy and panel identities.
        config: Original curation config.
        timeout: Maximum seconds per invocation.
        execute: Permit new native work; false for export and resume audit.
        repairing: Allow one explicitly recorded same-pair recovery attempt.

    Returns:
        Revalidated decision and evidence key.
    """
    key = pair_key(query, reference, binding)
    if candidate not in query["candidate_ids"]:
        raise ValueError("Classifier candidate does not own this query")
    directory = Path(output) / "evidence" / key
    expected_binding = {"query": query, "reference": reference,
                        "binding": binding, "candidate": candidate}
    if directory.exists():
        if read_json(directory / "binding.json") != expected_binding:
            raise ValueError("Cached task/classifier identity mismatch")
    else:
        directory.mkdir(parents=True)
        write_json(directory / "binding.json", expected_binding)
    source = {"query_id": query["query_id"],
              "query_sequence": query["sequence"],
              "query_sha256": query["sequence_sha256"], "query_length": 21,
              "reference_id": reference["reference_id"],
              "reference_sha256": reference["sequence_sha256"],
              "target_length": reference["length"],
              "target_sequence": sequence}
    if (sequence_sha(sequence) != reference["sequence_sha256"]
            or len(sequence) != reference["length"]):
        raise ValueError("Reference identity does not match native task")
    recovery_error = check_recovery(directory, expected_binding, repairing)
    previous = None
    for cap in (5, native.report_ceiling(len(sequence))):
        if recovery_error:
            result = {"state": "unavailable_required_recovery", "flag": None,
                      "blocks": [], "error": recovery_error}
            break
        arm = directory / str(cap)
        if not arm.exists():
            if not execute:
                result = {"state": "unavailable_interrupted", "flag": None,
                          "blocks": []}
                break
            try:
                invoke(arm, source, cap, timeout)
            except OSError as exc:
                result = {"state": "unavailable_launch", "flag": None,
                          "blocks": [], "error": str(exc)}
                break
        if not (arm / "call.json").exists():
            result = {"state": "unavailable_interrupted", "flag": None,
                      "blocks": []}
        else:
            try:
                result = parse_call(arm, source, config, candidate)
            except (ValueError, KeyError, OSError) as exc:
                result = {"state": "unavailable_verification", "flag": None,
                          "blocks": [], "error": str(exc)}
        write_json(arm / "assessment.json", {
            "state": result["state"], "flag": result["flag"],
            "reported_blocks": len(result["blocks"]),
            "greedy_complete": result.get("greedy_complete"),
            "error": result.get("error")})
        if previous is not None:
            prefix = [b["raw_block"] for b in previous["blocks"]]
            new_prefix = [b["raw_block"]
                          for b in result["blocks"][:len(prefix)]]
            if prefix != new_prefix:
                result = {"state": "unavailable_prefix_mismatch", "flag": None,
                          "blocks": result["blocks"]}
        if result["state"] != "unavailable_cap":
            break
        previous = result
    result.update({"evidence_key": key, "query_id": query["query_id"],
                   "reference_id": reference["reference_id"],
                   "organism": reference["organism"],
                   "evidence_binding": expected_binding})
    write_json(directory / "result.json", result)
    return result


def check_recovery(directory, expected_binding, repairing=False):
    """Keep failed attempts unavailable until the exact task is completed.

    Args:
        directory: One query/reference evidence directory.
        expected_binding: Exact input, evaluator and classification binding.
        repairing: Permit a pending recorded retry during its own execution.

    Returns:
        An explicit blocking reason, or None for verified ordinary/repaired work.
    """
    path = directory / "recovery.json"
    if not path.exists():
        return ("Missing required recovery provenance" if list(
            directory.glob("*-failed-attempt")) else None)
    recovery = read_json(path)
    if recovery["binding"] != expected_binding:
        return "Recovery task binding changed"
    for name, sha in recovery["archived_sha256"].items():
        if (not (directory / name).is_file()
                or digest(directory / name) != sha):
            return "Required failed-attempt evidence changed"
    if repairing and recovery["status"] == "pending_exact_task_retry":
        return None
    if recovery["status"] != "completed_exact_task_retry":
        return "Required task retry has not completed"
    for name, sha in recovery["replacement_calls_sha256"].items():
        if digest(directory / name) != sha:
            return "Required replacement evidence changed"
    return None


def retry_pair(output, query, reference, sequence, candidate, binding, config):
    """Explicitly retry an execution failure, retaining every original byte.

    Args:
        output: Frozen run directory; its runner must already be stopped.
        query: Original oriented query.
        reference: Original whole reference.
        sequence: Exact whole-reference sequence.
        candidate: Original classifier candidate.
        binding: Unchanged panel/evaluator/classification binding.
        config: Original candidate registry.

    Returns:
        Verified replacement decision, or unavailable after a failed retry.
    """
    key = pair_key(query, reference, binding)
    directory = Path(output) / "evidence" / key
    expected = {"query": query, "reference": reference,
                "candidate": candidate, "binding": binding}
    if read_json(directory / "binding.json") != expected:
        raise ValueError("Cannot retry a different scientific task")
    if (sequence_sha(sequence) != reference["sequence_sha256"]
            or len(sequence) != reference["length"]):
        raise ValueError("Cannot retry changed reference bytes")
    if (directory / "recovery.json").exists():
        raise ValueError("A recovery is already recorded; inspect it first")
    failures = []
    valid = {"failed", "unavailable_timeout", "unavailable_interrupted",
             "interrupted_or_running", "unavailable_launch"}
    for cap in (5, native.report_ceiling(reference["length"])):
        arm = directory / str(cap)
        if (arm / "call.json").exists() and read_json(
                arm / "call.json")["status"] in valid:
            failures.append(arm)
    if not failures:
        raise ValueError("No execution failure eligible for explicit retry")
    recovery = {"binding": expected, "status": "pending_exact_task_retry",
                "archived_sha256": {}, "replacement_calls_sha256": {}}
    moves = []
    for arm in failures:
        archive = directory / (arm.name + "-failed-attempt")
        if archive.exists():
            raise ValueError("Existing failure archive; refusing overwrite")
        recovery["archived_sha256"].update({
            str(archive.relative_to(directory) / p.name): digest(p)
            for p in arm.iterdir() if p.is_file()})
        moves.append((arm, archive))
    write_json(directory / "recovery.json", recovery)
    for arm, archive in moves:
        arm.rename(archive)
    result = assess_pair(output, query, reference, sequence, candidate,
                         binding, config, repairing=True)
    if result["state"] in ("positive_witness", "negative_reported_pair"):
        recovery["status"] = "completed_exact_task_retry"
        recovery["replacement_calls_sha256"] = {
            str(p.relative_to(directory)): digest(p) for p in
            directory.glob("*/call.json") if p.parent.name.isdigit()}
    write_json(directory / "recovery.json", recovery)
    return assess_pair(output, query, reference, sequence, candidate,
                       binding, config, execute=False)


def repair(output, maximum=8):
    """Complete stopped execution failures without discarding their history.

    Args:
        output: Frozen run directory, with no active assessment process.
        maximum: Hard maximum same-pair recoveries in this invocation.

    Returns:
        Explicit recovery states, never assumed no-hit measurements.
    """
    output = Path(output)
    manifest = validate_manifest(output)
    binding = read_json(output / "execution_binding.json")
    if build_identity()["evaluator_id"] != binding["evaluator_id"]:
        raise ValueError("Evaluator changed before recovery")
    records = catalog(manifest)
    queries = {q["query_id"]: q for q in
               read_rows(output / "selected_queries.jsonl")}
    config = read_json(ROOT / "configuracion.json")
    recovered = []
    for directory in sorted((output / "evidence").glob("*")):
        failed = [p for p in directory.glob("*/call.json")
                  if p.parent.name.isdigit() and read_json(p)["status"]
                  not in ("completed",)]
        if not failed:
            continue
        if len(recovered) >= maximum:
            break
        saved = read_json(directory / "binding.json")
        if saved["binding"] != binding:
            raise ValueError("Changed failed-task policy")
        qid = saved["query"]["query_id"]
        rid = saved["reference"]["reference_id"]
        source = read_json(failed[0].parent / "source.json")
        result = retry_pair(output, queries[qid], records[rid],
                            source["target_sequence"], saved["candidate"],
                            binding, config)
        recovered.append({"evidence_key": directory.name,
                          "state": result["state"]})
    write_json(output / "recoveries.json", recovered)
    return recovered


def position_state(results, total):
    """Apply Boolean OR without converting failed/omitted work into zero.

    Args:
        results: Verified attempted pair decisions across both orientations.
        total: Required query/reference tasks for this position.

    Returns:
        Explicit position state and nullable flag.
    """
    blockers = [r for r in results if r["state"] not in (
        "positive_witness", "negative_reported_pair")]
    if blockers:
        return "unavailable_required_dependency", None
    if any(r["state"] == "positive_witness" for r in results):
        return "positive_complete", True
    if len(results) == total and all(r["flag"] is False for r in results):
        return "negative_complete", False
    return "untested_or_incomplete", None


def priority_tasks(window, choices, records, config):
    """Prefer nonmember references without excluding any intended record.

    Args:
        window: Original position and both oriented query occurrences.
        choices: Heuristic whole-reference suggestions per query.
        records: Complete required reference catalog.
        config: Candidate-specific intended-feature registry.

    Returns:
        All suggested tasks, ordered to seek an unintended witness first.
    """
    gene = window["parent_id"].split("@")[0]
    intended = {r["reference_id"] for r in config["candidates"][gene][
        "curation"]["intended_targets"]}
    tasks = [(qid, rid, score)
             for qid in {o["query_id"] for o in window["occurrences"]}
             for rid, score in choices[qid]]
    return sorted(tasks, key=lambda item: (
        item[1] in intended, -item[2], records[item[1]]["length"],
        item[1], item[0]))


def run(output, wall_seconds=3600, max_pairs=20000):
    """Run simple priority rounds then exhaustive fallback for undecided work.

    Args:
        output: Prepared frozen panel directory.
        wall_seconds: Wall budget for this resumable execution segment.
        max_pairs: Attempted pair ceiling for this segment.

    Returns:
        Execution summary with complete/incomplete position counts.
    """
    output = Path(output)
    manifest = validate_manifest(output)
    binding = read_json(output / "execution_binding.json")
    priority_stamp = read_json(output / "priority_identity.json")
    if priority_stamp["priority_sha256"] != digest(output / "priority.tsv"):
        raise ValueError("Changed priority artifact")
    if build_identity()["evaluator_id"] != binding["evaluator_id"]:
        raise ValueError("Changed native evaluator")
    records = catalog(manifest)
    queries = {r["query_id"]: r for r in
               read_rows(output / "selected_queries.jsonl")}
    windows = list(read_rows(output / "selected_windows.jsonl"))
    config = read_json(ROOT / "configuracion.json")
    choices = {qid: [] for qid in queries}
    for line in (output / "priority.tsv").read_text().splitlines():
        qid, rid, score = line.split("\t")
        choices[qid].append((rid, int(score)))
    cached = [(directory, read_json(directory / "binding.json"))
              for directory in sorted((output / "evidence").glob("*"))]
    # A resume loads sparse priority and prior fallback inputs in one pass.
    needed = {rid for rows in choices.values() for rid, score in rows}
    needed.update(saved["reference"]["reference_id"] for _, saved in cached)
    sequences = {rid: seq for rid, seq in references(manifest)
                 if rid in needed}
    results = {w["window_id"]: [] for w in windows}
    attempted = set()
    for directory, saved in cached:
        qid = saved["query"]["query_id"]
        rid = saved["reference"]["reference_id"]
        if saved["binding"] != binding or qid not in queries:
            raise ValueError("Unrelated cached evidence in run directory")
        result = assess_pair(
            output, queries[qid], records[rid], sequences[rid],
            saved["candidate"], binding, config, execute=False)
        attempted.add((qid, rid))
        for window in windows:
            if any(o["query_id"] == qid for o in window["occurrences"]):
                results[window["window_id"]].append(result)
    started, count = time.monotonic(), 0
    last_save = started
    schedule = {}
    for window in windows:
        schedule[window["window_id"]] = priority_tasks(
            window, choices, records, config)

    def checkpoint(stage):
        """Persist verified progress and run accounting.

        Args:
            stage: Current scheduling phase.

        Returns:
            Summary.
        """
        rows = []
        for window in windows:
            wid = window["window_id"]
            state, flag = position_state(results[wid], len(records) * 2)
            rows.append({"window_id": wid, "state": state, "flagged": flag,
                         "pairs": [{k: r.get(k) for k in (
                             "evidence_key", "state", "query_id",
                             "reference_id",
                             "organism", "witnesses", "error")}
                                   for r in results[wid]],
                         "unattempted_pair_state": (
                             "not_required_after_positive"
                             if state == "positive_complete"
                             else "required_untested"),
                         "required_pairs": len(records) * 2,
                         "attempted_pairs": len(results[wid])})
        from archivos import write_rows
        write_rows(output / "position_checkpoints.jsonl", rows)
        summary = {"stage": stage, "positions": dict(Counter(
            r["state"] for r in rows)), "segment_pairs": count,
            "segment_wall_seconds": time.monotonic() - started,
            "total_unique_attempted_pairs": len(attempted),
            "panel_id": manifest["panel_id"], "binding": binding}
        write_json(output / "execution.json", summary)
        print(summary, flush=True)
        return summary

    def attempt(window, qid, rid, sequence):
        """Execute one still-required task within segment limits.

        Args:
            window: Original position mapping.
            qid: Oriented query.
            rid: Whole reference.
            sequence: Exact reference bytes.

        Returns:
            Whether the segment budget allows continuing.
        """
        nonlocal count, last_save
        elapsed = time.monotonic() - started
        if elapsed >= wall_seconds or count >= max_pairs:
            return False
        if (qid, rid) in attempted:
            return True
        gene = window["parent_id"].split("@")[0]
        result = assess_pair(output, queries[qid], records[rid], sequence,
                             gene, binding, config,
                             timeout=min(120, max(1, wall_seconds - elapsed)))
        attempted.add((qid, rid))
        count += 1
        for other in windows:
            if any(o["query_id"] == qid for o in other["occurrences"]):
                results[other["window_id"]].append(result)
        if time.monotonic() - last_save >= 30:
            checkpoint("running")
            last_save = time.monotonic()
        return True

    def unsettled(window):
        """Check whether more native work can settle this position.

        Args:
            window: Original window mapping.

        Returns:
            True only for incomplete coverage without a prior failed task.
        """
        return position_state(results[window["window_id"]], len(records) * 2)[
            0] == "untested_or_incomplete"

    checkpoint("resumed")
    for depth in range(48):
        for window in windows:
            tasks = schedule[window["window_id"]]
            if not unsettled(window) or depth >= len(tasks):
                continue
            qid, rid, _ = tasks[depth]
            if not attempt(window, qid, rid, sequences[rid]):
                return checkpoint("segment_budget_exhausted")
        if not any(unsettled(w) for w in windows):
            return checkpoint("assessment_finished")
    checkpoint("fallback_required")
    for rid, seq in references(manifest):
        for window in windows:
            if not unsettled(window):
                continue
            for occ in window["occurrences"]:
                if not unsettled(window):
                    break
                if not attempt(window, occ["query_id"], rid, seq):
                    return checkpoint("segment_budget_exhausted")
        if not any(unsettled(w) for w in windows):
            break
    return checkpoint("assessment_finished")
