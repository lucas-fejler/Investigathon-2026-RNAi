"""Read scientific artifacts and persist deterministic identities."""

import gzip
import hashlib
import json
from pathlib import Path


def identity(value):
    """Hash a structured deterministic payload.

    Args:
        value: JSON-compatible payload.

    Returns:
        SHA-256 of canonical JSON bytes.
    """
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode()).hexdigest()


def digest(path):
    """Hash one required artifact.

    Args:
        path: Existing file.

    Returns:
        Hexadecimal SHA-256.
    """
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_json(path):
    """Read a JSON checkpoint.

    Args:
        path: JSON file.

    Returns:
        Parsed value.
    """
    return json.loads(Path(path).read_text())


def read_rows(path):
    """Read plain or gzipped JSONL records.

    Args:
        path: JSONL file.

    Yields:
        Parsed rows in source order.

    Returns:
        None.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            yield json.loads(line)


def write_json(path, value):
    """Replace a derived checkpoint atomically.

    Args:
        path: Destination, never a source artifact.
        value: JSON payload.

    Returns:
        None.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_rows(path, rows):
    """Write a derived JSONL export atomically.

    Args:
        path: Destination file.
        rows: JSON-compatible records.

    Returns:
        None.
    """
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    temporary.replace(path)


def fasta(path):
    """Stream full records without altering sequence bytes or orientation.

    Args:
        path: Plain or compressed FASTA.

    Yields:
        Original first-token identifier and concatenated sequence.

    Returns:
        None.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def sequence_sha(sequence):
    """Hash exact sequence bytes.

    Args:
        sequence: Preserved ASCII sequence.

    Returns:
        SHA-256.
    """
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def reverse_complement(sequence):
    """Return the opposite oriented query without changing the parent.

    Args:
        sequence: ACGT DNA.

    Returns:
        Reverse-complement DNA.
    """
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
