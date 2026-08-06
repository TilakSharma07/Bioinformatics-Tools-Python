"""Format validation with actionable error messages.

Every check reports the record ID and line number, so a failure in a
10-million-record file points at the offending record instead of just
saying "invalid FASTA".
"""
from __future__ import annotations

from pathlib import Path

from biotools.fasta import read_fasta
from biotools.fastq import PHRED_OFFSET, read_fastq

VALID_DNA = set("ACGTNRYKMSWBDHV-")
VALID_PROTEIN = set("ACDEFGHIKLMNPQRSTVWY*X-")


class ValidationError(Exception):
    """Raised when a file violates the format spec."""


def validate_fasta(path: str | Path, alphabet: str = "dna", strict: bool = True) -> dict:
    """Check a FASTA file. Returns a report dict; raises if strict and errors exist."""
    allowed = VALID_DNA if alphabet == "dna" else VALID_PROTEIN
    errors: list[str] = []
    seen: set[str] = set()
    n, total_len = 0, 0

    for rec in read_fasta(path):
        n += 1
        total_len += len(rec)
        if not rec.id:
            errors.append(f"record {n}: empty ID")
        if rec.id in seen:
            errors.append(f"record {n} ({rec.id}): duplicate ID")
        seen.add(rec.id)
        if not rec.sequence:
            errors.append(f"record {n} ({rec.id}): empty sequence")
        bad = set(rec.sequence.upper()) - allowed
        if bad:
            errors.append(
                f"record {n} ({rec.id}): invalid {alphabet} character(s) "
                f"{sorted(bad)[:5]}"
            )

    report = {
        "path": str(path),
        "records": n,
        "total_length": total_len,
        "errors": errors,
        "valid": not errors,
    }
    if strict and errors:
        raise ValidationError(f"{path}: {len(errors)} error(s); first: {errors[0]}")
    return report


def validate_fastq(path: str | Path, strict: bool = True) -> dict:
    """Check a FASTQ file: record structure, quality range, duplicate IDs."""
    errors: list[str] = []
    seen: set[str] = set()
    n, total_len, min_q, max_q = 0, 0, 127, -1

    try:
        for rec in read_fastq(path):
            n += 1
            total_len += len(rec)
            if rec.id in seen:
                errors.append(f"record {n} ({rec.id}): duplicate ID")
            seen.add(rec.id)
            for ch in rec.quality:
                q = ord(ch) - PHRED_OFFSET
                min_q, max_q = min(min_q, q), max(max_q, q)
                if not 0 <= q <= 93:
                    errors.append(
                        f"record {n} ({rec.id}): Phred {q} outside 0-93 "
                        f"(wrong quality encoding?)"
                    )
                    break
    except ValueError as exc:
        errors.append(str(exc))

    report = {
        "path": str(path),
        "records": n,
        "total_length": total_len,
        "phred_min": min_q if n else None,
        "phred_max": max_q if n else None,
        "errors": errors,
        "valid": not errors,
    }
    if strict and errors:
        raise ValidationError(f"{path}: {len(errors)} error(s); first: {errors[0]}")
    return report
