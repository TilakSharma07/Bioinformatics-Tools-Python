"""Command-line entry point: `biotools <command> <file>`."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from biotools import __version__
from biotools.fasta import read_fasta, write_fasta
from biotools.fastq import mean_quality, read_fastq
from biotools.seqstats import gc_content, n50, nucleotide_counts, translate
from biotools.validate import ValidationError, validate_fasta, validate_fastq


def cmd_stats(args: argparse.Namespace) -> int:
    """Per-record and summary statistics for a FASTA file."""
    lengths, rows = [], []
    for rec in read_fasta(args.input):
        lengths.append(len(rec))
        rows.append({"id": rec.id, "length": len(rec), "gc_percent": round(gc_content(rec.sequence), 2)})

    if not rows:
        print("No records found.", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps({"records": rows, "summary": _summary(lengths, rows)}, indent=2))
    else:
        print(f"{'id':<28}{'length':>12}{'GC%':>8}")
        print("-" * 48)
        for r in rows[: args.limit]:
            print(f"{r['id'][:27]:<28}{r['length']:>12,}{r['gc_percent']:>8.2f}")
        if len(rows) > args.limit:
            print(f"... and {len(rows) - args.limit:,} more (use --limit 0 for all)")
        s = _summary(lengths, rows)
        print("-" * 48)
        for k, v in s.items():
            print(f"{k:<28}{v:>20,}" if isinstance(v, int) else f"{k:<28}{v:>20}")
    return 0


def _summary(lengths: list[int], rows: list[dict]) -> dict:
    total = sum(lengths)
    return {
        "records": len(lengths),
        "total_bp": total,
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": round(total / len(lengths), 1),
        "n50": n50(lengths),
        "mean_gc_percent": round(sum(r["gc_percent"] for r in rows) / len(rows), 2),
    }


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate a FASTA or FASTQ file; exit 1 if invalid."""
    fn = validate_fastq if args.format == "fastq" else validate_fasta
    kwargs = {} if args.format == "fastq" else {"alphabet": args.alphabet}
    try:
        report = fn(args.input, strict=False, **kwargs)
    except ValidationError as exc:
        print(f"FAIL {exc}", file=sys.stderr)
        return 1

    print(json.dumps(report, indent=2) if args.json else _human_report(report))
    return 0 if report["valid"] else 1


def _human_report(report: dict) -> str:
    lines = [
        f"file      {report['path']}",
        f"records   {report['records']:,}",
        f"total bp  {report['total_length']:,}",
        f"status    {'VALID' if report['valid'] else 'INVALID'}",
    ]
    if report["errors"]:
        lines.append(f"errors    {len(report['errors'])}")
        lines += [f"  - {e}" for e in report["errors"][:10]]
        if len(report["errors"]) > 10:
            lines.append(f"  ... and {len(report['errors']) - 10} more")
    return "\n".join(lines)


def cmd_filter(args: argparse.Namespace) -> int:
    """Write out only the records passing length / GC thresholds."""
    kept, seen = [], 0
    for rec in read_fasta(args.input):
        seen += 1
        if len(rec) < args.min_length:
            continue
        gc = gc_content(rec.sequence)
        if not (args.min_gc <= gc <= args.max_gc):
            continue
        kept.append(rec)

    n = write_fasta(kept, args.output)
    pct = 100 * n / seen if seen else 0
    print(f"kept {n:,} / {seen:,} records ({pct:.1f}%) -> {args.output}")
    return 0


def cmd_translate(args: argparse.Namespace) -> int:
    """Six-frame or single-frame translation of every record."""
    frames = [0, 1, 2] if args.all_frames else [args.frame]
    out = []
    from biotools.fasta import FastaRecord
    from biotools.seqstats import reverse_complement

    for rec in read_fasta(args.input):
        for fr in frames:
            out.append(FastaRecord(f"{rec.id}_frame+{fr + 1}", rec.description,
                                   translate(rec.sequence, fr, to_stop=args.to_stop)))
            if args.all_frames:
                rc = reverse_complement(rec.sequence)
                out.append(FastaRecord(f"{rec.id}_frame-{fr + 1}", rec.description,
                                       translate(rc, fr, to_stop=args.to_stop)))
    n = write_fasta(out, args.output)
    print(f"wrote {n:,} translated records -> {args.output}")
    return 0


def cmd_fastq_stats(args: argparse.Namespace) -> int:
    """Read-length and quality summary for a FASTQ file."""
    lengths, quals, n_pass = [], [], 0
    for rec in read_fastq(args.input):
        lengths.append(len(rec))
        q = mean_quality(rec)
        quals.append(q)
        n_pass += q >= args.min_quality
        if args.limit and len(lengths) >= args.limit:
            break

    if not lengths:
        print("No reads found.", file=sys.stderr)
        return 1

    print(f"reads              {len(lengths):,}")
    print(f"total bases        {sum(lengths):,}")
    print(f"read length        {min(lengths)}-{max(lengths)} "
          f"(mean {sum(lengths) / len(lengths):.1f})")
    print(f"mean quality       Q{sum(quals) / len(quals):.2f}")
    print(f"reads >= Q{args.min_quality:<8} {n_pass:,} ({100 * n_pass / len(lengths):.1f}%)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="biotools",
                                description="FASTA/FASTQ utilities (no heavy dependencies).")
    p.add_argument("--version", action="version", version=f"biotools {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("stats", help="per-record + summary stats for FASTA")
    s.add_argument("input", type=Path)
    s.add_argument("--json", action="store_true")
    s.add_argument("--limit", type=int, default=20, help="rows to print (0 = all)")
    s.set_defaults(func=cmd_stats)

    v = sub.add_parser("validate", help="check file format integrity")
    v.add_argument("input", type=Path)
    v.add_argument("--format", choices=["fasta", "fastq"], default="fasta")
    v.add_argument("--alphabet", choices=["dna", "protein"], default="dna")
    v.add_argument("--json", action="store_true")
    v.set_defaults(func=cmd_validate)

    fl = sub.add_parser("filter", help="filter records by length and GC")
    fl.add_argument("input", type=Path)
    fl.add_argument("-o", "--output", type=Path, required=True)
    fl.add_argument("--min-length", type=int, default=0)
    fl.add_argument("--min-gc", type=float, default=0.0)
    fl.add_argument("--max-gc", type=float, default=100.0)
    fl.set_defaults(func=cmd_filter)

    t = sub.add_parser("translate", help="DNA -> protein")
    t.add_argument("input", type=Path)
    t.add_argument("-o", "--output", type=Path, required=True)
    t.add_argument("--frame", type=int, default=0, choices=[0, 1, 2])
    t.add_argument("--all-frames", action="store_true", help="all six reading frames")
    t.add_argument("--to-stop", action="store_true", help="stop at first stop codon")
    t.set_defaults(func=cmd_translate)

    q = sub.add_parser("fastq-stats", help="read length + quality summary")
    q.add_argument("input", type=Path)
    q.add_argument("--min-quality", type=float, default=20.0)
    q.add_argument("--limit", type=int, default=0, help="stop after N reads (0 = all)")
    q.set_defaults(func=cmd_fastq_stats)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "limit", None) == 0 and args.command == "stats":
        args.limit = 10**9
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError, ValidationError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
