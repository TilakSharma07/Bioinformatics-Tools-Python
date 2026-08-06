# biotools

FASTA/FASTQ parsing, validation and sequence statistics with **no third-party
dependencies** — standard library only. Importable as a package and usable as a
CLI.

```bash
pip install -e .
biotools stats examples/human_transcripts.fasta
```

```
id                                length     GC%
------------------------------------------------
NM_007294.4                        7,088   41.77
NM_000546.6                        2,512   53.38
NM_005228.5                        9,905   47.78
NM_000207.3                          465   63.87
NM_000518.5                          628   51.27
------------------------------------------------
records                                        5
total_bp                                  20,598
n50                                        7,088
mean_gc_percent                            51.61
```

(Real BRCA1, TP53, EGFR, INS and HBB transcripts fetched from NCBI by
`examples/fetch_examples.sh`.)

## Why no dependencies

Biopython is excellent and I use it elsewhere. This package exists for the case
where a parser has to run inside a container, a cluster job, or a CI step
without a package install — and where the parsing has to stream rather than
load a multi-gigabyte genome into RAM.

## Commands

| Command | Purpose |
|---|---|
| `biotools stats FILE` | per-record length + GC, plus N50 and totals |
| `biotools validate FILE` | format integrity; exits 1 on invalid input |
| `biotools filter FILE -o OUT` | keep records passing length / GC thresholds |
| `biotools translate FILE -o OUT` | DNA → protein, single or all six frames |
| `biotools fastq-stats FILE` | read length and quality distribution |

```bash
biotools validate reads.fastq.gz --format fastq --json
biotools filter genome.fa -o filtered.fa --min-length 1000 --min-gc 40 --max-gc 60
biotools translate cds.fa -o proteins.fa --all-frames --to-stop
```

## Library use

```python
from biotools import read_fasta, gc_content, n50, validate_fastq

lengths = []
for rec in read_fasta("genome.fa.gz"):      # generator: one record in RAM
    lengths.append(len(rec))
    print(rec.id, f"{gc_content(rec.sequence):.1f}%")

print("N50:", n50(lengths))

report = validate_fastq("reads.fastq.gz", strict=False)
if not report["valid"]:
    for err in report["errors"][:5]:
        print(err)     # 'record 8,213 (SRR1039508.8213): seq/qual length mismatch'
```

## Correctness details

These are the places where a naive implementation quietly gives wrong answers:

**Quality averaging.** Phred scores are logarithmic, so the arithmetic mean of
the integers overstates read quality. `mean_quality()` converts to error
probability, averages that, and converts back. A read that is half Q40 and half
Q0 averages to Q3 — not Q20.

**GC denominators.** Ns are excluded from the denominator. Counting them drags
the GC of an N-padded scaffold toward zero and makes assemblies look
artificially AT-rich.

**Streaming, not slurping.** `read_fasta` and `read_fastq` are generators. Peak
memory is one record, so a 40 GB genome costs the same RAM as a 4 KB one.

**Gzip by magic bytes**, not by file extension — `.fa` files that are actually
gzipped are common in the wild.

**Errors name the record.** A validation failure reports the record number and
ID (`record 8,213 (SRR1039508.8213): seq/qual length mismatch (63 vs 61)`), so
you can find the problem in a 10-million-record file.

**Ambiguity codes translate to X** rather than raising. Real sequence data has
Ns in it; a parser that crashes on them is not usable on real data.

## Tests

```bash
pip install -e ".[dev]"
pytest -v
```

29 tests covering parsing (wrapped lines, gzip, soft-masking, malformed input),
statistics (GC with Ns, N50 edge cases, reverse-complement involution),
translation (all frames, ambiguity, stop handling), and validation (invalid
alphabets, duplicate IDs, truncated records, Phred range).

## Requirements

Python ≥ 3.10. No runtime dependencies.

## License

MIT
