"""Streaming FASTQ reader and Phred quality helpers."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from biotools.fasta import _open_text

PHRED_OFFSET = 33  # Sanger / Illumina 1.8+


@dataclass(slots=True)
class FastqRecord:
    id: str
    sequence: str
    quality: str

    def __len__(self) -> int:
        return len(self.sequence)

    def phred_scores(self) -> list[int]:
        return [ord(c) - PHRED_OFFSET for c in self.quality]


def read_fastq(path: str | Path) -> Iterator[FastqRecord]:
    """Yield FastqRecord objects, four lines at a time.

    Raises ValueError on a truncated final record or a seq/qual length
    mismatch - both are silent-corruption bugs if left unchecked.
    """
    with _open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            header = header.rstrip("\n\r")
            seq = fh.readline().rstrip("\n\r")
            plus = fh.readline().rstrip("\n\r")
            qual = fh.readline().rstrip("\n\r")

            if not qual:
                raise ValueError(f"{path}: truncated record at '{header[:40]}'")
            if not header.startswith("@"):
                raise ValueError(f"{path}: expected '@' header, got '{header[:40]}'")
            if not plus.startswith("+"):
                raise ValueError(f"{path}: expected '+' separator at '{header[:40]}'")
            if len(seq) != len(qual):
                raise ValueError(
                    f"{path}: seq/qual length mismatch at '{header[:40]}' "
                    f"({len(seq)} vs {len(qual)})"
                )
            yield FastqRecord(id=header[1:].split()[0], sequence=seq, quality=qual)


def mean_quality(record: FastqRecord) -> float:
    """Error-probability mean, not the arithmetic mean of Phred scores.

    Phred is a log scale, so averaging the integers overstates quality. The
    correct average converts to error probability, means those, converts back.
    """
    import math

    scores = record.phred_scores()
    if not scores:
        return 0.0
    mean_p = sum(10 ** (-q / 10) for q in scores) / len(scores)
    if mean_p <= 0:
        return 60.0
    return -10 * math.log10(mean_p)
