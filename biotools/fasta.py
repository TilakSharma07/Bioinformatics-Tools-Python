"""Streaming FASTA reader/writer.

`read_fasta` is a generator: it holds one record in memory at a time, so a
40 GB genome FASTA costs the same RAM as a 4 KB one. Handles gzip
transparently and tolerates the usual real-world messiness (blank lines,
CRLF line endings, wrapped or unwrapped sequence lines).
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


@dataclass(slots=True)
class FastaRecord:
    id: str
    description: str
    sequence: str

    @property
    def header(self) -> str:
        return f"{self.id} {self.description}".strip()

    def __len__(self) -> int:
        return len(self.sequence)


def _open_text(path: str | Path) -> io.TextIOBase:
    """Open plain or gzipped text, chosen by magic bytes rather than suffix."""
    path = Path(path)
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fasta(path: str | Path, upper: bool = True) -> Iterator[FastaRecord]:
    """Yield FastaRecord objects one at a time.

    Args:
        path:  .fasta / .fa / .fna, optionally gzipped.
        upper: uppercase the sequence (soft-masked genomes use lowercase for
               repeats; uppercasing makes downstream base counting simple).
    """
    header: str | None = None
    chunks: list[str] = []

    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _build(header, chunks, upper)
                header, chunks = line[1:].strip(), []
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence data before any '>' header")
                chunks.append(line)

    if header is not None:
        yield _build(header, chunks, upper)


def _build(header: str, chunks: list[str], upper: bool) -> FastaRecord:
    seq = "".join(chunks)
    if upper:
        seq = seq.upper()
    rec_id, _, desc = header.partition(" ")
    return FastaRecord(id=rec_id, description=desc.strip(), sequence=seq)


def write_fasta(records: Iterable[FastaRecord], path: str | Path, wrap: int = 60) -> int:
    """Write records with sequence lines wrapped at `wrap` columns.

    Returns the number of records written.
    """
    n = 0
    with open(path, "w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(f">{rec.header}\n")
            seq = rec.sequence
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    fh.write(seq[i : i + wrap] + "\n")
            else:
                fh.write(seq + "\n")
            n += 1
    return n
