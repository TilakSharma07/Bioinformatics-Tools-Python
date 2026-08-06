"""Sequence statistics: GC, composition, N50, translation."""
from __future__ import annotations

from collections import Counter
from typing import Iterable

# Standard genetic code (NCBI translation table 1).
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTUNRYKMSWBDHVacgtunrykmswbdhv",
                            "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb")


def gc_content(seq: str) -> float:
    """GC as a percentage, ignoring Ns and gaps.

    Ns are excluded from the denominator - counting them would drag the GC
    of an N-padded scaffold artificially toward zero.
    """
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    if acgt == 0:
        return 0.0
    return 100.0 * (s.count("G") + s.count("C")) / acgt


def nucleotide_counts(seq: str) -> dict[str, int]:
    """Counts for A/C/G/T/N plus 'other' for anything else present."""
    c = Counter(seq.upper())
    out = {b: c.get(b, 0) for b in "ACGTN"}
    out["other"] = sum(v for k, v in c.items() if k not in "ACGTN")
    return out


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def n50(lengths: Iterable[int]) -> int:
    """Assembly N50: the length L where contigs >= L cover half the assembly."""
    lens = sorted((x for x in lengths if x > 0), reverse=True)
    if not lens:
        return 0
    half, running = sum(lens) / 2, 0
    for length in lens:
        running += length
        if running >= half:
            return length
    return lens[-1]


def translate(seq: str, frame: int = 0, to_stop: bool = False) -> str:
    """Translate DNA to protein in the given frame (0, 1 or 2).

    Codons containing N or any non-ACGT character become 'X' rather than
    raising - real sequence data has ambiguity codes in it.
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0, 1 or 2 (got {frame})")
    s = seq.upper().replace("U", "T")[frame:]
    aa = []
    for i in range(0, len(s) - len(s) % 3, 3):
        residue = CODON_TABLE.get(s[i : i + 3], "X")
        if residue == "*" and to_stop:
            break
        aa.append(residue)
    return "".join(aa)
