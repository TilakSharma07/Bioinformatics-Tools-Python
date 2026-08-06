"""biotools - small, dependency-light utilities for FASTA/FASTQ handling."""

from biotools.fasta import FastaRecord, read_fasta, write_fasta
from biotools.fastq import FastqRecord, read_fastq, mean_quality
from biotools.seqstats import gc_content, nucleotide_counts, n50, translate, reverse_complement
from biotools.validate import ValidationError, validate_fasta, validate_fastq

__version__ = "0.2.0"
__all__ = [
    "FastaRecord", "read_fasta", "write_fasta",
    "FastqRecord", "read_fastq", "mean_quality",
    "gc_content", "nucleotide_counts", "n50", "translate", "reverse_complement",
    "ValidationError", "validate_fasta", "validate_fastq",
]
