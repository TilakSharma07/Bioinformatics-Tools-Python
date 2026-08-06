"""Test suite. Run: pytest -v"""
from __future__ import annotations

import gzip

import pytest

from biotools.fasta import FastaRecord, read_fasta, write_fasta
from biotools.fastq import mean_quality, read_fastq
from biotools.seqstats import (gc_content, n50, nucleotide_counts,
                               reverse_complement, translate)
from biotools.validate import ValidationError, validate_fasta, validate_fastq

FASTA = ">seq1 first record\nATGCATGC\nATGC\n>seq2 second\nGGGGCCCC\n"
FASTQ = ("@read1\nATGCATGC\n+\nIIIIIIII\n"
         "@read2\nGGGGCCCC\n+\n!!!!!!!!\n")


@pytest.fixture
def fasta_file(tmp_path):
    p = tmp_path / "t.fasta"
    p.write_text(FASTA)
    return p


@pytest.fixture
def fastq_file(tmp_path):
    p = tmp_path / "t.fastq"
    p.write_text(FASTQ)
    return p


class TestFasta:
    def test_parses_all_records(self, fasta_file):
        recs = list(read_fasta(fasta_file))
        assert len(recs) == 2
        assert recs[0].id == "seq1"
        assert recs[0].description == "first record"

    def test_joins_wrapped_lines(self, fasta_file):
        assert list(read_fasta(fasta_file))[0].sequence == "ATGCATGCATGC"

    def test_reads_gzip(self, tmp_path):
        p = tmp_path / "t.fasta.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(FASTA)
        assert len(list(read_fasta(p))) == 2

    def test_lowercase_preserved_when_asked(self, tmp_path):
        p = tmp_path / "soft.fasta"
        p.write_text(">s\natgcATGC\n")
        assert list(read_fasta(p, upper=False))[0].sequence == "atgcATGC"
        assert list(read_fasta(p, upper=True))[0].sequence == "ATGCATGC"

    def test_rejects_data_before_header(self, tmp_path):
        p = tmp_path / "bad.fasta"
        p.write_text("ATGC\n>s\nATGC\n")
        with pytest.raises(ValueError, match="before any"):
            list(read_fasta(p))

    def test_roundtrip(self, tmp_path, fasta_file):
        out = tmp_path / "out.fasta"
        n = write_fasta(read_fasta(fasta_file), out, wrap=4)
        assert n == 2
        assert [r.sequence for r in read_fasta(out)] == ["ATGCATGCATGC", "GGGGCCCC"]


class TestFastq:
    def test_parses_records(self, fastq_file):
        recs = list(read_fastq(fastq_file))
        assert len(recs) == 2
        assert recs[0].id == "read1"
        assert recs[0].phred_scores() == [40] * 8

    def test_rejects_length_mismatch(self, tmp_path):
        p = tmp_path / "bad.fastq"
        p.write_text("@r\nATGC\n+\nII\n")
        with pytest.raises(ValueError, match="length mismatch"):
            list(read_fastq(p))

    def test_rejects_truncated(self, tmp_path):
        p = tmp_path / "trunc.fastq"
        p.write_text("@r\nATGC\n+\n")
        with pytest.raises(ValueError, match="truncated"):
            list(read_fastq(p))

    def test_mean_quality_is_error_weighted(self, fastq_file):
        """A Q40/Q0 mix must average well below 20 - error probability, not
        arithmetic mean of Phred scores."""
        recs = list(read_fastq(fastq_file))
        assert mean_quality(recs[0]) == pytest.approx(40, abs=0.01)
        assert mean_quality(recs[1]) == pytest.approx(0, abs=0.01)


class TestSeqStats:
    @pytest.mark.parametrize("seq,expected", [
        ("GGCC", 100.0), ("ATAT", 0.0), ("ATGC", 50.0), ("", 0.0),
        ("ATGCNNNN", 50.0),   # Ns excluded from the denominator
    ])
    def test_gc_content(self, seq, expected):
        assert gc_content(seq) == pytest.approx(expected)

    def test_nucleotide_counts(self):
        c = nucleotide_counts("ATGCNX")
        assert (c["A"], c["N"], c["other"]) == (1, 1, 1)

    def test_reverse_complement(self):
        assert reverse_complement("ATGC") == "GCAT"
        assert reverse_complement(reverse_complement("ATGCGTA")) == "ATGCGTA"

    @pytest.mark.parametrize("lengths,expected", [
        ([100, 200, 300, 400, 500], 400),
        ([], 0),
        ([50], 50),
    ])
    def test_n50(self, lengths, expected):
        assert n50(lengths) == expected

    def test_translate_start_and_stop(self):
        assert translate("ATGGCCTAA") == "MA*"
        assert translate("ATGGCCTAA", to_stop=True) == "MA"

    def test_translate_frames_differ(self):
        seq = "AATGGCCTAA"
        assert translate(seq, frame=1) == translate("ATGGCCTAA")

    def test_ambiguous_codon_becomes_X(self):
        assert translate("ATGNNNGCC") == "MXA"

    def test_rejects_bad_frame(self):
        with pytest.raises(ValueError, match="frame must be"):
            translate("ATGC", frame=3)


class TestValidate:
    def test_accepts_clean_fasta(self, fasta_file):
        assert validate_fasta(fasta_file)["valid"] is True

    def test_flags_invalid_characters(self, tmp_path):
        p = tmp_path / "bad.fasta"
        p.write_text(">s\nATGCZZZ\n")
        with pytest.raises(ValidationError):
            validate_fasta(p)
        assert validate_fasta(p, strict=False)["valid"] is False

    def test_flags_duplicate_ids(self, tmp_path):
        p = tmp_path / "dup.fasta"
        p.write_text(">s\nATGC\n>s\nGGCC\n")
        report = validate_fasta(p, strict=False)
        assert any("duplicate" in e for e in report["errors"])

    def test_protein_alphabet(self, tmp_path):
        p = tmp_path / "prot.fasta"
        p.write_text(">p\nMAKVLW*\n")
        assert validate_fasta(p, alphabet="protein")["valid"] is True

    def test_fastq_report_has_phred_range(self, fastq_file):
        r = validate_fastq(fastq_file)
        assert (r["phred_min"], r["phred_max"], r["records"]) == (0, 40, 2)
