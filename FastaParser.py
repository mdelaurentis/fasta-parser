#!/usr/bin/env python

import unittest
from itertools import dropwhile, takewhile

class Entry:
    def __init__(self, gi, accession, description, sequence=""):
        self.gi = gi
        self.accession = accession
        self.description = description
        self.sequence = sequence

    def __str__(self):
        header_parts = ("gi", str(self.gi),
                        "ref", self.accession,
                        self.description)
        return ">" + "|".join(header_parts) + "\n" + self.sequence
            
class FastaParser:

    def __init__(self, filename=None):
        """Create a parser that operates on the given file."""
        self.filename = filename

    def _parse_header_line(self, text):
        """Attempt to parse the given text as a FASTA header line and
        create an Entry out of it. If the line does not appear to be a
        FASTA header line, returns None, otherwise returns the new
        Entry."""

        if len(text) == 0 or text[0] != ">":
            return None
        parts = text[1:].split("|")
        (_, gi, _, accession, description) = parts
        return Entry(long(gi), accession, description)

    def entries(self, start=None, stop=None):
        """Returns a generator of the sequence of entries in the file."""
        entry = None
        with open(self.filename) as infile:
            for line in infile:
                line = line.rstrip("\n")

                # Attempt to parse the line as a header; if we get
                # None back that means it's not a header line.
                thisentry = self._parse_header_line(line)

                # If this line is an entry, we need to yield the
                # previous entry (unless this is the first entry) and
                # then set this line as our current entry.
                if thisentry is not None:
                    if entry is not None:
                        yield entry
                    entry = thisentry

                # Otherwise if it's not an entry assume it's a
                # sequence line, and append it to the current entry's
                # sequence
                else:
                    entry.sequence += line

        yield entry
        infile.close()
        return

    def entry(self, i):
        """Returns the ith entry (starting numbering at 1), or None if
        i is outside the acceptable range."""
        
        entries = list(self.entry_range(i, i+1))
        return None if len(entries) == 0 else entries[0]

    def first(self):
        """Returns the first entry in the file or None if there are no
        entries."""
        return self.entry(1)

    def last(self):
        """Returns the last entry in the file or None if there are no
        entries."""

        # Don't just get the entry where the index is the same as the
        # count of entries; that would require iterating through the
        # entries twice.
        result = None
        for e in self.entries():
            result = e
        return result

    def entry_range(self, start, stop):
        """Returns the subsequence of entries in the range from start
        to stop (both inclusive), starting at 1 for the first
        entry. If start or stop is outside the actual range of
        indices, simply ignores the parts of the range that are
        invalid."""

        i = 1
        for e in self.entries():
            if i > stop:
                return
            if i >= start:
                yield e
            i += 1

    def count(self):
        """Returns the number of entries in the file."""
        return sum(1 for e in self.entries())

    def __len__(self):
        """Returns the number of entries in the file."""
        return self.count()

class TestFastaParser(unittest.TestCase):

    test_file = "first100.fna"
    sample_gis = [197313646L, 197313649L, 197313647L, 215983060L, 32452934L]

    valid_header=">gi|355477125|ref|NW_001493874.3| Bos taurus breed Hereford chromosome 1 genomic scaffold, alternate assembly Btau_4.6.1 Chr1.scaffold45"

    first_entry = Entry(197313646L, "NR_001588.2", " Homo sapiens Shwachman-Bodian-Diamond syndrome pseudogene 1 (SBDSP1), transcript variant 3, non-coding RNA", "CCTTTTTGGGCGTGGAAAGATGGCGGTAAAAGCCACAATGCGCAGGCGTCATCGCTCACTTCTCCCCTCCCGGCTTCTGCTCCACCTGACGCCTGCGCAGTAAGTAAGCCTGCCAGACACGCTGTGGCGGCTGCCTGAAGCTAGTGAGTCGCGGCGCCGCGCACTTGTGGTTGGGTCAGTGCCGCGCGCCGCTCGGTCGTTACCGCGAGGCGCTGGTGGCCTTCAGGCTGGACGGCGCGGGTCAGCCCTGGTTTGCCGGCTTCTGGGTCTTTGAACAGCCGCGATGTCGATCTTCACCCCCACCAACCAGATCCGCCTAACCAATGTGGCCGTGGTACGGATGAAGCGCGCCAGGAAGCGCTTCGAAATCGCCTGCTACAGAAACAAGGTCGTCGGCTGGCGGAGCGGCTTGGAAAAAGACCTTGATGAAGTTCTGCAGACCCACTCAGTGTTTGTAAATGTTTCCTAAGGTCAGGTTGCCAAGAAGGAAGATCTCATCAGTGCGTTTGGAACAGATGACCAAACTGAAATCTATTTTGACTAAAGGAGAAGTTCAAGTATCAGATAAAGACACACACAACTGGAGCAGATGTTTAGGGACATTGCAATTATTGTGGCAGACAAATGTGTGACTCCTGAAACAAAGAGACCATACACCGTGATCCTTATTGAGAGAGCCATGAAGGACATCCACTATTTGGTGAAAACCAACAGGAGTACAAAACAGCAGGCTTTGGAAGTGATAAAGCAGTTAAAAGAGAAAATGAAGATAGAACGTGCTCACATGAGGCTTCAGTTCATCCTTCCAGTGAATGAAGGCAAGAAGCTGAAAGAAAAGCTCAAGCCACTGATCAAGGTCATAGAAAGTAAAGATTATGGCCAACAGTTAGAAATCGTAAGAGTCAAATATTTTCTTTGCTTCATGTTACCTAAATATTGTATTCTCTAGTAATAAATTTGTAGCAAACATTCAAAAAAAAAAAAAAAAAAAA")

    def parser(self):
        return FastaParser(self.test_file)

    def test_parse_header_line(self):
        fp = self.parser()
        entry = fp._parse_header_line(self.valid_header)
        self.assertEquals(355477125, entry.gi)
        self.assertEquals("NW_001493874.3", entry.accession)
        self.assertEquals(" Bos taurus breed Hereford chromosome 1 genomic scaffold, alternate assembly Btau_4.6.1 Chr1.scaffold45", entry.description)

    def test_parse_header_line_returns_none(self):
        fp = self.parser()
        self.assertIsNone(fp._parse_header_line("foobar"))

    def test_len(self):
        fp = self.parser()
        self.assertEquals(5, len(fp))
        self.assertEquals(5, len(fp))

    def test_entries(self):
        fp = self.parser()
        entries = list(fp.entries())
        self.assertEquals(self.first_entry.gi, entries[0].gi)
        self.assertEquals(self.first_entry.accession, entries[0].accession)
        self.assertEquals(self.first_entry.description, entries[0].description)
        self.assertEquals(self.first_entry.sequence, entries[0].sequence)

    def test_entry_range(self):
        entries = self.parser().entry_range(2, 4)
        expected = self.sample_gis[1:4]
        got = [e.gi for e in entries]
        self.assertEquals(expected, got)

    def test_entry_range_ignores_invalid_start(self):
        expected = self.sample_gis[0:3]
        got = [e.gi for e in self.parser().entry_range(-10, 3)]
        self.assertEquals(expected, got)

    def test_entry_range_ignores_invalid_stop(self):
        expected = self.sample_gis[2:]
        got = [e.gi for e in self.parser().entry_range(3, 10)]
        self.assertEquals(expected, got)

    def test_entry(self):
        expected = self.sample_gis[2]
        got = self.parser().entry(3).gi
        self.assertEquals(expected, got)

    def test_invalid_entry(self):
        self.assertIsNone(self.parser().entry(17))

    def test_first(self):
        expected = self.sample_gis[0]
        got = self.parser().first().gi
        self.assertEquals(expected, got)

    def test_last(self):
        expected = self.sample_gis[4]
        got = self.parser().last().gi
        self.assertEquals(expected, got)
        
if __name__ == '__main__':
    unittest.main()
