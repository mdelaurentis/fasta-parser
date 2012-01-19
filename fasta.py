#!/usr/bin/env python

import unittest
import pickle
from itertools import dropwhile, takewhile
import os
import sys

class IndexException(Exception):
    pass

class Entry:
    """Represents an entry in a FASTA file, including the header
    information and the sequence."""
    
    def __init__(self, gi, accession, description, sequence="", pos=None):
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
    """Parses fasta files. Once we reach the end of a file foo, we
    index the entries in the file and save the index to foo.idx so
    that we can quickly access random entries in the file later."""    

    def __init__(self, filename=None, index_filename=None):
        """Create a parser that operates on the given file."""
        self.filename = filename
        self.index_filename = index_filename
        if self.index_filename is None:
            self.index_filename = filename + ".idx"
        self.index = None
        try:
            with open(self.index_filename) as infile:
                self.index = pickle.load(infile)
            infile.close()
        except IOError:
            sys.stderr.write(
                "Couldn't load index at %s, proceeding without index for %s\n"
                % (self.index_filename, self.filename))

    def save_index(self, force=False):
        """Saves the index to the filename specified in my
        index_filename property."""
        self.index = [e.pos for e in self.entries()]
        with open(self.index_filename, "w") as outfile:
            pickle.dump(self.index, outfile)
        outfile.close()

    def clear_index(self):
        """Removes my index file."""
        try:
            os.remove(self.index_filename)
        except:
            raise IndexException("Couldn't clear index")
        
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

    def entries(self):
        """Returns a generator of the sequence of all
        the entries in the file."""
        return self._parse()

    def _parse(self, start=None, stop=None):
        """Returns a generator of the sequence of entries in the
        file. If start is provided, it must be the byte offset into
        the file where I should start parsing. If stop is provided, it
        must be the byte offset where I should stop parsing. If I
        reach the end of the file, I will set my index property so
        that the ith entry in the index points to the position of the
        ith entry in the file."""
        
        entry = None
        counter = 0

        with open(self.filename) as infile:
            if start is not None:
                infile.seek(start)
            while True:
                pos = infile.tell()

                line = infile.readline()
                if line == "" or (stop is not None and pos >= stop):
                    break
                line = line.rstrip("\n")

                # Attempt to parse the line as a header; if we get
                # None back that means it's not a header line.
                thisentry = self._parse_header_line(line)

                # If this line is an entry, we need to yield the
                # previous entry (unless this is the first entry) and
                # then set this line as our current entry.
                if thisentry is not None:
                    thisentry.pos = pos
                    if entry is not None:
                        yield entry
                    entry = thisentry

                # Otherwise if it's not an entry assume it's a
                # sequence line, and append it to the current entry's
                # sequence
                else:
                    if entry is None:
                        raise Exception(
                            "%s does not appear to be a valid FASTA file" % self.filename)
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
        return self.entry(len(self))


    def _entry_range_from_index(self, start, stop):
        """Parses the entries in the given range from the file using
        the index (both start and stop are inclusive). This should be
        much faster than _entry_range_from_stream for ranges that skip
        over lots of entries in the top of the file."""

        # If start is after the acceptable range or stop is before it,
        # we shouldn't do anything at all.
        if start > self.count() or stop < 1:
            return ()

        startpos = None
        stoppos = None
        # Set startpos and stoppos if they are within the range.
        if start > 0:
            startpos = self.index[start-1] 
        if stop <= self.count():
            stoppos = self.index[stop]
        
        return self._parse(startpos, stoppos)

    def _entry_range_from_stream(self, start, stop):
        """Returns the entries in the given range (inclusive) from the
        file without using an index. This will only be called on files
        that are not yet indexed."""
        
        i = 1
        for e in self.entries():
            if i > stop:
                return
            if i >= start:
                yield e
            i += 1
            
    def entry_range(self, start, stop):
        """Returns the subsequence of entries in the range from start
        to stop (both inclusive), starting at 1 for the first
        entry. If start or stop is outside the actual range of
        indices, simply ignores the parts of the range that are
        invalid."""

        if self.index is None:
            return self._entry_range_from_stream(start, stop)
        return self._entry_range_from_index(start, stop)

    def count(self):
        """Returns the number of entries in the file."""
        if self.index is not None:
            return len(self.index)
        return sum(1 for e in self.entries())

    def __len__(self):
        """Returns the number of entries in the file."""
        return self.count()

################################################################################
# Tests


class TestFastaParser(unittest.TestCase):

    test_file = "test.fna"
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
        e = self.parser().entry(17)
        self.assertIsNone(self.parser().entry(17))

    def test_first(self):
        expected = self.sample_gis[0]
        got = self.parser().first().gi
        self.assertEquals(expected, got)

    def test_last(self):
        expected = self.sample_gis[4]
        got = self.parser().last().gi
        self.assertEquals(expected, got)

    def test_index(self):
        self.parser().save_index()
        
if __name__ == '__main__':
    unittest.main()
