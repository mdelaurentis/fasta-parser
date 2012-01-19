#!/usr/bin/env python

import unittest

class Entry:
    def __init__(self, gi, accession, description):
        self.gi = gi
        self.accession = accession
        self.description = description
        self.sequence = ""

    def __str__(self):
        header_parts = ("gi", self.gi,
                        "ref", self.accession,
                        self.description)
        return ">" + "|".join(header_parts) + "\n" + self.sequence
        
class FastaParser:

    def __init__(self, filename=None):
        self.filename = filename

    def _parse_header_line(self, text):
        if len(text) == 0 or text[0] != ">":
            return None
        parts = text[1:].split("|")
        (_, gi, _, accession, description) = parts
        return Entry(gi, accession, description)

    def entries(self):
        entry = None
        with open(self.filename) as infile:
            for line in infile:
                thisentry = self._parse_header_line(line)
                if thisentry is None:
                    entry.sequence += line
                else:
                    yield entry
                    entry = thisentry
        infile.close()

    def entry(self, n):
        return self.entries()[n+1]

    def count(self):
        count = 0
        for e in self.entries():
            count += 1
        return count

    def __len__(self):
        return self.count()
        

class TestFastaParser(unittest.TestCase):

    fp = FastaParser()

    valid_header=">gi|355477125|ref|NW_001493874.3| Bos taurus breed Hereford chromosome 1 genomic scaffold, alternate assembly Btau_4.6.1 Chr1.scaffold45"

    def test_parse_header_line(self):
        entry = self.fp._parse_header_line(self.valid_header)
        self.assertEquals("355477125", entry.gi)
        self.assertEquals("NW_001493874.3", entry.accession)
        self.assertEquals(" Bos taurus breed Hereford chromosome 1 genomic scaffold, alternate assembly Btau_4.6.1 Chr1.scaffold45", entry.description)

    def test_parse_header_line_returns_none(self):
        self.assertIsNone(self.fp._parse_header_line("foobar"))
        
if __name__ == '__main__':
    unittest.main()
