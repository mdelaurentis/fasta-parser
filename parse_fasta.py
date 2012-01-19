#!/usr/bin/env python

import unittest

class Entry:
    def __init__(self, text):
        self.text = text

    def gi():
        return "Foo"

    def accession(self):
        return "Bar"

    def description(self):
        return "Bar"

    def sequence(self):
        return "Foo"
        
class FastaParser:

    def __init__(self, filename=None):
        self.filename = filename

    def _is_header_line(self, text):
        return len(text) > 0 and text[0] == ">"

    def _parse_header_line(self, text):
        if len(text) == 0 or text[0] != ">":
            return None
        parts = text[1:].split("|")
        (_, gi, _, accession, description) = parts
        return (gi, accession, description)

class TestFastaParser(unittest.TestCase):

    fp = FastaParser()

    valid_header=">gi|355477125|ref|NW_001493874.3| Bos taurus breed Hereford chromosome 1 genomic scaffold, alternate assembly Btau_4.6.1 Chr1.scaffold45"

    def test_is_header_line(self):

        self.assertTrue(self.fp._is_header_line(">foobar"))
        self.assertFalse(self.fp._is_header_line("foobar"))


    def test_parse_header_line(self):
        (gi, accession, description) = \
            self.fp._parse_header_line(self.valid_header)
        self.assertEquals("355477125", gi)
        self.assertEquals("NW_001493874.3", accession)
        self.assertEquals(" Bos taurus breed Hereford chromosome 1 genomic scaffold, alternate assembly Btau_4.6.1 Chr1.scaffold45", description)
        
if __name__ == '__main__':
    unittest.main()
