#!/usr/bin/env python

from fasta import FastaParser, Entry
import unittest

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
