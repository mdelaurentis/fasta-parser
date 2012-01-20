#!/usr/bin/env python

from fasta.parser import FastaParser, Entry, OutOfBoundsException
import unittest

class TestFastaParser(unittest.TestCase):

    test_file = "fasta/test.fna"
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
        fp.clear_index()
        self.assertEquals(5, len(fp))

    def test_indexed_len(self):
        fp = self.parser()
        fp.save_index()
        self.assertEquals(5, len(fp))

    def test_entries(self):
        fp = self.parser()
        entries = list(fp.entries())
        self.assertEquals(self.first_entry.gi, entries[0].gi)
        self.assertEquals(self.first_entry.accession, entries[0].accession)
        self.assertEquals(self.first_entry.description, entries[0].description)
        self.assertEquals(self.first_entry.sequence, entries[0].sequence)

    def test_entry(self):
        expected = self.sample_gis[2]
        fp = self.parser()
        fp.clear_index()
        self.assertEquals(expected, fp.entry(3).gi)

    def test_indexed_entry(self):
        expected = self.sample_gis[2]
        fp = self.parser()
        fp.save_index()
        self.assertEquals(expected, fp.entry(3).gi)

    def test_invalid_entry(self):
        fp = self.parser()
        fp.clear_index()
        self.assertRaises(OutOfBoundsException, fp.entry, 17)
        self.assertRaises(OutOfBoundsException, fp.entry, -1)

    def test_indexed_invalid_entry(self):
        fp = self.parser()
        fp.save_index()
        self.assertRaises(OutOfBoundsException, fp.entry, 17)
        self.assertRaises(OutOfBoundsException, fp.entry, -1)        

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

    def test_format_sequence(self):
        format4 = lambda seq: Entry(0, None, None, seq)._format_sequence(4)
        self.assertEquals("abcd\nefgh\nijkl", format4("abcdefghijkl"))
        self.assertEquals("abcd\nefgh\nijk", format4("abcdefghijk"))
        self.assertEquals("abcd\nefgh\ni", format4("abcdefghi"))
        
        
if __name__ == '__main__':
    unittest.main()
