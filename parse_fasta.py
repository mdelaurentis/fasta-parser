#!/usr/bin/env python

import unittest
from FastaParser import FastaParser
import sys
import argparse

parser = argparse.ArgumentParser(description='Browse a fasta file')
parser.add_argument("--start", help="The index of the first entry to show")
parser.add_argument("--stop", help="The index of the last entry to show")
parser.add_argument("--first", help="Show just the first entry")
parser.add_argument("--last", help="Show just the last entry")
parser.add_argument("--count",
                    help="Just show the number of entries in the file",
                    action="store_const",
                    const=True)
parser.add_argument("infile")
args = parser.parse_args()

parser = FastaParser(args.infile)

if args.count:
    print "I found %d entries" % len(parser)


#for e in fp.entries():
#    print e

