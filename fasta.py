#!/usr/bin/env python

import unittest
import pickle
from itertools import islice
import os
import sys
import logging

# Set up logging. TODO: I'm not to familiar with Python logging; there
# may be a more idiomatic way to do it, in a config file or something.
logger = logging.getLogger("fasta")
hdlr = logging.FileHandler('fasta.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.WARNING)

class OutOfBoundsException(Exception):
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
        return ">" + "|".join(header_parts) + "\n" + self._format_sequence()

    def _format_sequence(self, max_len=80):
        subsequences = []
        for i in range(0, len(self.sequence), max_len):
            subsequences.append(self.sequence[i:i+max_len])
        return "\n".join(subsequences)
            
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
            logger.warn(
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
            logger.warn("Warning: couldn't clear index")
        
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

    def _position(self, entry_num):
        if entry_num < 1:
            raise OutOfBoundsException("Index must be greater than 0: " + str(entry_num))
        if self.index is None:
            return None
        if entry_num > len(self.index):
            raise OutOfBoundsException("No entry with number " + str(entry_num))
        return self.index[entry_num - 1]

    def entries(self, offset=None):
        """Returns a generator of the sequence of entries in the
        file. If offset is provided, it must be the byte offset into
        the file where I should start parsing."""

        entry = None
        
        with open(self.filename) as infile:
            if offset is not None:
                infile.seek(offset)
            while True:
                pos = infile.tell()
                line = infile.readline()
                                
                if line == "":
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

        if entry is not None:
            yield entry
        infile.close()
        return

    def entry(self, i):
        """Returns the ith entry (starting numbering at 1), or None if
        i is outside the acceptable range."""

        if i < 1:
            raise OutOfBoundsException("Index must be greater than 0: " + str(i))
        entries = None

        if self.index is not None:
            if i > len(self):
                raise OutOfBoundsException("No entry with number " + str(i))
            offset = self.index[i - 1]
            entries = islice(self.entries(offset), 0, 1)
        else:
            entries = islice(self.entries(), i - 1, i)
        entries = list(entries)
        if len(entries) == 0:
            raise OutOfBoundsException("No entry with number " + str(i))
        return list(entries)[0]

    def first(self):
        """Returns the first entry in the file or None if there are no
        entries."""
        return self.entry(1)

    def last(self):
        """Returns the last entry in the file or None if there are no
        entries."""
        return self.entry(len(self))

    def count(self):
        """Returns the number of entries in the file."""
        if self.index is not None:
            return len(self.index)
        return sum(1 for e in self.entries())

    def __len__(self):
        """Returns the number of entries in the file."""
        return self.count()
