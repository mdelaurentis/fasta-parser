#!/usr/bin/env python

import pickle
from itertools import islice
import os
import logging

# Set up logging. TODO: I'm not too familiar with Python logging; there
# may be a more idiomatic way to do it, in a config file or something.
logger = logging.getLogger("fasta")
logging.basicConfig(filename="fasta.log")

class OutOfBoundsException(Exception):
    """An exception thrown when someone attempts to access by index an
    entry that does not exist in the file."""
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
        """Format the sequence so that each line is at most max_len
        characters long."""
        subsequences = []
        for i in range(0, len(self.sequence), max_len):
            subsequences.append(self.sequence[i:i+max_len])
        return "\n".join(subsequences)
            
class FastaParser:
    """Parses fasta files."""

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

    def save_index(self):
        """Save the index to the filename specified in the
        index_filename property."""
        self.index = [e.pos for e in self.entries()]
        logger.info("Saving index for %s" % self.filename)
        with open(self.index_filename, "w") as outfile:
            pickle.dump(self.index, outfile)
        outfile.close()

    def clear_index(self):
        """Remove the index file."""
        try:
            os.remove(self.index_filename)
        except:
            logger.warn("Warning: couldn't clear index")
        
    def _parse_header_line(self, text):
        
        """Attempt to parse the given text as a FASTA header line and
        create an Entry out of it. If the line does not appear to be a
        FASTA header line, return None, otherwise return the new
        Entry."""

        if len(text) == 0 or text[0] != ">":
            return None
        parts = text[1:].split("|")
        (_, gi, _, accession, description) = parts
        return Entry(long(gi), accession, description)

    def entries(self, offset=None):
        
        """Return a generator of the sequence of entries in the
        file. If offset is provided, it must be the byte offset into
        the file where I should start parsing."""

        entry = None
        
        with open(self.filename) as infile:

            # If an offset was supplied, start parsing there.
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
                    # If there's no entry then that means we saw a
                    # sequence line before a header line, which is
                    # invalid.
                    if entry is None:
                        raise Exception(
                            "%s does not appear to be a valid FASTA file" % self.filename)
                    entry.sequence += line
                    
        infile.close()
        if entry is not None:
            yield entry
        return

    def entry(self, i):
        
        """Return the ith entry (starting numbering at 1), or raise an
        OutOfBoundsException if i is outside the acceptable range."""

        if i < 1:
            raise OutOfBoundsException(
                "Index must be greater than 0: " + str(i))
        
        entries = None

        # If I am using an index, I can check to make sure that i is
        # in the appropriate range. Then find the offset of the entry
        # and skip to that point in the file.
        if self.index is not None:
            if i > len(self):
                raise OutOfBoundsException("No entry with number " + str(i))
            offset = self.index[i - 1]
            entries = islice(self.entries(offset), 0, 1)

        # Otherwise I need to just skip over i - 1 entries
        else:
            entries = islice(self.entries(), i - 1, i)
        entries = list(entries)
        if len(entries) == 0:
            raise OutOfBoundsException("No entry with number " + str(i))
        return list(entries)[0]

    def first(self):
        """Return the first entry in the file or raise an
        OutOfBoundsException if there are no entries."""
        return self.entry(1)

    def last(self):
        """Return the last entry in the file or raise an
        OutOfBoundsException if there are no entries."""
        return self.entry(len(self))

    def count(self):
        """Return the number of entries in the file."""
        if self.index is not None:
            return len(self.index)
        return sum(1 for e in self.entries())

    def __len__(self):
        """Return the number of entries in the file."""
        return self.count()
