
Contents
========

fasta.py is a module with two classes:

* Entry - Represents a single entry in a fasta file.
* FastaParser - A parser that is initialized with the name of a fasta
  file and knows how to parse it.

fasta_test.py contains the unit tests for fasta.

parse_fasta.py is a command line program that you can use to
experiment with the parser.

Strategy
========

Python has a "generator", which provides lazy iteration over
values. I'm using generators to lazily read in the file so that we
never store the entire contents of the file in memory.

As FastaParser reads the entries in the file, it keeps track of the
start position of each entry. Once it reaches the end of the file, it
stores this data as an index in a *.idx file. Before this index is
created, any requests to access an entry in the file will cause
FastaParser to read and parse entries starting at the beginning until
it gets to the requested entry. After the index is created, when
handling a request for a particular entry or range of entries, it will
seek to the appropriate point in the file and start parsing from
there. This makes random access much faster.

Command-Line Interface
======================

parse_fasta is a command line interface that you can use to play with
the parser. If you run ./parse_fasta -h, it will print out some help
information.