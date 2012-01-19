
Contents
========

fasta.py is a module with two classes:

* Entry - Represents a single entry in a fasta file.
* FastaParser - A parser that is initialized with the name of a fasta
  file and knows how to parse it.

fasta_test.py contains the unit tests for fasta.

parse_fasta is a command line program that you can use to
experiment with the parser.


Strategy
========

Lazy Iteration
--------------

Python has a "generator", which provides lazy iteration over
values. I'm using generators to lazily read in the file so that we
never store the entire contents of the file in memory.

Indexing
--------

FastaParser has two modes of operation: index and non indexed.

In non-indexed mode, all accesses require the parser to read sequentially through the file.

If you ask it to index a file (by calling save_index()) it will
traverse the file and record the start position of each entry. It then
saves that index in an *.idx file. For all subsequent random accesses
it will seek to the appropriate position in the file and start reading
entries from there.

Command-Line Interface
======================

parse_fasta is a command line interface that you can use to play with
the parser. If you run ./parse_fasta -h, it will print out some help
information.

Index
-----

`parse_fasta index <infile>` will force the index for infile to be created.

Count
-----

`parse_fasta count <infile>` while print the number of entries in infile. It will also create the index if it doesn't already exist. You can use this to see the speedup that the index provides:

  > time ./parse_fasta count vertebrate_mammalian.2.rna.fna 
  Couldn't load index at vertebrate_mammalian.2.rna.fna.idx
  I have 30935 entries

  real	0m2.374s
  user	0m2.315s
  sys	0m0.055s
  mike@trc459 ~/src/fasta-parser
  > time ./parse_fasta count vertebrate_mammalian.2.rna.fna 
  I have 30935 entries

  real	0m0.171s
  user	0m0.151s
  sys	0m0.018s

Show
----

You can see an individual entry:

    mike@trc459 ~/src/fasta-parser
    > ./parse_fasta show --nth 1000 vertebrate_mammalian.2.rna.fna
    >gi|291384767|ref|XM_002709208.1| PREDICTED: Oryctolagus cuniculus hypothetical protein LOC100345993 (LOC100345993), mRNA
    ATGGCGGGAGGGAAGGCCACCTTGGAATTTCTTCCCGAGTCACCCCCGGGGACGACGAGAGGGGAAAG
    -- snip --
    mike@trc459 ~/src/fasta-parser

  Or a range of entries

    ./parse_fasta show --start 30 --stop 40 vertebrate_mammalian.2.rna.fna

  Or all the entries:

    ./parse_fasta show vertebrate_mammalian.2.rna.fna
