
Contents
========

fasta/parser.py is a module with two classes:

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

In non-indexed mode, all accesses require the parser to read
sequentially through the file.

You can index the file by calling save_index() on the parser (or
`running parse_fasta index` at the command line). To index the file,
the parser simply traverses through the file and records the byte
offset of the start of each entry. It stores on disk an array mapping
the entry's sequence number in the file to its byte offset. Once a
file is indexed, the parser handles requests for a particular entry by
seeking to that point in the file rather than reading in the whole file.        

Command-Line Interface
======================

parse_fasta is a command line interface that you can use to play with
the parser. If you run ./parse_fasta -h, it will print out some help
information.

Index
-----

`parse_fasta index <infile>` will force the index for infile to be created.

    > ./parse_fasta index vertebrate_mammalian.2.rna.fna 

Count
-----

`parse_fasta count <infile>` while print the number of entries in infile. It will also create the index if it doesn't already exist. You can use this to see the speedup that the index provides:

    > time ./parse_fasta count vertebrate_mammalian.2.rna.fna 
    I have 30935 entries
    
    real    0m2.214s
    user    0m2.158s
    sys     0m0.053s
    mike@trc459 ~/src/fasta-parser
    > ./parse_fasta index vertebrate_mammalian.2.rna.fna 
    mike@trc459 ~/src/fasta-parser
    > time ./parse_fasta count vertebrate_mammalian.2.rna.fna 
    I have 30935 entries
    
    real    0m0.173s
    user    0m0.151s
    sys     0m0.017s

Show
----

You can see an individual entry:

    > ./parse_fasta show --nth 1000 vertebrate_mammalian.2.rna.fna 
    >gi|291384767|ref|XM_002709208.1| PREDICTED: Oryctolagus cuniculus hypothetical protein LOC100345993 (LOC100345993), mRNA
    ATGGCGGGAGGGAAGGCCACCTTGGAATTTCTTCCCGAGTCACCCCCGG...
    -- snip --

Or all entries:

    > ./parse_fasta show vertebrate_mammalian.2.rna.fna 

Possible Enhancements
=====================

Validation
----------

I am currently not doing much to validate the file. For header lines,
I'm just stripping off the initial '>' character, splitting on '|' and
expecting the fields to be in the proper order. It might be good to
validate the gi number and accession number. For sequence lines, I'm
accepting any text. It might be good to ensure that the sequence
consists only of valid characters.

Compression
-----------

It would be interesting (and not hard) to allow the parser to read directly from a compressed file.