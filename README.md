# Razor
Razor is a fast and simple adapter/quality trimming tool for fastq sequences.
It uses the ERNE-FILTER[1] algorthm for quality trimming, which operates on both ends of a read. 
The 3' adapter-trimming algorithm works by moving from left to right in a read until a match is found, and then removing all bases after that point.
Currently, Razor only works on unix-like systems, and requires the zlib compression library (www.zlib.net).

Razor supports reading and writing both uncompressed and gzipped fastq files; the format is determined by the filename suffix.
## Installation
    git clone https://raw.githubusercontent.com/abardill/razor
    cd razor
    gcc -o razor razor.c -Wall -lz -O2
## Usage
Perform quality trimming with a threshold of 20, removing reads whose mean quality post-trimming is less than 15:

	$ razor -q 20 -m 15 -o out.fq.gz in.fq.gz
	82408952 reads were read
     1871332 reads underwent quality trimming
      715766 of these were filtered
    10734758 bases were quality trimmed

Trim a 3' adapter, requiring its first 16 bases to be matched:

	$ razor -a CTGTCTCTTATACACATCT --min-adapter-match 16 -o out.fq in.fq
	13584188 reads were read
   	   22675 reads underwent adapter trimming
          35 of these were filtered
 	 1099363 adapter bases were trimmed

Read from standard input (denoted by '-') and trim adapters contained in a fasta file, discarding reads shorter than 30 bases. Write to standard output:

	$ cat in.fq | razor -f adapters.fa -l 30 - > out.fq

Combine adapter and quality trimming (in that order), but retain reads which do not pass the --min-length and/or --min-mean-quality filters, replacing their sequence and quality strings with empty strings:

	$ razor -f adapters.fa -q 30 --keep-empty -o out.fq.gz in.fq.gz
	
This is useful if you want to maintain pairing when trimming paired-end reads.

To print the full help information use

	$ razor --help
	Usage: razor [options] [in.fq]
	Options: 
	  -o, --output=FILE		The file to output processed reads to. If unspecified, reads are written to stdout.
	  -q, --quality-threshold=INT	Trim bases from 5' and 3' ends using the ERNE-FILTER algorithm.
	  -l, --min-length=INT		After processing, discard sequences below a minimum length [20].
	  -m, --min-mean-quality=INT	After quality trimming, discard reads below a minimum mean quality.
	  -a, --adapter=STR		Adapter to be trimmed from the 3' ends of the reads.
	  -f, --adapter-file=FILE	Fasta file containing adapter sequences.
	  -M, --min-adapter-match=INT	The minimum number of leading bases a subsequence must share with an adapter
					for a match to occur [12].
	  -k, --keep-empty		Sequences filtered according to the min-length and min-mean-quality options
					have their sequence and quality strings truncated to empty and are retained.
	  -p, --phred=[33|64]		The phred encoding of the quality scores in the input fastq [33].
	      --version			Display version information and exit.
	      --help			Display this help text and exit.

 
## References
[1] Fabbro, Cristian Del, Simone Scalabrin, Michele Morgante, and Federico M. Giorgi. “An Extensive Evaluation of Read Trimming Effects on Illumina NGS Data Analysis.” PLOS ONE 8, no. 12 (December 23, 2013): e85024. https://doi.org/10.1371/journal.pone.0085024.

