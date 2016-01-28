

<!--
 ! This manual is written in "markdown" format and thus contains some
 ! distracting formatting clutter.  See 'MANUAL' for an easier-to-read version
 ! of this text document, or see the HTML manual online.
 ! -->

Introduction
============

What is Centrifuge?
-----------------

[Centrifuge] is a novel microbial classification engine that enables
rapid, accurate and sensitive labeling of reads and quantification of
species on desktop computers.  The system uses a novel indexing scheme
based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini
(FM) index, optimized specifically for the metagenomic classification
problem. Centrifuge requires a relatively small index (4.7 GB for all
complete bacterial and viral genomes plus the human genome) and
classifies sequences at very high speed, allowing it to process the
millions of reads from a typical high-throughput DNA sequencing run
within a few minutes.  Together these advances enable timely and
accurate analysis of large metagenomics data sets on conventional
desktop computers

[Centrifuge]:     http://www.ccb.jhu.edu/software/centrifuge

[Burrows-Wheeler Transform]: http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[FM Index]:        http://en.wikipedia.org/wiki/FM-index

[GPLv3 license]:   http://www.gnu.org/licenses/gpl-3.0.html


Obtaining Centrifuge
==================

Download Centrifuge and binaries from the Releases sections on the right side.
Binaries are available for Intel architectures (`x86_64`) running Linux, and Mac OS X.

Building from source
--------------------

Building Centrifuge from source requires a GNU-like environment with GCC, GNU Make
and other basics.  It should be possible to build Centrifuge on most vanilla Linux
installations or on a Mac installation with [Xcode] installed.  Centrifuge can
also be built on Windows using [Cygwin] or [MinGW] (MinGW recommended). For a 
MinGW build the choice of what compiler is to be used is important since this
will determine if a 32 or 64 bit code can be successfully compiled using it. If 
there is a need to generate both 32 and 64 bit on the same machine then a multilib 
MinGW has to be properly installed. [MSYS], the [zlib] library, and depending on 
architecture [pthreads] library are also required. We are recommending a 64 bit
build since it has some clear advantages in real life research problems. In order 
to simplify the MinGW setup it might be worth investigating popular MinGW personal 
builds since these are coming already prepared with most of the toolchains needed.

First, download the [source package] from the Releases secion on the right side.
Unzip the file, change to the unzipped directory, and build the
Centrifuge tools by running GNU `make` (usually with the command `make`, but
sometimes with `gmake`) with no arguments.  If building with MinGW, run `make`
from the MSYS environment.

Centrifuge is using the multithreading software model in order to speed up 
execution times on SMP architectures where this is possible. On POSIX 
platforms (like linux, Mac OS, etc) it needs the pthread library. Although
it is possible to use pthread library on non-POSIX platform like Windows, due
to performance reasons Centrifuge will try to use Windows native multithreading
if possible.

[Cygwin]:   http://www.cygwin.com/
[MinGW]:    http://www.mingw.org/
[MSYS]:     http://www.mingw.org/wiki/msys
[zlib]:     http://cygwin.com/packages/mingw-zlib/
[pthreads]: http://sourceware.org/pthreads-win32/
[GnuWin32]: http://gnuwin32.sf.net/packages/coreutils.htm
[Xcode]:    http://developer.apple.com/xcode/
[Github site]: https://github.com/infphilo/centrifuge

Running Centrifuge
=============

Adding to PATH
--------------

By adding your new Centrifuge directory to your [PATH environment variable], you
ensure that whenever you run `centrifuge`, `centrifuge-build` or `centrifuge-inspect`
from the command line, you will get the version you just installed without
having to specify the entire path.  This is recommended for most users.  To do
this, follow your operating system's instructions for adding the directory to
your [PATH].

If you would like to install Centrifuge by copying the Centrifuge executable files
to an existing directory in your [PATH], make sure that you copy all the
executables, including `centrifuge`, `centrifuge-class`, `centrifuge-build`, `centrifuge-build-bin`, `centrifuge-inspect`
and `centrifuge-inspect-bin`. Furthermore you need the files
in the scripts/ folder for genome downloading and compression.

[PATH environment variable]: http://en.wikipedia.org/wiki/PATH_(variable)
[PATH]: http://en.wikipedia.org/wiki/PATH_(variable)

Database building
-----------------

Centrifuge can be build with arbritary sequences. Standard choices are
all complete bacterial and viral genomes, or using the sequences that
are part of the BLAST nt database. Centrifuge always needs the
nodes.dmp file from the NCBI taxonomy dump to build the taxonomy tree,
as well as a sequence ID to taxonomy ID map. The map is a tab-separated
file with the sequence ID to taxonomy ID map.

### Getting taxonomy

    mkdir taxonomy && cd taxonomy
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar xvvf taxdump.tar.gz
    cd ..

### nt database

NCBI BLAST's nt database contains all spliced non-redundant coding
sequences from multiplpe databases, infered from genommic
sequences. Traditionally used with BLAST, a download of the FASTA is
provided on the NCBI homepage. Further we need a sequence ID to taxonomy ID
map that can be generated from a GI taxid dump:

    wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
    gunzip nt.gz && mv -v nt nt.fa

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
    gunzip -c gi_taxid_nucl.dmp.gz | sed 's/^/gi|/' > gi_taxid_nucl.map

    # build index using 16 cores and a small bucket size, which will require less memory
    centrifuge-build -p 16 --bmax 1342177280 --conversion-table gi_taxid_nucl.map --taxonomy-tree taxonomy/nodes.dmp nt.fa nt

### Refseq database

TODO

    scripts/DownloadTaxonomy.sh -l taxonomy
    scripts/DownloadGenomes.sh -g refseq -a 'Complete Genome' -l library bacteria archaea viral

    centrifuge-build --ftabchars 14 -bmax 2147483648 --conversion-table seqid-to-taxid.map --taxonomy-tree nodes.dmp sequences.fa 


Reporting
---------
TODO


Inspecting the database
-----------------------

The database can be inspected with `centrifuge-inspect`.  To extract raw sequences:

    centrifuge-inspect <centrifuge database>

Extract the sequence ID to taxonomy ID conversion table from the database

    centrifuge-inspect --conversion-table <centrifuge database>

Extract the taxonomy tree from the index:

    centrifuge-inspect --taxonomy-tree <centrifuge database>

Extract the lengths of the sequences from the index (each row has two columns: taxonomic ID and length):

    centrifuge-inspect --size-table <centrifuge database>


Wrapper
-------

The `centrifuge`, `centrifuge-build` and `centrifuge-inspect` executables are actually 
wrapper scripts that call binary programs as appropriate. Also, the `centrifuge` wrapper
provides some key functionality, like the ability to handle compressed inputs,
and the functionality for [`--un`], [`--al`] and related options.

It is recommended that you always run the centrifuge wrappers and not run the
binaries directly.

Performance tuning
------------------

1.  If your computer has multiple processors/cores, use `-p NTHREADS`

    The [`-p`] option causes Centrifuge to launch a specified number of parallel
    search threads.  Each thread runs on a different processor/core and all
    threads find alignments in parallel, increasing alignment throughput by
    approximately a multiple of the number of threads (though in practice,
    speedup is somewhat worse than linear).

Command Line
------------


### Usage

    centrifuge [options]* -x <centrifuge-idx> {-1 <m1> -2 <m2> | -U <r>}

### Main arguments

<table><tr><td>

[`-x`]: #centrifuge-options-x

    -x <centrifuge-idx>

</td><td>

The basename of the index for the reference genomes.  The basename is the name of
any of the index files up to but not including the final `.1.cf` / `.rev.1.cf`
/ etc.  `centrifuge` looks for the specified index first in the current directory,
then in the directory specified in the `CENTRIFUGE_INDEXES` environment variable.

</td></tr><tr><td>

[`-1`]: #centrifuge-options-1

    -1 <m1>

</td><td>

Comma-separated list of files containing mate 1s (filename usually includes
`_1`), e.g. `-1 flyA_1.fq,flyB_1.fq`.  Sequences specified with this option must
correspond file-for-file and read-for-read with those specified in `<m2>`. Reads
may be a mix of different lengths. If `-` is specified, `centrifuge` will read the
mate 1s from the "standard in" or "stdin" filehandle.

</td></tr><tr><td>

[`-2`]: #centrifuge-options-2

    -2 <m2>

</td><td>

Comma-separated list of files containing mate 2s (filename usually includes
`_2`), e.g. `-2 flyA_2.fq,flyB_2.fq`.  Sequences specified with this option must
correspond file-for-file and read-for-read with those specified in `<m1>`. Reads
may be a mix of different lengths. If `-` is specified, `centrifuge` will read the
mate 2s from the "standard in" or "stdin" filehandle.

</td></tr><tr><td>

[`-U`]: #centrifuge-options-U

    -U <r>

</td><td>

Comma-separated list of files containing unpaired reads to be aligned, e.g.
`lane1.fq,lane2.fq,lane3.fq,lane4.fq`.  Reads may be a mix of different lengths.
If `-` is specified, `centrifuge` gets the reads from the "standard in" or "stdin"
filehandle.

</td></tr><tr><td>

[`-S`]: #centrifuge-options-S

    -S <hit>

</td><td>

File to write SAM alignments to.  By default, alignments are written to the
"standard out" or "stdout" filehandle (i.e. the console).

</td></tr></table>

### Options

#### Input options

<table>
<tr><td id="centrifuge-options-q">

[`-q`]: #centrifuge-options-q

    -q

</td><td>

Reads (specified with `<m1>`, `<m2>`, `<s>`) are FASTQ files.  FASTQ files
usually have extension `.fq` or `.fastq`.  FASTQ is the default format.  See
also: [`--solexa-quals`] and [`--int-quals`].

</td></tr>
<tr><td id="centrifuge-options-qseq">

[`--qseq`]: #centrifuge-options-qseq

    --qseq

</td><td>

Reads (specified with `<m1>`, `<m2>`, `<s>`) are QSEQ files.  QSEQ files usually
end in `_qseq.txt`.  See also: [`--solexa-quals`] and [`--int-quals`].

</td></tr>
<tr><td id="centrifuge-options-f">

[`-f`]: #centrifuge-options-f

    -f

</td><td>

Reads (specified with `<m1>`, `<m2>`, `<s>`) are FASTA files.  FASTA files
usually have extension `.fa`, `.fasta`, `.mfa`, `.fna` or similar.  FASTA files
do not have a way of specifying quality values, so when `-f` is set, the result
is as if `--ignore-quals` is also set.

</td></tr>
<tr><td id="centrifuge-options-r">

[`-r`]: #centrifuge-options-r

    -r

</td><td>

Reads (specified with `<m1>`, `<m2>`, `<s>`) are files with one input sequence
per line, without any other information (no read names, no qualities).  When
`-r` is set, the result is as if `--ignore-quals` is also set.

</td></tr>
<tr><td id="centrifuge-options-c">

[`-c`]: #centrifuge-options-c

    -c

</td><td>

The read sequences are given on command line.  I.e. `<m1>`, `<m2>` and
`<singles>` are comma-separated lists of reads rather than lists of read files.
There is no way to specify read names or qualities, so `-c` also implies
`--ignore-quals`.

</td></tr>
<tr><td id="centrifuge-options-s">

[`-s`/`--skip`]: #centrifuge-options-s
[`-s`]: #centrifuge-options-s

    -s/--skip <int>

</td><td>

Skip (i.e. do not align) the first `<int>` reads or pairs in the input.

</td></tr>
<tr><td id="centrifuge-options-u">

[`-u`/`--qupto`]: #centrifuge-options-u
[`-u`]: #centrifuge-options-u

    -u/--qupto <int>

</td><td>

Align the first `<int>` reads or read pairs from the input (after the
[`-s`/`--skip`] reads or pairs have been skipped), then stop.  Default: no limit.

</td></tr>
<tr><td id="centrifuge-options-5">

[`-5`/`--trim5`]: #centrifuge-options-5
[`-5`]: #centrifuge-options-5

    -5/--trim5 <int>

</td><td>

Trim `<int>` bases from 5' (left) end of each read before alignment (default: 0).

</td></tr>
<tr><td id="centrifuge-options-3">

[`-3`/`--trim3`]: #centrifuge-options-3
[`-3`]: #centrifuge-options-3

    -3/--trim3 <int>

</td><td>

Trim `<int>` bases from 3' (right) end of each read before alignment (default:
0).

</td></tr><tr><td id="centrifuge-options-phred33-quals">

[`--phred33`]: #centrifuge-options-phred33-quals

    --phred33

</td><td>

Input qualities are ASCII chars equal to the [Phred quality] plus 33.  This is
also called the "Phred+33" encoding, which is used by the very latest Illumina
pipelines.

[Phred quality]: http://en.wikipedia.org/wiki/Phred_quality_score

</td></tr>
<tr><td id="centrifuge-options-phred64-quals">

[`--phred64`]: #centrifuge-options-phred64-quals

    --phred64

</td><td>

Input qualities are ASCII chars equal to the [Phred quality] plus 64.  This is
also called the "Phred+64" encoding.

</td></tr>
<tr><td id="centrifuge-options-solexa-quals">

[`--solexa-quals`]: #centrifuge-options-solexa-quals

    --solexa-quals

</td><td>

Convert input qualities from [Solexa][Phred quality] (which can be negative) to
[Phred][Phred quality] (which can't).  This scheme was used in older Illumina GA
Pipeline versions (prior to 1.3).  Default: off.

</td></tr>
<tr><td id="centrifuge-options-int-quals">

[`--int-quals`]: #centrifuge-options-int-quals

    --int-quals

</td><td>

Quality values are represented in the read input file as space-separated ASCII
integers, e.g., `40 40 30 40`..., rather than ASCII characters, e.g., `II?I`....
 Integers are treated as being on the [Phred quality] scale unless
[`--solexa-quals`] is also specified. Default: off.

</td></tr></table>

#### Alignment options

<table>

<tr><td id="centrifuge-options-n-ceil">

[`--n-ceil`]: #centrifuge-options-n-ceil

    --n-ceil <func>

</td><td>

Sets a function governing the maximum number of ambiguous characters (usually
`N`s and/or `.`s) allowed in a read as a function of read length.  For instance,
specifying `-L,0,0.15` sets the N-ceiling function `f` to `f(x) = 0 + 0.15 * x`,
where x is the read length.  See also: [setting function options].  Reads
exceeding this ceiling are [filtered out].  Default: `L,0,0.15`.

[filtered out]: #filtering

</td></tr>

<tr><td id="centrifuge-options-ignore-quals">

[`--ignore-quals`]: #centrifuge-options-ignore-quals

    --ignore-quals

</td><td>

When calculating a mismatch penalty, always consider the quality value at the
mismatched position to be the highest possible, regardless of the actual value. 
I.e. input is treated as though all quality values are high.  This is also the
default behavior when the input doesn't specify quality values (e.g. in [`-f`],
[`-r`], or [`-c`] modes).

</td></tr>
<tr><td id="centrifuge-options-nofw">

[`--nofw`]: #centrifuge-options-nofw

    --nofw/--norc

</td><td>

If `--nofw` is specified, `centrifuge` will not attempt to align unpaired reads to
the forward (Watson) reference strand.  If `--norc` is specified, `centrifuge` will
not attempt to align unpaired reads against the reverse-complement (Crick)
reference strand. In paired-end mode, `--nofw` and `--norc` pertain to the
fragments; i.e. specifying `--nofw` causes `centrifuge` to explore only those
paired-end configurations corresponding to fragments from the reverse-complement
(Crick) strand.  Default: both strands enabled. 

</td></tr>

</table>

#### Paired-end options

<table>

<tr><td id="centrifuge-options-fr">

[`--fr`/`--rf`/`--ff`]: #centrifuge-options-fr
[`--fr`]: #centrifuge-options-fr
[`--rf`]: #centrifuge-options-fr
[`--ff`]: #centrifuge-options-fr

    --fr/--rf/--ff

</td><td>

The upstream/downstream mate orientations for a valid paired-end alignment
against the forward reference strand.  E.g., if `--fr` is specified and there is
a candidate paired-end alignment where mate 1 appears upstream of the reverse
complement of mate 2 and the fragment length constraints ([`-I`] and [`-X`]) are
met, that alignment is valid.  Also, if mate 2 appears upstream of the reverse
complement of mate 1 and all other constraints are met, that too is valid.
`--rf` likewise requires that an upstream mate1 be reverse-complemented and a
downstream mate2 be forward-oriented. ` --ff` requires both an upstream mate 1
and a downstream mate 2 to be forward-oriented.  Default: `--fr` (appropriate
for Illumina's Paired-end Sequencing Assay).

</td></tr></table>

#### Output options

<table>

<tr><td id="centrifuge-options-t">

[`-t`/`--time`]: #centrifuge-options-t
[`-t`]: #centrifuge-options-t

    -t/--time

</td><td>

Print the wall-clock time required to load the index files and align the reads. 
This is printed to the "standard error" ("stderr") filehandle.  Default: off.

</td></tr>

<!--
<tr><td id="centrifuge-options-un">

[`--un`]: #centrifuge-options-un
[`--un-gz`]: #centrifuge-options-un
[`--un-bz2`]: #centrifuge-options-un

    --un <path>
    --un-gz <path>
    --un-bz2 <path>

</td><td>

Write unpaired reads that fail to align to file at `<path>`.  These reads
correspond to the SAM records with the FLAGS `0x4` bit set and neither the
`0x40` nor `0x80` bits set.  If `--un-gz` is specified, output will be gzip
compressed. If `--un-bz2` is specified, output will be bzip2 compressed.  Reads
written in this way will appear exactly as they did in the input file, without
any modification (same sequence, same name, same quality string, same quality
encoding).  Reads will not necessarily appear in the same order as they did in
the input.

</td></tr>
<tr><td id="centrifuge-options-al">

[`--al`]: #centrifuge-options-al
[`--al-gz`]: #centrifuge-options-al
[`--al-bz2`]: #centrifuge-options-al

    --al <path>
    --al-gz <path>
    --al-bz2 <path>

</td><td>

Write unpaired reads that align at least once to file at `<path>`.  These reads
correspond to the SAM records with the FLAGS `0x4`, `0x40`, and `0x80` bits
unset.  If `--al-gz` is specified, output will be gzip compressed. If `--al-bz2`
is specified, output will be bzip2 compressed.  Reads written in this way will
appear exactly as they did in the input file, without any modification (same
sequence, same name, same quality string, same quality encoding).  Reads will
not necessarily appear in the same order as they did in the input.

</td></tr>
<tr><td id="centrifuge-options-un-conc">

[`--un-conc`]: #centrifuge-options-un-conc
[`--un-conc-gz`]: #centrifuge-options-un-conc
[`--un-conc-bz2`]: #centrifuge-options-un-conc

    --un-conc <path>
    --un-conc-gz <path>
    --un-conc-bz2 <path>

</td><td>

Write paired-end reads that fail to align concordantly to file(s) at `<path>`.
These reads correspond to the SAM records with the FLAGS `0x4` bit set and
either the `0x40` or `0x80` bit set (depending on whether it's mate #1 or #2).
`.1` and `.2` strings are added to the filename to distinguish which file
contains mate #1 and mate #2.  If a percent symbol, `%`, is used in `<path>`,
the percent symbol is replaced with `1` or `2` to make the per-mate filenames.
Otherwise, `.1` or `.2` are added before the final dot in `<path>` to make the
per-mate filenames.  Reads written in this way will appear exactly as they did
in the input files, without any modification (same sequence, same name, same
quality string, same quality encoding).  Reads will not necessarily appear in
the same order as they did in the inputs.

</td></tr>
<tr><td id="centrifuge-options-al-conc">

[`--al-conc`]: #centrifuge-options-al-conc
[`--al-conc-gz`]: #centrifuge-options-al-conc
[`--al-conc-bz2`]: #centrifuge-options-al-conc

    --al-conc <path>
    --al-conc-gz <path>
    --al-conc-bz2 <path>

</td><td>

Write paired-end reads that align concordantly at least once to file(s) at
`<path>`. These reads correspond to the SAM records with the FLAGS `0x4` bit
unset and either the `0x40` or `0x80` bit set (depending on whether it's mate #1
or #2). `.1` and `.2` strings are added to the filename to distinguish which
file contains mate #1 and mate #2.  If a percent symbol, `%`, is used in
`<path>`, the percent symbol is replaced with `1` or `2` to make the per-mate
filenames. Otherwise, `.1` or `.2` are added before the final dot in `<path>` to
make the per-mate filenames.  Reads written in this way will appear exactly as
they did in the input files, without any modification (same sequence, same name,
same quality string, same quality encoding).  Reads will not necessarily appear
in the same order as they did in the inputs.

</td></tr>
-->

<tr><td id="centrifuge-options-quiet">

[`--quiet`]: #centrifuge-options-quiet

    --quiet

</td><td>

Print nothing besides alignments and serious errors.

</td></tr>
<tr><td id="centrifuge-options-met-file">

[`--met-file`]: #centrifuge-options-met-file

    --met-file <path>

</td><td>

Write `centrifuge` metrics to file `<path>`.  Having alignment metric can be useful
for debugging certain problems, especially performance issues.  See also:
[`--met`].  Default: metrics disabled.

</td></tr>
<tr><td id="centrifuge-options-met-stderr">

[`--met-stderr`]: #centrifuge-options-met-stderr

    --met-stderr

</td><td>

Write `centrifuge` metrics to the "standard error" ("stderr") filehandle.  This is
not mutually exclusive with [`--met-file`].  Having alignment metric can be
useful for debugging certain problems, especially performance issues.  See also:
[`--met`].  Default: metrics disabled.

</td></tr>
<tr><td id="centrifuge-options-met">

[`--met`]: #centrifuge-options-met

    --met <int>

</td><td>

Write a new `centrifuge` metrics record every `<int>` seconds.  Only matters if
either [`--met-stderr`] or [`--met-file`] are specified.  Default: 1.

</td></tr>
</table>

#### Performance options

<table><tr>

<td id="centrifuge-options-o">

[`-o`/`--offrate`]: #centrifuge-options-o
[`-o`]: #centrifuge-options-o
[`--offrate`]: #centrifuge-options-o

    -o/--offrate <int>

</td><td>

Override the offrate of the index with `<int>`.  If `<int>` is greater
than the offrate used to build the index, then some row markings are
discarded when the index is read into memory.  This reduces the memory
footprint of the aligner but requires more time to calculate text
offsets.  `<int>` must be greater than the value used to build the
index.

</td></tr>
<tr><td id="centrifuge-options-p">

[`-p`/`--threads`]: #centrifuge-options-p
[`-p`]: #centrifuge-options-p

    -p/--threads NTHREADS

</td><td>

Launch `NTHREADS` parallel search threads (default: 1).  Threads will run on
separate processors/cores and synchronize when parsing reads and outputting
alignments.  Searching for alignments is highly parallel, and speedup is close
to linear.  Increasing `-p` increases Centrifuge's memory footprint. E.g. when
aligning to a human genome index, increasing `-p` from 1 to 8 increases the
memory footprint by a few hundred megabytes.  This option is only available if
`bowtie` is linked with the `pthreads` library (i.e. if `BOWTIE_PTHREADS=0` is
not specified at build time).

</td></tr>
<tr><td id="centrifuge-options-reorder">

[`--reorder`]: #centrifuge-options-reorder

    --reorder

</td><td>

Guarantees that output SAM records are printed in an order corresponding to the
order of the reads in the original input file, even when [`-p`] is set greater
than 1.  Specifying `--reorder` and setting [`-p`] greater than 1 causes Centrifuge
to run somewhat slower and use somewhat more memory then if `--reorder` were
not specified.  Has no effect if [`-p`] is set to 1, since output order will
naturally correspond to input order in that case.

</td></tr>
<tr><td id="centrifuge-options-mm">

[`--mm`]: #centrifuge-options-mm

    --mm

</td><td>

Use memory-mapped I/O to load the index, rather than typical file I/O.
Memory-mapping allows many concurrent `bowtie` processes on the same computer to
share the same memory image of the index (i.e. you pay the memory overhead just
once).  This facilitates memory-efficient parallelization of `bowtie` in
situations where using [`-p`] is not possible or not preferable.

</td></tr></table>

#### Other options

<table>
<tr><td id="centrifuge-options-qc-filter">

[`--qc-filter`]: #centrifuge-options-qc-filter

    --qc-filter

</td><td>

Filter out reads for which the QSEQ filter field is non-zero.  Only has an
effect when read format is [`--qseq`].  Default: off.

</td></tr>
<tr><td id="centrifuge-options-seed">

[`--seed`]: #centrifuge-options-seed

    --seed <int>

</td><td>

Use `<int>` as the seed for pseudo-random number generator.  Default: 0.

</td></tr>
<tr><td id="centrifuge-options-non-deterministic">

[`--non-deterministic`]: #centrifuge-options-non-deterministic

    --non-deterministic

</td><td>

Normally, Centrifuge re-initializes its pseudo-random generator for each read.  It
seeds the generator with a number derived from (a) the read name, (b) the
nucleotide sequence, (c) the quality sequence, (d) the value of the [`--seed`]
option.  This means that if two reads are identical (same name, same
nucleotides, same qualities) Centrifuge will find and report the same alignment(s)
for both, even if there was ambiguity.  When `--non-deterministic` is specified,
Centrifuge re-initializes its pseudo-random generator for each read using the
current time.  This means that Centrifuge will not necessarily report the same
alignment for two identical reads.  This is counter-intuitive for some users,
but might be more appropriate in situations where the input consists of many
identical reads.

</td></tr>
<tr><td id="centrifuge-options-version">

[`--version`]: #centrifuge-options-version

    --version

</td><td>

Print version information and quit.

</td></tr>
<tr><td id="centrifuge-options-h">

    -h/--help

</td><td>

Print usage information and quit.

</td></tr></table>

Centrifuge output
----------

Following is a brief description of the [SAM] format as output by `centrifuge`. 
For more details, see the [SAM format specification][SAM].

The `centrifuge-build` indexer
===========================

`centrifuge-build` builds a Centrifuge index from a set of DNA sequences.
`centrifuge-build` outputs a set of 6 files with suffixes `.1.cf`, `.2.cf`, and
`.3.cf`.  These files together
constitute the index: they are all that is needed to align reads to that
reference.  The original sequence FASTA files are no longer used by Centrifuge
once the index is built.

Use of Karkkainen's [blockwise algorithm] allows `centrifuge-build` to trade off
between running time and memory usage. `centrifuge-build` has two options
governing how it makes this trade: [`--bmax`]/[`--bmaxdivn`],
and [`--dcv`].  By default, `centrifuge-build` will automatically search for the
settings that yield the best running time without exhausting memory.  This
behavior can be disabled using the [`-a`/`--noauto`] option.

The indexer provides options pertaining to the "shape" of the index, e.g.
[`--offrate`](#centrifuge-build-options-o) governs the fraction of [Burrows-Wheeler]
rows that are "marked" (i.e., the density of the suffix-array sample; see the
original [FM Index] paper for details).  All of these options are potentially
profitable trade-offs depending on the application.  They have been set to
defaults that are reasonable for most cases according to our experiments.  See
[Performance tuning] for details.

The Centrifuge index is based on the [FM Index] of Ferragina and Manzini, which in
turn is based on the [Burrows-Wheeler] transform.  The algorithm used to build
the index is based on the [blockwise algorithm] of Karkkainen.

[Blockwise algorithm]: http://portal.acm.org/citation.cfm?id=1314852
[Burrows-Wheeler]: http://en.wikipedia.org/wiki/Burrows-Wheeler_transform
[Performance tuning]: #performance-tuning

Command Line
------------

Usage:

    centrifuge-build [options]* --conversion-table <table_in> --taxonomy-tree <taxonomy_in> <reference_in> <cf_base>

### Main arguments

<table><tr><td>

    <reference_in>

</td><td>

A comma-separated list of FASTA files containing the reference sequences to be
aligned to, or, if [`-c`](#centrifuge-build-options-c) is specified, the sequences
themselves. E.g., `<reference_in>` might be `chr1.fa,chr2.fa,chrX.fa,chrY.fa`,
or, if [`-c`](#centrifuge-build-options-c) is specified, this might be
`GGTCATCCT,ACGGGTCGT,CCGTTCTATGCGGCTTA`.

</td></tr><tr><td>

    <cf_base>

</td><td>

The basename of the index files to write.  By default, `centrifuge-build` writes
files named `NAME.1.cf`, `NAME.2.cf`, and `NAME.3.cf`, where `NAME` is `<cf_base>`.

</td></tr></table>

### Options

<table><tr><td>

    -f

</td><td>

The reference input files (specified as `<reference_in>`) are FASTA files
(usually having extension `.fa`, `.mfa`, `.fna` or similar).

</td></tr><tr><td id="centrifuge-build-options-c">

    -c

</td><td>

The reference sequences are given on the command line.  I.e. `<reference_in>` is
a comma-separated list of sequences rather than a list of FASTA files.

</td></tr>
<tr><td id="centrifuge-build-options-a">

[`-a`/`--noauto`]: #centrifuge-build-options-a

    -a/--noauto

</td><td>

Disable the default behavior whereby `centrifuge-build` automatically selects
values for the [`--bmax`], [`--dcv`] and [`--packed`] parameters according to
available memory.  Instead, user may specify values for those parameters.  If
memory is exhausted during indexing, an error message will be printed; it is up
to the user to try new parameters.

</td></tr><tr><td id="centrifuge-build-options-p">

[`-p`]: #centrifuge-build-options-p

    -p <int>

</td><td>

Launch `NTHREADS` parallel search threads (default: 1).

</td></tr><tr><td id="centrifuge-build-options-conversion-table">

[`--conversion-table`]: #centrifuge-build-options-conversion-table

    --conversion-table <file>

</td><td>

List of UIDs (unique ID) and corresponding taxonomic IDs.

</td></tr><tr><td id="centrifuge-build-options-taxonomy-tree">

[`--taxonomy-tree`]: #centrifuge-build-options-taxonomy-tree

    --taxonomy-tree <file>

</td><td>

Taxonomic tree.

</td></tr><tr><td id="centrifuge-build-options-taxonomy-tree">

[`--size-table`]: #centrifuge-build-options-size-table

    --size-table <file>

</td><td>

List of taxonomic IDs and lengths of the sequences belonging to the same taxonomic IDs.

</td></tr><tr><td id="centrifuge-build-options-bmax">

[`--bmax`]: #centrifuge-build-options-bmax

    --bmax <int>

</td><td>

The maximum number of suffixes allowed in a block.  Allowing more suffixes per
block makes indexing faster, but increases peak memory usage.  Setting this
option overrides any previous setting for [`--bmax`], or [`--bmaxdivn`]. 
Default (in terms of the [`--bmaxdivn`] parameter) is [`--bmaxdivn`] 4.  This is
configured automatically by default; use [`-a`/`--noauto`] to configure manually.

</td></tr><tr><td id="centrifuge-build-options-bmaxdivn">

[`--bmaxdivn`]: #centrifuge-build-options-bmaxdivn

    --bmaxdivn <int>

</td><td>

The maximum number of suffixes allowed in a block, expressed as a fraction of
the length of the reference.  Setting this option overrides any previous setting
for [`--bmax`], or [`--bmaxdivn`].  Default: [`--bmaxdivn`] 4.  This is
configured automatically by default; use [`-a`/`--noauto`] to configure manually.

</td></tr><tr><td id="centrifuge-build-options-dcv">

[`--dcv`]: #centrifuge-build-options-dcv

    --dcv <int>

</td><td>

Use `<int>` as the period for the difference-cover sample.  A larger period
yields less memory overhead, but may make suffix sorting slower, especially if
repeats are present.  Must be a power of 2 no greater than 4096.  Default: 1024.
 This is configured automatically by default; use [`-a`/`--noauto`] to configure
manually.

</td></tr><tr><td id="centrifuge-build-options-nodc">

[`--nodc`]: #centrifuge-build-options-nodc

    --nodc

</td><td>

Disable use of the difference-cover sample.  Suffix sorting becomes
quadratic-time in the worst case (where the worst case is an extremely
repetitive reference).  Default: off.

</td></tr><tr><td>

    -r/--noref

</td><td>

Do not build the `NAME.3.cf` and `NAME.4.cf` portions of the index, which
contain a bitpacked version of the reference sequences and are used for
paired-end alignment.

</td></tr><tr><td>

    -3/--justref

</td><td>

Build only the `NAME.3.cf` and `NAME.4.cf` portions of the index, which
contain a bitpacked version of the reference sequences and are used for
paired-end alignment.

</td></tr><tr><td id="centrifuge-build-options-o">

    -o/--offrate <int>

</td><td>

To map alignments back to positions on the reference sequences, it's necessary
to annotate ("mark") some or all of the [Burrows-Wheeler] rows with their
corresponding location on the genome. 
[`-o`/`--offrate`](#centrifuge-build-options-o) governs how many rows get marked:
the indexer will mark every 2^`<int>` rows.  Marking more rows makes
reference-position lookups faster, but requires more memory to hold the
annotations at runtime.  The default is 4 (every 16th row is marked; for human
genome, annotations occupy about 680 megabytes).  

</td></tr><tr><td>

    -t/--ftabchars <int>

</td><td>

The ftab is the lookup table used to calculate an initial [Burrows-Wheeler]
range with respect to the first `<int>` characters of the query.  A larger
`<int>` yields a larger lookup table but faster query times.  The ftab has size
4^(`<int>`+1) bytes.  The default setting is 10 (ftab is 4MB).

</td></tr><tr><td>

    --seed <int>

</td><td>

Use `<int>` as the seed for pseudo-random number generator.

</td></tr><tr><td>

    --cutoff <int>

</td><td>

Index only the first `<int>` bases of the reference sequences (cumulative across
sequences) and ignore the rest.

</td></tr><tr><td>

    -q/--quiet

</td><td>

`centrifuge-build` is verbose by default.  With this option `centrifuge-build` will
print only error messages.

</td></tr><tr><td>

    -h/--help

</td><td>

Print usage information and quit.

</td></tr><tr><td>

    --version

</td><td>

Print version information and quit.

</td></tr></table>

The `centrifuge-inspect` index inspector
=====================================

`centrifuge-inspect` extracts information from a Centrifuge index about what kind of
index it is and what reference sequences were used to build it. When run without
any options, the tool will output a FASTA file containing the sequences of the
original references (with all non-`A`/`C`/`G`/`T` characters converted to `N`s).
 It can also be used to extract just the reference sequence names using the
[`-n`/`--names`] option or a more verbose summary using the [`-s`/`--summary`]
option.

Command Line
------------

Usage:

    centrifuge-inspect [options]* <cf_base>

### Main arguments

<table><tr><td>

    <cf_base>

</td><td>

The basename of the index to be inspected.  The basename is name of any of the
index files but with the `.X.cf` suffix omitted.
`centrifuge-inspect` first looks in the current directory for the index files, then
in the directory specified in the `Centrifuge_INDEXES` environment variable.

</td></tr></table>

### Options

<table><tr><td>

    -a/--across <int>

</td><td>

When printing FASTA output, output a newline character every `<int>` bases
(default: 60).

</td></tr><tr><td id="centrifuge-inspect-options-n">

[`-n`/`--names`]: #centrifuge-inspect-options-n

    -n/--names

</td><td>

Print reference sequence names, one per line, and quit.

</td></tr><tr><td id="centrifuge-inspect-options-s">

[`-s`/`--summary`]: #centrifuge-inspect-options-s

    -s/--summary

</td><td>

Print a summary that includes information about index settings, as well as the
names and lengths of the input sequences.  The summary has this format: 

    Colorspace	<0 or 1>
    SA-Sample	1 in <sample>
    FTab-Chars	<chars>
    Sequence-1	<name>	<len>
    Sequence-2	<name>	<len>
    ...
    Sequence-N	<name>	<len>

Fields are separated by tabs.  Colorspace is always set to 0 for Centrifuge.

</td></tr><tr><td id="centrifuge-inspect-options-conversion-table">

[`--conversion-table`]: #centrifuge-inspect-options-conversion-table

    --conversion-table

</td><td>

Print a list of UIDs (unique ID) and corresponding taxonomic IDs.

</td></tr><tr><td id="centrifuge-inspect-options-taxonomy-tree">

[`--taxonomy-tree`]: #centrifuge-inspect-options-taxonomy-tree

    --taxonomy-tree

</td><td>

Print a taxonomic tree.

</td></tr><tr><td id="centrifuge-inspect-options-taxonomy-tree">

[`--size-table`]: #centrifuge-inspect-options-size-table

    --size-table

</td><td>

Print a list of taxonomic IDs and lengths of the sequences belonging to the same taxonomic IDs.

</td></tr><tr><td>

    -v/--verbose

</td><td>

Print verbose output (for debugging).

</td></tr><tr><td>

    --version

</td><td>

Print version information and quit.

</td></tr><tr><td>

    -h/--help

</td><td>

Print usage information and quit.

</td></tr></table>

Getting started with Centrifuge
===================================================

Centrifuge comes with some example files to get you started.  The example files
are not scientifically significant; these files will simply let you start running Centrifuge and
downstream tools right away.

First follow the manual instructions to [obtain Centrifuge].  Set the `CENTRIFUGE_HOME`
environment variable to point to the new Centrifuge directory containing the
`centrifuge`, `centrifuge-build` and `centrifuge-inspect` binaries.  This is important,
as the `CENTRIFUGE_HOME` variable is used in the commands below to refer to that
directory.

[obtain Centrifuge]: #obtaining-centrifuge

Indexing a reference genome
---------------------------

To create an index for two small sequences included with Centrifuge, create a new temporary directory (it doesn't matter where), change into that directory, and run:

    $CENTRIFUGE_HOME/centrifuge-build --conversion-table $CENTRIFUGE_HOME/example/reference/gi_to_tid.dmp --taxonomy-tree $CENTRIFUGE_HOME/example/reference/taxonomy.dmp $CENTRIFUGE_HOME/example/reference/test.fa test

The command should print many lines of output then quit. When the command
completes, the current directory will contain ten new files that all start with
`test` and end with `.1.cf`, `.2.cf`, `.3.cf`.  These files constitute the index - you're done!

You can use `centrifuge-build` to create an index for a set of FASTA files obtained
from any source, including sites such as [UCSC], [NCBI], and [Ensembl]. When
indexing multiple FASTA files, specify all the files using commas to separate
file names.  For more details on how to create an index with `centrifuge-build`,
see the [manual section on index building].  You may also want to bypass this
process by obtaining a pre-built index.

[UCSC]: http://genome.ucsc.edu/cgi-bin/hgGateway
[NCBI]: http://www.ncbi.nlm.nih.gov/sites/genome
[Ensembl]: http://www.ensembl.org/
[manual section on index building]: #the-centrifuge-build-indexer
[using a pre-built index]: #using-a-pre-built-index

Classifying example reads
----------------------

Stay in the directory created in the previous step, which now contains the
`test` index files.  Next, run:

    $CENTRIFUGE_HOME/centrifuge -f -x test $CENTRIFUGE_HOME/example/reads/input.fa

This runs the Centrifuge classifier, which classifies a set of unpaired reads to the
the genomes using the index generated in the previous step.
The classification results are reported to stdout, and a
short classification summary is written to centrifuge-species_report.csv.

You will see something like this:

    readID  uniqueID taxID     score	2ndBestScore	hitLength	numMatches
    C_1 gi|7     9913      4225	4225		80.00		2
    C_1 gi|4     9646      4225	4225		80.00		2
    C_2 gi|4     9646      4225	4225		80.00		2
    C_2 gi|7     9913      4225	4225		80.00		2
    C_3 gi|7     9913      4225	4225		80.00		2
    C_3 gi|4     9646      4225	4225		80.00		2
    C_4 gi|4     9646      4225	4225		80.00		2
    C_4 gi|7     9913      4225	4225		80.00		2
    1_1 gi|4     9646      4225	0		80.00		1
    1_2 gi|4     9646      4225	0		80.00		1
    2_1 gi|7     9913      4225	0		80.00		1
    2_2 gi|7     9913      4225	0		80.00		1
    2_3 gi|7     9913      4225	0		80.00		1
    2_4 gi|7     9913      4225	0		80.00		1
    2_5 gi|7     9913      4225	0		80.00		1
    2_6 gi|7     9913      4225	0		80.00		1
