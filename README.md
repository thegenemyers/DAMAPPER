# The Damapper Library

## _Author:  Gene Myers_
## _First:   July 11, 2016_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/damapper-commands).

The commands below permit one to find the best, or several good alternative,
locations that a read may have been sampled from in a reference genome.  The
code is designed for PACBIO reads with an assumption of a 15% average error
rate and occational large drop out regions where the quality is exceedingly
low.

```
1. damapper [-vbpCN] [-k<int(20)>] [-t<int>] [-M<int>] [-T<int(4)>]
                     [-e<double(.85)] [-s<int(100)>] [-n<double(1.)>]
                     [-m<track>]+ <ref:db|dam> <reads:db|dam> ...
```

Search the reference data base \<ref\> for the best matches to each read in the list of
databases or database blocks \<reads\>.  The parameters -v, -b, -k, -t, -M, -T -e, -s,
and -m are exactly as for the daligner (see here).  Matches at the expected correlation
given by -e are sought and unlike daligner they can be of any length (i.e. short).
Alignments at much lower rates are not forced in order to report a single matching
interval, but rather when a read has bad stretches, the result is a chain of local
alignments over the good segments of the read.  Damapper reports the best chain
covering a segment of a read, and may report several chains covering disjoint segments
(e.g. the case of a chimeric read).  If the -n option is given then all chains that
are within given fraction of optimal are also reported, e.g. -n.95 reports all matches
within 95% of the top match.

For a given reference data base X and read block Y, damapper produces the single file
Y.X..las.
Note the order of X and Y is reversed
from that on the command line (which is different than the conventions for
daligner).
The reversal of the order of the blocks emphasizes that in the result, the A-reads
are the mapped reads, and the B-reads are contigs of the reference.  Each file is
sorted in order of the A-reads, and if a match is a chain of local alignments the
la's in the chain occur in increasing order of A-coordinate.  Moreover, if several
near optimal matches are reported, matches appear in order of score starting
with the best.

Itf the -C option is set, then damapper also outputs a file X.Y.las for a given
block pair that contains all
the same matches as in X.Y.las but where the A-read is a contig of the
reference and the B-read is a mapped read.  And if the -N options is set, then the
files Y.X.las is not produced.

The -p option requests that damapper produce a repeat profile track for each read.
For each trace-point sized interval (established by the -s parameter) the number of
times, c, this segment is involved in a distinct alignment to the reference is
estimated.  0 is recorded for segments that don't match anything in the reference,
1 for segments that match uniquely, and otherwise floor(log<sub>10</sub>c/10) up to a cap of
40 corresponding to 10,000 copies.  Observe carefully that the track has the same
form as intrinsic quality values.  They can be output by DBdump and DaViewer is able
to graphically display them.  These profiles should make it obvious when a read does
not have a unique location in a reference sequence due to its being entirely or almost
entirely repetitive.

```
2. HPC.damapper [-vbpCN]
                [-k<int(20)>] [-t<int>] [-M<int>] [-e<double(.85)] [-s<int(100)]
                [-n<double(1.)>] [-m<track>]+ [-B<int(4)>] [-T<int(4)>] [-f<name>]
                <ref:db|dam> <reads:db|dam> [<first:int>[-<last:int>]]
```

HPC.damapper writes a UNIX shell script to the standard output that maps every read in
blocks \<first\> to \<last\> of database \<reads\> to a reference sequence \<ref\>.  If \<last\>
is missing then only the single block \<first\> is mapped, and if \<first\> is also missing
then all blocks of the database are mapped.  Except for the -B and -f options, all
other options are passed through to damapper, save of -v which is passed to every
program in the script.  The -B option determines the maximum number of blocks mapped
per call to damapper in the script.

The command script output by HPC.damapper and other HPC.\<x\> programs consists of
command blocks each of which begins with a comment line (begins with #) followed by a
potentially long list of lines each containing a shell command.  Command blocks whose
comment mentions "jobs" and gives the number of said in parenthesis, we call parallel
blocks because each command line in the block can be sent to a node in a cluster for
independent execution, i.e. none of the commands in a block depend on another in the
block.  The remaining command blocks we call house-keeping blocks because they can be
executed by the shell on the launch/server node and the commands are either checking
the integrity of .las files with LAcheck, or removing intermediate files with rm. Each
block should be performed in the order given and should complete before the next
block is performed.

If the -f option is set, then each command block is written to a file with a name of
the form \<name\>.#.\<description\> where \<name\> is specified by the user in the -f option
argument, # gives the order in which the command block in the given file is to be
performed in relation to other command block files, and \<description\> is a (very)
short symbolic reminder of what the block is doing.  For example,
"HPC.damapper -fJOBS REF DB" would produce the files:

```
    JOBS.01.OVL
    JOBS.02.CHECK.OPT
```

The files with the suffix .OPT are optional and need not be executed albeit we highly
recommend that one run all the CHECK blocks.
