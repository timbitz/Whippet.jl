# Whippet

[![Build Status](https://travis-ci.com/timbitz/Whippet.svg?token=R7mZheNGhsReQ7hn2gdf&branch=master)](https://travis-ci.com/timbitz/Whippet)
[![codecov](https://codecov.io/gh/timbitz/Whippet/branch/master/graph/badge.svg?token=RKE7BSr47v)](https://codecov.io/gh/timbitz/Whippet)


##### Ultra fast & lightweight quantification of gene expression and event-specific splicing levels from RNA-seq.

## Features
- High performance PolyA+ Spliced Read Alignment
  - Repetitive read assignment for gene families
- Robust quantification of the expression and transcriptome structure of model and non-model organisms
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
- Accurate splice graph representations of high complexity event types (splicing and alt-3'/5' end usage)
  - Pseudo _de novo_ event discovery
  - Entropic measurements of Splicing-event Complexity
  - Circular splicing discovery



### 1) Install 
Install most recent [julia release here](http://julialang.org/downloads/), which must be v0.4.x! (v0.5 will be supported in future releases).  If you are new to julia, or installing programs via command line, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started)

### 2) Clone Whippet
Make sure dependencies are satisfied. Executables are in bin/
```
git clone https://github.com/timbitz/Whippet.git
cd Whippet
julia dependencies.jl
```
NOTE: `julia dependencies.jl` may be noisy with deprecated syntax warnings.  This is due to the rapid pace at which base julia is being developed and does not actually mean that there was/is a fatal problem with Whippet or its dependencies. 

### 3) Build an index.  
You need your genome sequence in fasta, and a gene annotation file in GTF or Refflat format. Default examples are supplied for hg19.

```bash
$ julia whippet-index.jl --fasta hg19.fa.gz --flat anno/refseq_hg19.flat.gz
```

### 4) Quantify FASTQ files.
```bash
$ julia whippet-quant.jl file.fastq.gz
```

Or if you have paired-end RNA-seq data...
```bash
$ julia whippet-quant.jl fwd_file.fastq.gz rev_file.fastq.gz
```

You can output the alignments in SAM format with the `--sam` flag and convert to bam with a pipe:
```bash
$ julia whippet-quant.jl fwd_file.fastq.gz --sam | samtools view -bS - > fwd_file.bam
```

It is also possible to pool fastq files at runtime using shell commands, and the optional (`--force-gz`) for pooled gz files (files without .gz suffix)
```bash
$ julia whippet-quant.jl <( cat SRR208080{1,2,3,4,5,6,7,8,9}.fastq.gz ) --force-gz -o SRR208080_1-9
```

### 5) Compare multiple psi files
```bash
$ ls *.psi.gz
sample1-r1.psi.gz sample1-r2.psi.gz sample2-r1.psi.gz sample2-r2.psi.gz
$ julia whippet-delta.jl -a sample1 -b sample2
OR
$ julia whippet-delta.jl -a sample1-r1.psi.gz,sample1-r2.psi.gz -b sample2-r1.psi.gz,sample2-r2.psi.gz
```

---

## Output Formats

The output format for `whippet-quant.jl` is saved into two files a `.psi.gz` and a `.tpm.gz`.

The `.tpm.gz` file contains a simple format compatible with many downstream tools:

Gene | TpM | Read Counts
---- | --- | -----------
NFIA | 2897.11 | 24657.0

Meanwhile the `.psi.gz` file is a bit more complex and requires more explanation:

Gene | Node | Coord | Strand | Type | Psi | CI Width | CI Lo,Hi | Total Reads | Complexity | Entropy | Inc Paths | Exc Paths
---- | ---- | ----- | ------ | ---- | --- | -------- | -------- | ----------- | ---------- | ------- | --------- | ---------
NFIA | 2 | chr1:61547534-61547719 | + | AF | 0.782 | 0.191 | 0.669,0.86 | 49.0 | K1 | 0.756 | IntSet([2, 4, 5]) | IntSet([1, 5])
NFIA | 4 | chr1:61548433-61548490 | + | CE | 0.8329 | 0.069 | 0.795,0.864 | 318.0 | K2 | 1.25 | IntSet([2, 4, 5]),IntSet([3, 4, 5]) | IntSet([1, 5])
NFIA | 5 | chr1:61553821-61554352 | + | CE | 0.99 | NA | NA | NA | NA | NA | NA | NA
NFIA | 6 | chr1:61743192-61743257 | + | CE | 0.99 | NA | NA | NA | NA | NA | NA | NA

In contrast to many other splicing quantification tools, Whippet allows for dynamic quantification of observed splicing patterns.  Therfore in order to maintain consistent output from different Whippet runs on various samples, the basic unit of quantification is a SpliceGraph `node`.  It is possible (and even likely) that many nodes are never spliced entirely on their own as is the case with alternative 5' and 3' splice sites, and core exon nodes whose neighboring alt 5' or 3' splice sites are used.  Therefore this must be taken into account when intersecting `.psi.gz` coordinate output with other formats that represent full exons (which can be one or more adjacent nodes combined).


Type | Interpretation
---- | --------------
 CE  | Core exon, which may be bounded by one or more alternative AA/AD nodes
 AA  | Alternative Acceptor splice site
 AD  | Alternative Donor splice site
 RI  | Retained intron
 TS  | Tandem transcription start site
 TE  | Tandem alternative polyadenylation site
 AF  | Alternative First exon
 AL  | Alternative Last exon
 
Each node is defined by a type (above) and has a corresponding value for `Psi` or the Percent-Spliced-In followed by the 90% confidence interval (both the width as well as lower and higher boundaries).

## Advanced Index Building

If you are building an index for a non-model organism or an index for a custom purpose, there are some general guidelines that can help to ensure that the index you build is as effective as it can be.  For example, a whippet index that is missing many annotated splice sites that are frequently used, may not be able to align all reads well. Similarly, a whippet index that is built with a huge number of decoy splice sites from indels in EST or mRNA annotations, may have too many `splice sites`, which will divide exons into many tiny nodes, making seeding to those segments more difficult. In general you should seek to:
  * Increase the number of true splice sites that `whippet index` is given. 
  * Avoid giving annotations with 'indels' such as ESTs or mRNAs without filtering for valid splice sites first.
  * Choose a Kmer size appropriate for the node size and number in the transcriptome and the read length you plan to align.

```bash
$ julia whippet-index.jl -h
Whippet v0.3 loading and compiling... 
usage: whippet-index.jl [-k KMER] --fasta FASTA [--flat FLAT]
                        [--gtf GTF] [--index INDEX] [-h]

optional arguments:
  -k, --kmer KMER  Kmer size to use for exon-exon junctions (default
                   9) (type: Int64, default: 9)
  --fasta FASTA    File containg the genome in fasta, one entry per
                   chromosome [.gz]
  --flat FLAT      Gene annotation file in RefFlat format
  --gtf GTF        Gene anotation file in GTF format
  --index INDEX    Output prefix for saving index 'dir/prefix'
                   (default Whippet/index/graph) (default:
                   "/path/to/Whippet/src/../index/graph")
  -h, --help       show this help message and exit
```

## Custom Alignment Parameters

```bash
$ julia whippet-quant.jl -h
Whippet v0.3 loading and compiling... 
usage: whippet-quant.jl [-x INDEX] [-o OUT] [-s] [-L SEED-LEN]
                        [-M SEED-TRY] [-T SEED-TOL] [-B SEED-BUF]
                        [-I SEED-INC] [-P PAIR-RANGE] [-X MISMATCHES]
                        [-S SCORE-MIN] [--psi-body-read] [--stranded]
                        [--pair-same-strand] [--phred-33] [--phred-64]
                        [--no-circ] [--no-tpm] [--force-gz] [-h]
                        filename.fastq[.gz] [paired_mate.fastq[.gz]]

positional arguments:
  filename.fastq[.gz]
  paired_mate.fastq[.gz]


optional arguments:
  -x, --index INDEX     Output prefix for saving index 'dir/prefix'
                        (default Whippet/index/graph) (default:
                        "/path/to/Whippet/index/graph")
  -o, --out OUT         Where should the gzipped output go
                        'dir/prefix'? (default:
                        "/path/to/Whippet/output")
  -s, --sam             Should SAM format be sent to stdout?
  -L, --seed-len SEED-LEN
                        Seed length (type: Int64, default: 18)
  -M, --seed-try SEED-TRY
                        Number of failed seeds to try before giving up
                        (type: Int64, default: 3)
  -T, --seed-tol SEED-TOL
                        Number of seed hits to tolerate (type: Int64,
                        default: 4)
  -B, --seed-buf SEED-BUF
                        Ignore this many bases from beginning and end
                        of read for seed (type: Int64, default: 5)
  -I, --seed-inc SEED-INC
                        Number of bases to increment seed each
                        iteration (type: Int64, default: 18)
  -P, --pair-range PAIR-RANGE
                        Seeds for paired end reads must match within _
                        bases of one another (type: Int64, default:
                        2500)
  -X, --mismatches MISMATCHES
                        Allowable number of mismatches in alignment
                        (counted as 1-10^(-phred/10)) (type: Int64,
                        default: 3)
  -S, --score-min SCORE-MIN
                        Minimum percent matching (matches -
                        mismatches) / read_length (type: Float64,
                        default: 0.6)
  --psi-body-read       Allow exon-body reads in quantification of PSI
                        values
  --stranded            Is the data strand specific in fwd
                        orientation? If so, increase speed with this
                        flag
  --pair-same-strand    Whippet by default tries to align fwd/rev
                        pairs, if your data is fwd/fwd or rev/rev set
                        this flag
  --phred-33            Qual string is encoded in Phred+33 integers
                        (default)
  --phred-64            Qual string is encoded in Phred+64 integers
  --no-circ             Do not allow back/circular splicing
  --no-tpm              Should tpm file be sent to
                        output/prefix.tpm.gz? (default on)
  --force-gz            Regardless of suffix, consider read input as
                        gzipped
  -h, --help            show this help message and exit

```
