# Whippet

[![Latest Release](https://img.shields.io/github/release/timbitz/Whippet.jl.svg)](https://github.com/timbitz/Whippet.jl/releases/latest)
[![Whippet](http://pkg.julialang.org/badges/Whippet_0.6.svg)](http://pkg.julialang.org/detail/Whippet)
[![Build Status](https://travis-ci.com/timbitz/Whippet.jl.svg?token=R7mZheNGhsReQ7hn2gdf&branch=master)](https://travis-ci.com/timbitz/Whippet.jl)
[![codecov](https://codecov.io/gh/timbitz/Whippet.jl/branch/master/graph/badge.svg?token=RKE7BSr47v)](https://codecov.io/gh/timbitz/Whippet.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Gitter Chat](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/Whippet-jl/Lobby)

## Graphical Overview
![Whippet Schematic](https://timbitz.github.io/images/Whippet-Fig1.gif)

## Features
- Splice graph representations of transcriptome structure
  - Build an index for any species with a genome and annotation file
  - _de novo_ event discovery (splicing from/to any combination of annotated donor/acceptor splice sites)
- High speed PolyA+ Spliced Read Alignment (Read lengths <= 255)
  - Repetitive read assignment for gene families
- Fast and robust quantification of transcriptome structure and expression using EM
  - Dynamic building and entropic measurements of splicing events of any complexity
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
- Differential splicing comparisons
  - Probabilistic calculations of delta PSI leveraging multi-sample biological replicates

Paper pre-print: http://www.biorxiv.org/content/early/2017/07/03/158519

## Get started

### 1) Install 
Get latest version of [Julia](https://julialang.org/downloads) if you don't have it.  If you are new to julia, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started)

Install Whippet through the Julia REPL:
```julia
julia> Pkg.add("Whippet")
julia> using Whippet
```

The Whippet package directory can now be found in the default location:
```bash
$ cd ~/.julia/v0.6/Whippet
```

_Notes_:
* If you are having trouble finding the Whippet directory, you can ask the Julia REPL: `julia> Pkg.dir("Whippet")`
* For all executables in `Whippet/bin`, you can use the `-h` flag to get a list of the available command line options, their usage and defaults.
* You can always update to the latest version of Whippet using `Pkg.update()` in the Julia REPL!
* You should install Julia and its packages locally, if you absolutely have to install system-wide, there is some help [here](https://groups.google.com/forum/#!topic/julia-users/9lQZJlLs99M) 
* It might be convenient to add a link (`ln -s ~/.julia/v0.6/Whippet`) to this directory for easy access, or export `Whippet/bin` to your path.
 
### 2) Build an index.

You need your genome sequence in fasta, and a gene annotation file in GTF or Refflat format. Default annotation supplied for hg19 GENCODEv25 TSL1-level transcriptome in `Whippet/anno`. 

```bash
$ julia bin/whippet-index.jl --fasta hg19.fa.gz --gtf anno/gencode_hg19.v25.tsl1.gtf.gz
```

_Notes_: 
* Whippet only uses GTF `exon` lines (others are ignored). These must contain both `gene_id` and `transcript_id` attributes (which should not be the same as one another!).  This GTF file should be consistent with the [GTF2.2](http://mblab.wustl.edu/GTF22.html) specification, and should have all entries for a transcript in a continuous block. 
* You can specify the output name and location of the index to build using the `-x / --index` parameter. The default (for both whippet-index.jl and whippet-quant.jl) is a generic index named `graph` located at `~/.julia/v0.6/Whippet/index/graph.jls`, so you must have write-access to this location to use the default.


### 3) Quantify FASTQ files.

#### a) Single-end reads
```bash
$ julia bin/whippet-quant.jl file.fastq.gz
```
Or you can provide a link/s to the fastq file/s on the web using `--url`, or just the experiment accession id (SRR id) using `--ebi`!  These will align on the fly as the reads are downloaded from the web directly into alignment. For example:
```bash
$ julia bin/whippet-quant.jl --ebi SRR1199010
```
is equivalent to
```bash
$ julia bin/whippet-quant.jl --url ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/000/SRR1199010/SRR1199010.fastq.gz
```

#### b) Paired-end reads
```bash
$ julia bin/whippet-quant.jl fwd_file.fastq.gz rev_file.fastq.gz
```
To stream reads from ebi.ac.uk, just provide a paired-end SRR id to to the `--ebi` flag:
```bash
$ julia bin/whippet-quant.jl --ebi SRR1658024
```
is equivalent to
```bash
$ julia bin/whippet-quant.jl --url http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/004/SRR1658024/SRR1658024_1.fastq.gz http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/004/SRR1658024/SRR1658024_2.fastq.gz
```

#### c) Non-default input/output
To specify output location or a specific index:
```bash
$ julia bin/whippet-quant.jl fwd_file.fastq.gz -o outputname -x customindex.jls
```

You can also output the alignments in SAM format with the `--sam` flag and convert to bam with a pipe:
```bash
$ julia bin/whippet-quant.jl fwd_file.fastq.gz --sam | samtools view -bS - > fwd_file.bam
```

It is also possible to pool fastq files at runtime using shell commands, and the optional (`--force-gz`) for pooled gz files (files without .gz suffix)
```bash
$ julia bin/whippet-quant.jl <( cat time-series_{1,2,3,4,5}.fastq.gz ) --force-gz -o interval_1-5
```

### 4) Compare multiple psi files
Compare `.psi.gz` files from from two samples `-a` and `-b` with any number of replicates (comma delimited list of files or common pattern matching) per sample.
```bash
$ ls *.psi.gz
sample1-r1.psi.gz sample1-r2.psi.gz sample2-r1.psi.gz sample2-r2.psi.gz
$ julia bin/whippet-delta.jl -a sample1 -b sample2
OR
$ julia bin/whippet-delta.jl -a sample1-r1.psi.gz,sample1-r2.psi.gz -b sample2-r1.psi.gz,sample2-r2.psi.gz
```
Note: comparisons of single files still need a comma: `-a singlefile_a.psi.gz, -b singlefile_b.psi.gz,`

---

## Output Formats

The output format for `whippet-quant.jl` is saved into two core quant filetypes, `.psi.gz` and `.tpm.gz` files.

Each `.tpm.gz` file contains a simple format compatible with many downstream tools (one for the TPM of each annotated isoform, and another at the gene-level):

Gene/Isoform | TpM | Read Counts
---- | --- | -----------
NFIA | 2897.11 | 24657.0

---
Meanwhile the `.psi.gz` file is a bit more complex and requires more explanation. Quantified nodes can fall into a number of "node-type" categories:

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
 BS  | Circular back-splicing (output only when `--circ` flag is given)

Each node is defined by a type (above) and has a corresponding value for `Psi` or the Percent-Spliced-In.

The `.psi.gz` file itself is outputted in the form:

Gene | Node | Coord | Strand | Type | Psi | CI Width | CI Lo,Hi | Total Reads | Complexity | Entropy | Inc Paths | Exc Paths
---- | ---- | ----- | ------ | ---- | --- | -------- | -------- | ----------- | ---------- | ------- | --------- | ---------
NFIA | 2 | chr1:61547534-61547719 | + | AF | 0.782 | 0.191 | 0.669,0.86 | 49.0 | K1 | 0.756 | 2-4-5:0.782 | 1-5:0.218
NFIA | 4 | chr1:61548433-61548490 | + | CE | 0.8329 | 0.069 | 0.795,0.864 | 318.0 | K2 | 1.25 | 2-4-5:0.3342,3-4-5:0.4987 | 1-5:0.1671
NFIA | 5 | chr1:61553821-61554352 | + | CE | 0.99 | NA | NA | NA | NA | NA | NA | NA
NFIA | 6 | chr1:61743192-61743257 | + | CE | 0.99 | NA | NA | NA | NA | NA | NA | NA

Whippet outputs quantification at the node-level (one PSI value for each node, not exon). For example, Whippet CE ('core exon') nodes may be bounded by one or more AA or AD nodes. So when those AA or AD nodes are spliced in, the full exon consists of the nodes combined (AA+CE). 

We believe that this is (and should be) the correct generalization of event-level splicing quantification for the following reasons:

1. Whippet allows for dynamic quantification of observed splicing patterns. Therefore, in order to maintain consistent output from different Whippet runs on various samples, the basic unit of quantification needs to be a `node`. 
2. Since it is possible for the entire exon (CE+AA both with or without the AA) to be excluded from the mature message, the PSI of each node should be given independently. Whippet can handle highly complex events, in contrast to many other splicing quantification tools which only report binary event types. If Whippet output was enumerated for all such combinations, the output for complex events would quickly become exponential and approach uselessness.
3. The general definition of Percent-spliced-in for an exon (or node) should be the percentage of transcripts that 'have that exon spliced in', irrespective of the upstream or downstream splice sites that connect to it (those merely alter another node's PSI). Therefore, we feel that the output of PSI values should not change based on upstream or downstream splice-sites (as they do with some other programs).

In conclusion, this must be taken into account when intersecting `.psi.gz` coordinate output with other formats that represent full exons (which can be one or more adjacent nodes combined). 

Complexity refers to the discrete categories based-on the log2(number of paths) through each local splicing event (LSE). Entropy refers to the shannon-entropy of the relative expression of the paths through the LSE.

---
The output format of `.diff.gz` files from `whippet-delta.jl` is:

Gene | Node | Coord | Strand | Type | Psi_A | Psi_B | DeltaPsi | Probability | Complexity | Entropy
---- | ---- | ----- | ------ | ---- | ----- | ----- | -------- | ----------- | ---------- | -------
ENSG00000117448.13_2 | 9 | chr1:46033654-46033849 | + | CE | 0.95971 | 0.97876 | -0.019047 | 0.643 | K0 | 0.0   
ENSG00000117448.13_2 | 10 | chr1:46034157-46034356 | + | CE | 0.9115 | 0.69021 | 0.22129 | 0.966 | K1 | 0.874 

Between the set of replicates from -a and -b `whippet-delta.jl` outputs a mean Psi_A (from emperical sampling of the joint posterior distribution) and mean Psi_B.  It also outputs the mean deltaPsi (Psi_A - Psi_B) and the probability that there is some change in deltaPsi, such that P(|deltaPsi| > 0.0). To determine the significance of output in `.diff.gz` files, two filters are strongly encouraged: (1) significant nodes should have a high probability (>=0.95) suggesting there is likely to be a difference between the samples, and (2) significant nodes should have an absolute magnitude of change (|DeltaPsi| > x) where x is some value you consider to be a biologically relevant magnitude of change (we suggest |DeltaPsi| > 0.1).

---

## Index Building Strategies

If you are building an index for another organism, there are some general guidelines that can help to ensure that the index you build is as effective as it can be. In general you should seek to:
  * Use only the highest quality annotations you can find.
  * If a high-quality annotation set like GENCODE, contains both high and low quality transcripts, filter it! For human, we only use TSL-1 annotations.
  * Avoid giving annotations with 'indels' such as ESTs or mRNAs without filtering out invalid splice sites first.
  * If you plan to align very short reads (~36nt), decrease the Kmer size (we have used 6nt before), otherwise the default Kmer size (9nt) should be used.

---

## Troubleshooting

With all of the executables in `Whippet/bin`, you can use the `-h` flag to get a list of the available command line options and their usage.  If you are having trouble using or interpreting the output of `Whippet` then please ask a question in our [gitter chat](https://gitter.im/Whippet-jl/Lobby)!.  If you are having trouble running a Whippet executable in /bin, try running `Pkg.test("Whippet")` from the Julia REPL, Whippet should run cleanly if your environment is working correctly (though it is untested on Windows, it should still work).  If you still think you have found a bug feel free to open an issue in github or make a pull request! 
