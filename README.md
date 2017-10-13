# Whippet

[![Build Status](https://travis-ci.com/timbitz/Whippet.jl.svg?token=R7mZheNGhsReQ7hn2gdf&branch=master)](https://travis-ci.com/timbitz/Whippet.jl)
[![codecov](https://codecov.io/gh/timbitz/Whippet.jl/branch/master/graph/badge.svg?token=RKE7BSr47v)](https://codecov.io/gh/timbitz/Whippet.jl)
[![Gitter Chat](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/Whippet-jl/Lobby)

## Features
- Splice graph representations of transcriptome structure
  - Build an index for any species with a genome and annotation file
  - _de novo_ event discovery (splicing from/to any combination of annotated donor/acceptor splice sites)
- High speed PolyA+ Spliced Read Alignment (Read lengths <= 255)
  - Repetitive read assignment for gene families
  - On-the-fly alignment/analysis of SRR accession ids using ebi.ac.uk
- Fast and robust quantification of transcriptome structure and expression using EM
  - Dynamic building and entropic measurements of splicing events of any complexity
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
- Differential splicing comparisons
  - Probabilistic calculations of delta PSI leveraging multi-sample biological replicates

Paper pre-print: http://www.biorxiv.org/content/early/2017/07/03/158519

## Get started

### 1) Install 
Install [Julia v0.5.2](https://julialang.org/downloads/oldreleases.html) (NOTE: julia v0.6 is released and will be supported very soon).  If you are new to julia, or installing programs via command line, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started)

### 2) Clone Whippet

```
git clone https://github.com/timbitz/Whippet.jl.git
cd Whippet.jl
julia dependencies.jl
```
Notes: 
* Deprecated syntax `warnings` for dependent packages can be ignored.
* For all executables in Whippet.jl/bin, you can use the `-h` flag to get a list of the available command line options, their usage and defaults.

### 3) Build an index.  
You need your genome sequence in fasta, and a gene annotation file in GTF or Refflat format. Default annotation supplied for hg19 GENCODEv25 TSL1-level transcriptome. 

```bash
$ julia bin/whippet-index.jl --fasta hg19.fa.gz --gtf anno/gencode_hg19.v25.tsl1.gtf.gz
```

Notes: 
* Whippet only uses GTF `exon` lines (others are ignored). These must contain both `gene_id` and `transcript_id` attributes, consistent with the [GTF2.2](http://mblab.wustl.edu/GTF22.html) specification.
* You can specify the output name and location of the index to build using the `--index` parameter, the default (for both whippet-index.jl and whippet-quant.jl) is a generic index named `graph` located at `Whippet.jl/index/graph`.


### 4) Quantify FASTQ files.
```bash
$ julia bin/whippet-quant.jl file.fastq.gz
```

Or if you have paired-end RNA-seq data...
```bash
$ julia bin/whippet-quant.jl fwd_file.fastq.gz rev_file.fastq.gz
```

Or you can provide a link/s to the fastq file/s on the web using `--url`, or just the experiment accession id using `--ebi`!  These will align on the fly as the reads are downloaded from the web directly into alignment. For example:
```
$ julia bin/whippet-quant.jl --ebi SRR1199010
```
is equivalent to
```bash
$ julia bin/whippet-quant.jl --url ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/000/SRR1199010/SRR1199010.fastq.gz
```

You can output the alignments in SAM format with the `--sam` flag and convert to bam with a pipe:
```bash
$ julia bin/whippet-quant.jl fwd_file.fastq.gz --sam | samtools view -bS - > fwd_file.bam
```

It is also possible to pool fastq files at runtime using shell commands, and the optional (`--force-gz`) for pooled gz files (files without .gz suffix)
```bash
$ julia bin/whippet-quant.jl <( cat SRR208080{1,2,3,4,5,6,7,8,9}.fastq.gz ) --force-gz -o SRR208080_1-9
```

### 5) Compare multiple psi files
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

The output format for `whippet-quant.jl` is saved into two core quant files a `.psi.gz` and a `.tpm.gz`.

The `.tpm.gz` file contains a simple format compatible with many downstream tools:

Gene | TpM | Read Counts
---- | --- | -----------
NFIA | 2897.11 | 24657.0

---

Meanwhile the `.psi.gz` file is a bit more complex and requires more explanation:

Gene | Node | Coord | Strand | Type | Psi | CI Width | CI Lo,Hi | Total Reads | Complexity | Entropy | Inc Paths | Exc Paths
---- | ---- | ----- | ------ | ---- | --- | -------- | -------- | ----------- | ---------- | ------- | --------- | ---------
NFIA | 2 | chr1:61547534-61547719 | + | AF | 0.782 | 0.191 | 0.669,0.86 | 49.0 | K1 | 0.756 | 2-4-5:0.782 | 1-5:0.218
NFIA | 4 | chr1:61548433-61548490 | + | CE | 0.8329 | 0.069 | 0.795,0.864 | 318.0 | K2 | 1.25 | 2-4-5:0.3342,3-4-5:0.4987 | 1-5:0.1671
NFIA | 5 | chr1:61553821-61554352 | + | CE | 0.99 | NA | NA | NA | NA | NA | NA | NA
NFIA | 6 | chr1:61743192-61743257 | + | CE | 0.99 | NA | NA | NA | NA | NA | NA | NA

In contrast to many other splicing quantification tools, Whippet allows for dynamic quantification of observed splicing patterns.  Therefore, in order to maintain consistent output from different Whippet runs on various samples, the basic unit of quantification is a SpliceGraph `node`.  It is possible (and even likely) that many nodes are never spliced entirely on their own as is the case with alternative 5' and 3' splice sites, and core exon nodes whose neighboring alt 5' or 3' splice sites are used.  Therefore this must be taken into account when intersecting `.psi.gz` coordinate output with other formats that represent full exons (which can be one or more adjacent nodes combined).


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
 
Each node is defined by a type (above) and has a corresponding value for `Psi` or the Percent-Spliced-In followed by the 90% confidence interval (both the width as well as lower and higher boundaries).

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
  * Use only the highest quality annotations you can find (for human we use TSL1-level). 
  * Avoid giving annotations with 'indels' such as ESTs or mRNAs without filtering out invalid splice sites first.
  * If you plan to align very short reads (~36nt), decrease the Kmer size (we have used 6nt before), otherwise the default Kmer size (9nt) should be used.

---

## Troubleshooting

With all of the executables in Whippet.jl/bin, you can use the `-h` flag to get a list of the available command line options and their usage.  If you are having trouble using or interpreting the output of `Whippet.jl` then please ask a question in our [gitter chat](https://gitter.im/Whippet-jl/Lobby)!.  If you are having trouble running a Whippet executable in /bin, try running `julia test/runtests.jl` in the main `Whippet.jl` directory, Whippet should run cleanly if your environment is working correctly (though it is untested on Windows, it should still work).  If you still think you have found a bug feel free to open an issue in github or make a pull request! 
