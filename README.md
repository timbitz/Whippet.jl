# Whippet
##### Ultra fast & lightweight quantification of gene expression and event-specific splicing levels from RNA-seq.

### Why use Whippet?

### Features
- High performance
- Robust quantification of the expression and transcriptome structure of model and non-model organisms
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
- Accurate splice graph representations of high complexity event types (splicing and alt-3'/5' end usage)
  - Pseudo _de novo_ event discovery
  - Circular splicing discovery
- Accurate repetitive read assignment for gene families

### How to use Whippet

## 1) Install 
Install most recent [julia release here](http://julialang.org/downloads/), which must be >= v0.4.  If you are new to julia, or installing programs via command line, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started)

## 2) Clone Whippet
Make sure dependencies are satisfied. Executables are in bin/
```
git clone https://github.com/timbitz/Whippet.git
cd Whippet
julia dependencies.jl
```

## 3) Build an index.  
You need your genome sequence in fasta, and a gene annotation file in refflat. A default example is supplied for hg19 in anno/refseq_hg19.flat.gz
```bash
$ julia whippet-index.jl --fasta hg19.fa.gz --flat refseq_hg19.flat.gz
```

## 4) Quantify FASTQ files.
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

It is also possible to pool fastq files at runtime using shell commands, and the optional (`--force-gz`) for pooled gz files
```bash
$ julia whippet-quant.jl <( cat SRR208080{1,2,3,4,5,6,7,8,9}.fastq.gz ) --force-gz -o SRR208080_1-9
```

You can alter alignment parameters to match your needs based on read length and average quality etc.
For example: In order to align very short reads (lets say 35nt reads), the default --score-min/-S is 45 so nothing will match. So we can change this to 30 and adjust the seed increment to something more reasonable (10 from default of 18)
```bash
$ julia whippet-quant.jl SRR2080801.fastq.gz -x ../genomes/Strongylocentrotus/Strongylocentrotus -S 25 -I 10 -o SRR2080801
```

## 5) Compare multiple psi files
```bash
$ ls *.psi.gz
sample1-r1.psi.gz sample1-r2.psi.gz sample2-r1.psi.gz sample2-r2.psi.gz
$ julia whippet-delta.jl -a sample1 -b sample2
OR
$ julia whippet-delta.jl -a sample1-r1.psi.gz,sample1-r2.psi.gz -b sample2-r1.psi.gz,sample2-r2.psi.gz
```
