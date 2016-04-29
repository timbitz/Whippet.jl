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

1) Install most recent [julia release here](http://julialang.org/downloads/), which must be >= v0.4.  If you are new to julia, or installing programs via command line, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started)

2) Clone whippet and make sure dependencies are satisfied. Executables are in bin/
```
git clone https://github.com/timbitz/Whippet.git
cd Whippet
julia dependencies.jl
```

3) Build an index.  You need your genome sequence in fasta, and a gene annotation file in refflat. A default example is supplied for hg19 in anno/refseq_hg19.flat.gz
```
julia whippet-index --fasta hg19.fa.gz --flat refseq_hg19.flat.gz
```

4) Quantify FASTQ files.
```
julia whippet-quant file.fastq.gz
```
or if you have paired-end RNA-seq data...
```
julia whippet-quant fwd_file.fastq.gz rev_file.fastq.gz
```
