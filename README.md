# Whippet
##### Ultra fast & lightweight quantification of gene expression and event-specific splicing levels from RNA-seq.

### Why use Whippet?

### Features
- High performance
- Robust quantification in any species
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
- Accurate splice graph representations of high complexity event types (splicing and alt-3'/5' end usage)
  - Pseudo _de novo_ event discovery
  - Circular splicing discovery
- Accurate repetitive read assignment for gene families

### How to use Whippet

1) Install most recent julia release [here](http://julialang.org/downloads/), must be >= v0.4

2) Clone whippet and make sure dependencies are satisfied. Executables are in bin/
```
git clone https://github.com/timbitz/Whippet.git
cd Whippet
julia dependencies.jl
```

3) Build an index.  You need your genome sequence in fasta, and a gene annotation file in refflat.
```
whippet-index --fasta genome.fa.gz --flat gene_model.flat.gz
```

4) Quantify FASTQ files.
```
whippet-quant my-file.fastq.gz
```
