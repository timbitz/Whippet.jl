# Whippet
##### Ultra fast & lightweight quantification of gene expression and event-specific splicing levels from RNA-seq.

### Why use Whippet?

### Features
- High performance PolyA+ Spliced Read Alignment
- Robust quantification of the expression and transcriptome structure of model and non-model organisms
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
- Accurate splice graph representations of high complexity event types (splicing and alt-3'/5' end usage)
  - Pseudo _de novo_ event discovery
  - Circular splicing discovery
- Accurate repetitive read assignment for gene families

### How to use Whippet

## 1) Install 
Install most recent [julia release here](http://julialang.org/downloads/), which must be v0.4.x! (v0.5 will be supported in future releases).  If you are new to julia, or installing programs via command line, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started)

## 2) Clone Whippet
Make sure dependencies are satisfied. Executables are in bin/
```
git clone https://github.com/timbitz/Whippet.git
cd Whippet
julia dependencies.jl
```
NOTE: `julia dependencies.jl` may be noisy with deprecated syntax warnings.  This is due to the rapid pace at which base julia is being developed and does not actually mean that there was/is a fatal problem with Whippet or its dependencies.

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

## 5) Compare multiple psi files
```bash
$ ls *.psi.gz
sample1-r1.psi.gz sample1-r2.psi.gz sample2-r1.psi.gz sample2-r2.psi.gz
$ julia whippet-delta.jl -a sample1 -b sample2
OR
$ julia whippet-delta.jl -a sample1-r1.psi.gz,sample1-r2.psi.gz -b sample2-r1.psi.gz,sample2-r2.psi.gz
```

---

### Output Format


### Advanced Index Building

If you are building an index for a non-model organism or an index for a custom purpose, there are some general guidelines that can help to ensure that the index you build is as effectively as it can be.  For example, a whippet index that is missing many annotated splice sites that are frequently used, may not be able to align all reads well. Similarly, a whippet index that is built with a huge number of decoy splice sites from indels in EST or mRNA annotations, may have too many `splice sites`, which will divide exons into many tiny nodes, making seeding to those segments more difficult. In general you should seek to:
  * Increase the number of true splice sites that `whippet index` is given. 
  * Avoid giving annotations with 'indels' such as ESTs or mRNAs without filtering for valid splice sites first.
  * Choose a Kmer size appropriate for the node size and number in the transcriptome and the read length you plan to align.

```bash
$ bin/whippet-index.jl -h
Whippet v0.1-rc4 loading and compiling... 
usage: whippet-index.jl [-k KMER] --fasta FASTA --flat FLAT
                        [--index INDEX] [-h]

optional arguments:
  -k, --kmer KMER  Kmer size to use for exon-exon junctions (default
                   9) (type: Int64, default: 9)
  --fasta FASTA    File containg the genome in fasta, one entry per
                   chromosome [.gz]
  --flat FLAT      Gene annotation file in RefFlat format
  --index INDEX    Output prefix for saving index 'dir/prefix'
                   (default Whippet/index/graph) (default:
                   "/Users/timsw/Documents/git/Whippet/src/../index/graph")
  -h, --help       show this help message and exit
```

### Custom Alignment Parameters

```bash

$ julia bin/whippet-quant.jl -h
Whippet v0.1-rc4 loading and compiling... 
usage: whippet-quant.jl [-x INDEX] [-o OUT] [-s] [-L SEED-LEN]
                        [-M SEED-TRY] [-T SEED-TOL] [-B SEED-BUF]
                        [-I SEED-INC] [-P PAIR-RANGE] [-X MISMATCHES]
                        [-S SCORE-MIN] [-j] [--stranded] [--rev-pair]
                        [--no-circ] [--no-tpm] [--force-gz] [-h]
                        filename.fastq[.gz] [paired_mate.fastq[.gz]]

positional arguments:
  filename.fastq[.gz]
  paired_mate.fastq[.gz]


optional arguments:
  -x, --index INDEX     Output prefix for saving index 'dir/prefix'
                        (default Whippet/index/graph) (default:
                        "/Users/timsw/Documents/git/Whippet/index/graph")
  -o, --out OUT         Where should the gzipped output go
                        'dir/prefix'? (default:
                        "/Users/timsw/Documents/git/Whippet/output")
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
                        (type: Int64, default: 2)
  -S, --score-min SCORE-MIN
                        Minimum alignment score (matches - mismatches)
                        (type: Int64, default: 45)
  -j, --junc-only       Only use junction reads, no internal exon
                        reads will be considered.
  --stranded            Is the data strand specific? If so, increase
                        speed with this flag
  --rev-pair            Is the second mate the reverse complement of
                        the first? If so, increase speed with this
                        flag
  --no-circ             Do not allow back/circular splicing
  --no-tpm              Should tpm file be sent to
                        output/prefix.tpm.gz? (default on)
  --force-gz            Regardless of suffix, consider read input as
                        gzipped
  -h, --help            show this help message and exit

```

You can alter alignment parameters to match your needs based on read length and average quality etc.
For example: In order to align very short reads (lets say 35nt reads), the default --score-min/-S is 45 so nothing will match. So we can change this to 30 and adjust the seed increment to something more reasonable (10 from default of 18)
```bash
$ julia whippet-quant.jl SRR2080801.fastq.gz -x ../genomes/Strongylocentrotus/Strongylocentrotus -S 25 -I 10 -o SRR2080801
```


