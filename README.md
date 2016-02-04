# Whippet
##### Uber fast & lightweight quantification of gene expression and event-specific splicing levels from RNA-seq.

### Why use Whippet?
With the birth of several new high-performance RNA-seq tools (sailfish/salmon/kallisto), it has become clear that you can achieve accurate quantification of gene-expression without trading high-performance.  However what these programs lack is a means to provide *event-specific* quantification of alternative splicing and 5'/3' end usage.  

### Features
- Ultra-high performance
- Robust quantification in any species
  - Event-specific Percent-spliced-in (PSI)
  - Gene expression (TpM)
  - Single-cell friendly (low depth is OK unlike 'junction-only' tools)
- Differential splicing statistics based on experimental replicates
- Pseudo de novo event discovery
  - Any combination of annotated splice sites is fair game
  - Circular splicing discovery
  - Trans-splicing discovery
- Repetitive read assignment through Expectation Maximization

