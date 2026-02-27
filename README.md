## Overview of Kente
<img src="kente_logo_v1.png" width="250">

### What is Kente? 
Kente is a pangenome graph–based framework for detecting horizontal gene transfer (HGT) from assembled metagenomic contigs.  

#### How does Kente find HGT events ? 
1. Kente takes as input assembled contigs of varying lengths. But preferably 20kbp and longer for more HGT calls.
2. Kente employs kraken2 as a pre-filtering and profiling step to give us an idea of what is in the sample know whihc graphs to map to. 
3. After the kraken2 step, Kente then uses minigraph to map the contigs against a pre-built bacterial genus graph database. This produces a contig-to-graph alignment file that contains the best matches of regions of the contigs to the graphs

## Kente's Pre-Built Bacterial Pangenome Graph Database.
In order to use Kente, you should have the pangenome graph DB downloaded, here... 
This pangenome graph database was built using minigraph and thousands of NCBI genomes, grouped by genus. 

## Why “Kente”?
Kente is inspired by the traditional Ghanaian woven fabric characterized by vibrant mosaic patterns formed from interlaced strands. No single thread defines the cloth and meaning emerges from the interplay of many strands woven together.
Similarly, bacterial genomes are not purely vertical lineages. Horizontal gene transfer creates mosaic architectures in which genomic segments move across clades, producing structured yet non-linear patterns.
Just as Kente fabric weaves diverse strands into coherent designs, Kente the bioinformatics tool identifies structural lineage switches woven into microbial genomes.

## Installation 
1. The Kente Pangenome Graph Database is available for download ...
## Software Requirements 
1. Python3.10+ 
2. Kraken2
3. Minigraph

## Running Kente 
 

