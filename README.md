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
2. Numpy
3. Kraken2
4. Minigraph

## Running Kente 
Kente operates as 4 step pipeline.
### Step 1: Taxonomic Pre-filtering
First, we use Kraken2 to identify the genera present in your assembled contigs. This drastically reduces the graph alignment search space and runtime.
``
./kente_kraken2.sh <INPUT_DIR> <OUTPUT_DIR> <KRAKEN_DB> <THREADS>
``
```
./kente_kraken2.sh /path/to/folder_with_contigs /path/to/output_folder /path/to/kraken2_database 16
```
 
### Step 2: Graph Alignment to generate GAF files
Second, run minigraph after you have the top genera the contigs aligns to to get the GAF file which contains the paths your contigs traverse across the graphs of interest. 
``
./kente_minigraph_alignment.sh <INPUT_DIR> <OUTPUT_DIR> <GRAPH_DB_DIR> <THREADS> <CONCURRENT_JOBS>
``
```
./kente_minigraph_alignment.sh /path/to/folder_with_contigs /path/to/your/output_folder /path/to/pangenome_graphs 16 4
```
### Step 3: Genus Blocks and Clade Extraction
Parse the resulting .gaf alignment files and assign a specific best matching genus to every genomic block of some bp (for every 500bp default)
``
python3 Kente_extract_genusblocks.py --gaf_dir <INPUT_GAF_DIR> --out <OUTPUT_FILE> --win <WINDOW_SIZE>
``
```
python3 Kente_extract_genusblocks.py \
  --gaf_dir /path/to/output_folder/gaf_by_genus \
  --out /path/to/results/sample_genus_blocks.tsv \
  --win 500
```
### Step 4: Finding HGT (ABA) Patterns 
Finally, scan the extracted blocks for topological transitions. Kente identifies high-confidence A-B-A events, as well as tip and mosaic transfers. 
```
python3 Kente_findHGT.py \
  --input /path/to/results/sample_genus_blocks.tsv \
  --out_aba /path/to/results/ABA_events.tsv \
  --out_other /path/to/results/partial_events.tsv \
  --min_insert 1000 \
  --max_unlab_bp 2000 \
  --min_flank 1000
```


