# K-meriser 

Takes an input reference genome, refernce index, k, and number of threads and generates for each chromosome a tsv file containing the unique overlapping k-mers in it, their position relative to the chromsome, their position relative to the genome start, and their number of ocurrences. 

TODO: add a step in the end to merge the results from the different chromosomes. 

# Prerequisits
* Download and install rust - the easiest way is using brew: `brew install rust`
* Download the relevant reference genome 
* Download and install samtools - `brew install samtools`
* Index the reference genome `samtools faidx <path/to/downloaded/genome>`
  
# Setup
* To build it run: `cargo build`
* To run it execute: `RUST_LOG=info cargo run -- <path to reference genome> <path to reference genome index> <k, e.g. 25> <number of threads to use> <path/to/output/`

For example: `RUST_LOG=info cargo run -- /Users/zoharetzioni/Downloads/GCF_000001405.40_GRCh38.p14_genomic.fna /Users/zoharetzioni/Downloads/GCF_000001405.40_GRCh38.p14_genomic.fna.fai 25 8 /Users/zoharetzioni/Downloads/genome-k-meriser/hg38/` 

