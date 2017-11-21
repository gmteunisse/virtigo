# Virtigo.
Virtigo is an annotation pipeline for viral contigs. It is comparable to PROKKA, but specific to viral sequences. The pipeline takes viral contigs, predicts open reading frames (ORFs) and then maps these ORFs to RefSeq Viral / pVOGs to obtain taxonomic and functional annotations. Taxonomic annotations are determined using the LCA(\*) algorithm. The predicted ORF taxonomies are subsequently used to assign taxonomic annotations to contigs, also using LCA(\*).

Other features include:
* Outputting a tree-like representation of the taxonomy of each contig.
* Mapping reads to ORFs and contigs to determine the average coverage per ORF/contig for abundance estimation.

## Dependencies
Virtigo is dependent on the following software packages:
* Python (tested with versions 2.7.10 and 3.6.2)
* BEDTools (tested with version 2.26)
* BLAST+ suite (tested with version 2.6.0)
* Samtools (tested with version 1.5)
* HMMer (tested with version 3.1)
* SciPy (optional; tested with version 0.18.1)

Please add the path to the bins of each dependency to your $PATH variable and make sure that the following bins can be called:
* genomeCoverageBed
* blastp
* samtools
* hmmscan

## Installation
Download all scripts and bins from GitHub in your desired manner and store in a target directory. Once downloaded, place the Virtigo directory in your preferred path and the path to you $PATH variable:
```
export PATH=$PATH:path/to/virtigo/dir
```
Virtigo can then be run using Python version 2 or 3 by executing:
```
virtigo.py <options>
```

## Initialization
Before Virtigo can be used, the appropriate databases need to be downloaded and generated. To do so, make sure you have an active internet connection and run Virtigo in initialization / update mode:
```
virtigo.py -i
```

## Usage
Run Virtigo in its standard mode by executing:
```
virtigo.py -c <path/to/contigs_file.fasta> -o <path/to/output_file>
```
For more usage options and help, use the `-h` flag.


