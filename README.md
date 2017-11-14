# Virtigo.
An annotation pipeline for viral contigs.

# Usage
Download all scripts and bins from GitHub in your desired manner and store in a target directory. Once downloaded, first run get_pvogs.sh to download the databases by navigating to the directory in the Terminal and executing:
```
bash get_pvogs.sh
```
Virtigo can then be run using Python version 2 or 3 by executing:
```
python virtigo.py -c <contigs.fasta> -o <output_file>
```
Use the `-h` flag to get usage help and options. 

Virtigo requires the user to install the following dependencies before usage:
* Python (tested with versions 2.7.10 and 3.6.2)
* BEDTools (tested with version 2.26)
* BLAST+ suite (tested with version 2.6.0)
* Samtools (tested with version 1.5)
* HMMer (tested with version 3.1)
* SciPy (optional; tested with version 0.18.1)

Please add the path to the bins of each dependency to your $PATH variable and modify bin name as such that they can be called as:
* genomeCoverageBed
* blastp
* samtools
* hmmscan
