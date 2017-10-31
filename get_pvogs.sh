#!/bin/bash

#This script needs to be run before being able to run Virtigo

echo Downloading pVOGs database...
python scripts/1_get_pVOGs.py
python scripts/2_get_sequences.py
echo Creating BLAST database...
bash scripts/3_make_blast_db.sh
echo Creating HMMer database...
bash scripts/4_create_hmmer_db.sh
echo Finished creating database.
