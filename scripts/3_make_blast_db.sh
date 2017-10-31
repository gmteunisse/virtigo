#!/usr/bin/bash

cd db/pVOGs/
mkdir blast_db

cd sequences
cat VOG* > all_sequences.fasta
cd ..
makeblastdb -in sequences/all_sequences.fasta -title pVOGs -out blast_db/pVOGs_db -dbtype prot
