#!/bin/bash
cd db/pVOGs
mkdir hmmer_db
cat HMMs/*.hmm > hmmer_db/hmmer_db
hmmpress hmmer_db/hmmer_db