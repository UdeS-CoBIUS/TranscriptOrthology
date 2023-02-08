#!/bin/bash

python3 ./scripts/transcriptOrthology.py -talg ./execution/inputs/transcripts_alignments/ENSGT00390000000104.alg -galg ./execution/inputs/genes_alignments/ENSGT00390000000104.alg -gtot ./execution/inputs/mapping_gene_to_transcripts/ENSGT00390000000104.fasta -nhxt ./execution/inputs/NHX_trees/ENSGT00390000000104.nwk -lowb 0.7 -outf ./execution/outputs/

