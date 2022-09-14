#!/bin/bash

python3 clean_genes.py
python3 multifasta.py -g clean_genes.tsv -f brumae_dnaA.fasta -j last_genes
run-mummer3 brumae.fasta last_genes.genes.fa brumae_mummer -s -L 
python3 process_mummer_out.py -f brumae_mummer.out -n result_mummer