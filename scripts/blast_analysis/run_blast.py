#!/usr/bin/python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
import os
import argparse

def read_genes(genes_file:str):
    '''
    To store genes info of a file into a list
    '''
    list_genes = list()
    with open(genes_file,'r') as in_file:
        for line in in_file:
            if '#' not in line:
                list_genes.append(line.strip('\n').split('\t'))
    return list_genes

def main():
    parser = argparse.ArgumentParser(description = 'Script to blast proteins aa into a fasta nt')
    parser.add_argument('-g', dest = 'genes', required = True, help = 'file with gene coordinates')
    parser.add_argument('-f', dest = 'fasta', required = True, help = 'file fasta to blast')
    parser.add_argument('-p', dest = 'protein', required = True, help = 'file fasta to get proteins')
    parser.add_argument('-n', dest = 'jobname', required = True, help = 'name for job')
    args = parser.parse_args()
    lista_genes = read_genes(args.genes) # list of genes to analyze
if __name__ == '__main__':
    main()