#!/usr/bin/python3
# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.Seq import Seq
import os
import argparse

#command line python3 blast_proteins.py -g gene.txt -f h37rv.fasta -p h37rv.fasta -n prueba
#module load biotools
#module load python/3.8

def read_genes(genes_file:str):
    '''
    To store genes info of a file into a list
    '''
    list_genes = list()
    with open(genes_file,'r') as in_file:
        for line in in_file:
            if '#' not in line:
                line = line.strip('\n').split('\t')
                list_genes.append(line)
    return list_genes

def sequences_positive(fasta_file:str, start: int, end: int, gene: str):
    '''
    To extract gene from fasta file
    '''
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        sequence = Seq(str(fasta.seq)[start-1:end])
        header = str(fasta.id)
        to_write = '>' + gene + '\n' + str(sequence) + '\n'
        return to_write

def sequences_negative(fasta_file:str, start: int, end: int, gene: str):
    '''
    to extract gene from fasta file in negative strand
    '''
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        sequence = Seq(str(fasta.seq)[start-1:end]).reverse_complement()
        header = str(fasta.id)
        to_write = '>' + gene + '\n' + str(sequence) + '\n'
        return to_write

def main():
    parser = argparse.ArgumentParser(description = 'script to blast proteins into a fasta nt')
    parser.add_argument('-g', dest = 'genes', required = True, help = 'file with gene coordinates')
    parser.add_argument('-f', dest = 'fasta', required = True, help = 'file fasta to blast')
    parser.add_argument('-n', dest = 'jobname', required = True, help = 'name for job')
    args = parser.parse_args()

    lista_genes = read_genes(args.genes)
    with open(args.jobname+'.genes.fa','w') as out_results:

        for gene in lista_genes:
            name_gene = str(gene[0])
            or_gene = str(gene[1])
            start_gene = int(gene[2])
            end_gene = int(gene[3])
            if or_gene == '+':
                sentencia = sequences_positive(args.fasta, start_gene, end_gene, name_gene)
            elif or_gene == '-':
                sentencia = sequences_negative(args.fasta, start_gene, end_gene, name_gene)
            out_results.write(sentencia)
if __name__ == '__main__':
    main()