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

def sequences_positive(fasta_file:str, start: int, end: int, gene: str):
    '''
    To extract gene from fasta file
    '''
    with open(gene + '.fasta','w') as out_gene:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for fasta in fasta_sequences:
            sequence = Seq(str(fasta.seq)[start-1:end]).translate()
            header = str(fasta.id)
            to_write = '>' + header + '\n' + str(sequence) + '\n'
            out_gene.write(to_write)

def sequences_negative(fasta_file:str, start: int, end: int, gene: str):
    '''
    to extract gene from fasta file in negative strand
    '''
    with open(gene + '.fasta','w') as out_gene:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for fasta in fasta_sequences:
            sequence = Seq(str(fasta.seq)[start-1:end]).reverse_complement().translate()
            header = str(fasta.id)
            to_write = '>' + header + '\n' + str(sequence) + '\n'
            out_gene.write(to_write)


def main():
    parser = argparse.ArgumentParser(description = 'Script to blast proteins aa into a fasta nt')
    parser.add_argument('-g', dest = 'genes', required = True, help = 'file with gene coordinates')
    parser.add_argument('-f', dest = 'fasta', required = True, help = 'file fasta to blast')
    parser.add_argument('-p', dest = 'protein', required = True, help = 'file fasta to get proteins')
    parser.add_argument('-n', dest = 'jobname', required = True, help = 'name for job')
    args = parser.parse_args()
    lista_genes = read_genes(args.genes) # list of genes to analyze

    with open(args.jobname+'.blastout','w') as out_results:
        header = '#gene\tquery\tsubject\tpercentage identity\talignment length\tmismatches\tgap opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbit score\n'
        out_results.write(header)
        for gene in lista_genes:
            name_gene = str(gene[0])
            or_gene = str(gene[1])
            start_gene = int(gene[2])
            end_gene = int(gene[3])
            if or_gene == '+':
                sequences_positive(args.fasta, start_gene, end_gene, name_gene)
            elif or_gene == '-':
                sequences_negative(args.fasta, start_gene, end_gene, name_gene)
            blast_command ='tblastn -query ' + name_gene + '.fasta -subject '+ args.protein+' -outfmt 7 -max_target_seqs 2 -evalue 1e-6 -max_hsps 1 -out '+name_gene+'.blast.out'
            os.system(blast_command)
            with open(name_gene+'.blast.out','r') as in_blast:
                for line in in_blast:
                    if '#' not in line:
                        line = name_gene+'\t'+line
                        out_results.write(line)
            os.system('rm ' + name_gene + '.blast.out')
            os.system('rm ' + name_gene + '.fasta')
if __name__ == '__main__':
    main()