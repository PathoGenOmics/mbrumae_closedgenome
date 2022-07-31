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

def sequences_positive(fasta_file:str, start: int, end: int, gene: str, out_gene):
    '''
    To extract gene from fasta file
    '''
    
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        sequence = Seq(str(fasta.seq)[start-1:end]).translate()
        header = str(fasta.id)
        to_write = '>' + header + '\n' + str(sequence) + '\n'
        out_gene.write(to_write)

def sequences_negative(fasta_file:str, start: int, end: int, gene: str, out_gene):
    '''
    to extract gene from fasta file in negative strand
    '''

    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        sequence = Seq(str(fasta.seq)[start-1:end]).reverse_complement().translate()
        header = str(fasta.id)
        to_write = '>' + header + '\n' + str(sequence) + '\n'
        out_gene.write(to_write)

def write_threshold(jobname:str,number:int,lista_genes:list):
    with open(jobname+'.'+ str(number)+'.gene.blastout','w') as out_write:
        with open(jobname+'.blastout','r') as out_results2:
            
            for line in out_results2:
                if '#' not in line:
                    l2 = line.strip('\n')
                    line = line.strip('\n').split('\t')
                    name = line[0]
                    aa_align =int(line[4])
                    aa_miss =int(line[5])
                    aa = aa_align - aa_miss
                    
                    for gene in lista_genes:
                        if gene[0] == name:
                            gene_len = int(gene[4])
                            n_gene = gene[0]
                    
                    per = (aa/gene_len)*100
                    print(name,aa_align,aa_miss,gene_len,n_gene,str(per))
                    if per >= number:
                        sentencia = l2+'\t'+str(per)+'\n'
                        out_write.write(sentencia)

def main():
    parser = argparse.ArgumentParser(description = 'Script to blast proteins aa into a fasta nt')
    parser.add_argument('-g', dest = 'genes', required = True, help = 'file with gene coordinates')
    parser.add_argument('-f', dest = 'fasta', required = True, help = 'file fasta to blast')
    parser.add_argument('-p', dest = 'protein', required = True, help = 'file fasta to get proteins')
    parser.add_argument('-n', dest = 'jobname', required = True, help = 'name for job')
    args = parser.parse_args()
    lista_genes = read_genes(args.genes) # list of genes to analyze


    for gene in lista_genes:
        name_gene = str(gene[0])
        or_gene = str(gene[1])
        start_gene = int(gene[2])
        end_gene = int(gene[3])
        with open(args.jobname+'.multi_aa.fasta', 'a') as in_file:
            if or_gene == '+':
                sequences_positive(args.fasta, start_gene, end_gene, name_gene, in_file)
            elif or_gene == '-':
                sequences_negative(args.fasta, start_gene, end_gene, name_gene, in_file)
            
    blast_command ='tblastn -query ' + args.jobname + '.multi_aa.fasta -subject ' + args.protein + ' -outfmt 7 -max_target_seqs 2 -evalue 1e-6 -max_hsps 1 -out ' + args.jobname + '.blast.out'
    os.system(blast_command)



if __name__ == '__main__':
    main()