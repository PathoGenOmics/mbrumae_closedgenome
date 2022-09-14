#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import argparse


def read_file(in_file,name):
    with open(name+'mummer.out','w') as out_w:
        with open(in_file,'r') as in_r:
            for line in in_r:
                if ">" in line:
                    header = line.strip('\n').strip('>')
                    h2 = header.strip(' Reverse')
                    if "Reverse" in header:
                        stat = "Reverse"
                    else:
                        stat = "Forward"
                else:
                    l = line.replace('  ',' ')
                    l = l.replace('  ',' ')
                    l = l.replace('  ',' ')
                    l = l.replace('  ',' ')
                    l = l.replace(' ','\t')
                    l = l.strip('\n').split('\t')
                    print(l)
                    n1 = int(l[2])
                    n2 = int(l[3])
                    if n2 > n1:
                        pair = n2-n1
                    else: 
                        pair = n1-n2
                    escribe = h2 + '\t' + stat + '\t' + str(pair) + '\n'
                    out_w.write(escribe)

def main():
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-f', dest = 'file', required = True, help = 'file with gene mummer results')
    parser.add_argument('-n', dest = 'name', required = True, help = 'name')

    args = parser.parse_args()
    read_file(args.file,args.name)


if __name__ == '__main__':
    main()