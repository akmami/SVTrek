#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:20:07 2022

@author: akmami

    run python extract_simulated_data.py -c chromosom -i input.vcf [-o output.txt]
or
    run python extract_simulated_data.py --chr chromosom --input input.vcf [--output output.txt]

Default output file is named output.txt

"""

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input')
parser.add_argument('-c', '--chr', default='')
parser.add_argument('-o', '--output', default='output.txt')
args = parser.parse_args()

chr = args.chr

vcf_f = args.input
txt_f = open(args.output, 'w')

if vcf_f is None:
    print("Please enter vcf file name.")
    exit(-1)

with open(vcf_f, 'r') as f:
    intro = True
    invalid_len_count = 0
    for line in f:
        if intro and '#CHROM' not in line:
            continue
        elif intro:
            intro = False
            txt_f.write('\t'.join(['ID', 'CHROM', 'POS', 'END', 'CIPOS1', 'CIPOS2', 'CIEND1', 'CIEND2']) + '\n')
            continue
        if 'SVLEN' not in line:
            continue
        line = line.split('\t')
        pos = line[1]
        chrom = line[0]
        
        # check validity of chromosome
        if not chrom.isdigit() and 'chr' in chrom:
        	chrom = chrom[3:]
        
        if chr != '' and chrom != chr:
            continue
        
        # find sv length
        sv_id = line[2]
        svlen = line[7].find('SVLEN=', 0, len(line[7]))
        svlen_end = line[7].find(';', svlen, len(line[7]))
        
        if svlen_end == -1:
        	svlen_end = len(line[7])
        
        svlen = line[7][svlen+6:svlen_end]
        if ',' in svlen:
            continue
        if svlen[0] == '-':
            svlen = int(svlen[1:])
        else:
            svlen = int(svlen)
        
        # filter the sv according to its length
        if svlen < 50:
        	continue
        
        # set end
        end = str( int(pos) + svlen + 1 )
        
        # extract sv type
        sv_type = 'none'
        
        if any('DEL' in s for s in line):
            sv_type = 'del'
        elif any('INS' in s for s in line):
            sv_type = 'ins'
        elif any('INV' in s for s in line):
            sv_type = 'inv'
        elif (not any('DUP' in s for s in line)) and svlen > 0:
            sv_type = 'ins'
        else:
            print("invalid sv type")
            print(line)
            continue
        
        # if sv is insertion, then end = pos + 1
        if sv_type == 'ins':
            end = str(int(pos) + 1)
        
        # set cofidence intervals
        change = 25
        
        if svlen < -500 or svlen > 500:
        	change = 100
        elif svlen < -250 or svlen > 250:
        	change = 50
        
        outer_start = str(0-change)
        inner_start = str(change)
        
        inner_end = str(0-change)
        outer_end = str(change)
        
            
        txt_f.write('\t'.join([sv_id, chrom, pos, end, sv_type, outer_start, inner_start, inner_end, outer_end]) + '\n')
    print("invalid number of svlen lines", invalid_len_count)
txt_f.close()
txt_f.close()
