#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 12:41:43 2022

@author: akmami
"""
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output', default='output.txt')
args = parser.parse_args()

vcf_f = args.input
txt_f = open(args.output, 'w')

if vcf_f is None:
    print("Please enter vcf file name.")
    exit(-1)


# vcf_f = '/Users/akmami/Desktop/AJtrio_TARDIS.bilkentuniv.072319.vcf'
# txt_f = open('/Users/akmami/Desktop/HG002_sv_summary.txt', 'w')

with open(vcf_f, 'r') as f:
    invalid_SV_type = 0
    intro = True
    for line in f:
        if intro and '#CHROM' not in line:
            continue
        elif intro:
            intro = False
            txt_f.write('\t'.join(['ID', 'POS', 'END', 'ALT', 'CIPOS1', 'CIPOS2', 'CIEND1', 'CIEND2']) + '\n')
            continue
        if 'IMPRECISE' not in line:
            continue
        line = line.split('\t')
        pos = line[1]
        sv_id = line[2]
        alt = ''
        if 'INS' in line[4]:
            alt = 'INS'
        elif 'DEL' in line[4]:
            alt = 'DEL'
        elif 'DUP' in line[4]:
            alt = 'DUP'
        elif 'INV' in line[4]:
            alt = 'INV'
        else:
            invalid_SV_type += 1
            continue
        end = line[7].find(';', 0, len(line[7]))
        end = line[7][4:end]
        ciend = line[7].find('CIEND=', 0, len(line[7]))
        ciend = line[7][ciend+6:line[7].find(';', ciend, len(line[7]))]
        cipos = line[7].find('CIPOS=', 0, len(line[7]))
        cipos = line[7][cipos+6:line[7].find(';', cipos, len(line[7]))]
        ciend = ciend.split(',')
        cipos = cipos.split(',')
        
        txt_f.write('\t'.join([sv_id, pos, end, alt, '\t'.join(cipos), '\t'.join(ciend)]) + '\n')
    print('Invalid number of SV', invalid_SV_type) 
txt_f.close()
txt_f.close()
