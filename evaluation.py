#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  20 10:37:37 2022

    run python evaluation.py -vcf input.vcf -r result.vcf [-o evaluation.result.txt]
or
    run python evaluation.py --vcf input.vcf -result result.vcf [--output evaluation.result.txt]
    
@author: akmami
"""

import sys
import argparse
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-vcf', '--vcf')
parser.add_argument('-r', '--result')
parser.add_argument('-o', '--output')
args = parser.parse_args()

vcf = args.vcf
result = args.result
out = args.output


if vcf is None:
    print("Please enter vcf file path.")
    exit(-1)


if result is None:
    print("Please enter vcf result file path.")
    exit(-1)

if out is not None:
    output = open(out, "w")
    redirect_stdout(output)


vcf_sv = dict()

with open(vcf, 'r') as f1:
    intro = True
    invalid_len_count = 0
    for line in f1:
        if intro and '#CHROM' not in line:
            continue
        elif intro:
            intro = False
            continue
        if 'SVLEN' not in line:
            continue
        line = line.split('\t')
        pos = int(line[1])
        chrom = line[0]
        
        # check validity of chromosome
        if not chrom.isdigit() and 'chr' in chrom:
        	chrom = chrom[3:]
        
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
        
        end = int(pos) + svlen + 1
        
        if any('DEL' in s for s in line):
            vcf_sv[sv_id] = {"pos": pos, "end": end, "chrom": chrom, "svlen": svlen, "sv_type": "del"}
        elif any('INS' in s for s in line):
            vcf_sv[sv_id] = {"pos": pos, "chrom": chrom, "svlen": svlen, "sv_type": "ins"}
        elif any('INV' in s for s in line):
            vcf_sv[sv_id] = {"pos": pos, "end": end, "chrom": chrom, "svlen": svlen, "sv_type": "inv"}
        elif (not any('DUP' in s for s in line)) and svlen > 0:
            vcf_sv[sv_id] = {"pos": pos, "chrom": chrom, "svlen": svlen, "sv_type": "ins"}
        else:
            continue


ins_pos = dict()
del_pos = dict()
del_end = dict()
inv_pos = dict()
inv_end = dict()

ins_gr = 0
ins_ls = 0
del_pos_gr = 0
del_pos_ls = 0
del_end_gr = 0
del_end_ls = 0
inv_pos_gr = 0
inv_pos_ls = 0
inv_end_gr = 0
inv_end_ls = 0

with open(result, 'r') as f1:

    for line in f1:
        if '#' == line[0]:
            continue
        
        line = line.split('\t')
        chrom = line[0]
        sv_id = line[1]
        sv_type = line[2]
        sv_pos = int(line[3])

        if sv_id not in vcf_sv:
            print(sv_id + " does not exists")
            continue

        sv = vcf_sv[sv_id]

        if sv_type == 'del':
            sv_end = int(line[4])
            if sv_pos == -1:
                del_pos[float('inf')] = del_pos[float('inf')] + 1 if float('inf') in del_pos else 1
            else:
                del_pos[sv["pos"]-sv_pos] = del_pos[sv["pos"]-sv_pos] + 1 if sv["pos"]-sv_pos in del_pos else 1
        
            if sv_end == -1:
                del_end[float('inf')] = del_end[float('inf')] + 1 if float('inf') in del_end else 1
            else:
                del_end[sv["end"]-sv_end] = del_end[sv["end"]-sv_end] + 1 if sv["end"]-sv_end in del_end else 1

            if sv_pos == -1 or abs(sv["pos"]-sv_pos) > 10:
                del_pos_gr += 1
            else:
                del_pos_ls += 1
  
            if sv_end == -1 or abs(sv["end"]-sv_end) > 10:
                del_end_gr += 1
            else:
                del_end_ls += 1
            
        
        if sv_type == 'ins':
            if sv_pos == -1:
                ins_pos[float('inf')] = ins_pos[float('inf')] + 1 if float('inf') in ins_pos else 1
            else:
                ins_pos[sv["pos"]-sv_pos] = ins_pos[sv["pos"]-sv_pos] + 1 if sv["pos"]-sv_pos in ins_pos else 1

            if sv_pos == -1 or abs(sv["pos"]-sv_pos) > 10:
                ins_gr += 1
            else:
                ins_ls += 1

        if sv_type == 'inv':
            sv_end = int(line[4])
            if sv_pos == -1:
                inv_pos[float('inf')] = inv_pos[float('inf')] + 1 if float('inf') in inv_pos else 1
            else:
                inv_pos[sv["pos"]-sv_pos] = inv_pos[sv["pos"]-sv_pos] + 1 if sv["pos"]-sv_pos in inv_pos else 1

            if sv_end == -1:
                inv_end[float('inf')] = inv_end[float('inf')] + 1 if float('inf') in inv_end else 1
            else:
                inv_end[sv["end"]-sv_end] = inv_end[sv["end"]-sv_end] + 1 if sv["end"]-sv_end in inv_end else 1
            
            if sv_pos == -1 or abs(sv["pos"]-sv_pos) > 10:
                inv_pos_gr += 1
            else:
                inv_pos_ls += 1
  
            if sv_end == -1 or abs(sv["end"]-sv_end) > 10:
                inv_end_gr += 1
            else:
                inv_end_ls += 1

print("Results:")

print("Insertion pos:")
for key in sorted(ins_pos.keys()):
    print("{} : {}".format(key, ins_pos[key]))
print("")

print("Deletion pos:")
for key in sorted(del_pos.keys()):
    print("{} : {}".format(key, del_pos[key]))
print("Deletion end:")
for key in sorted(del_end.keys()):
    print("{} : {}".format(key, del_end[key]))
print("")

print("Inversion pos:")
for key in sorted(inv_pos.keys()):
    print("{} : {}".format(key, inv_pos[key]))
print("Inversion end:")
for key in sorted(inv_end.keys()):
    print("{} : {}".format(key, inv_end[key]))
print("")

if out is not None:
    output.close()

print("insertion pos")
print(">10\t<10")
print(str(ins_gr) + "\t" + str(ins_ls))
print("")

print("deletion pos")
print(">10\t<10")
print(str(del_pos_gr) + "\t" + str(del_pos_ls))
print("deletion end")
print(">10\t<10")
print(str(del_end_gr) + "\t" + str(del_end_ls))
print("")

print("inversion pos")
print(">10\t<10")
print(str(inv_pos_gr) + "\t" + str(inv_pos_ls))
print("inversion end")
print(">10\t<10")
print(str(inv_end_gr) + "\t" + str(inv_end_ls))
print("")
