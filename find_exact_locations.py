#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 10:37:37 2022

    run python find_exact_locations.py -i input.txt [-o output.txt]
or
    run python find_exact_locations.py --input input.txt [--output output.txt]
    
@author: akmami
"""

import sys
import argparse
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output', default='output.txt')
args = parser.parse_args()

soft_clips = args.input
output = open(args.output, 'w')

if soft_clips is None:
    print("Please enter soft clips file file name.")
    exit(-1)

sv_pos_list = []
sv_end_list = []


with open(soft_clips, 'r') as f:
    not_valid_interval_sv_count = 0
    for line in f:
        if '#' in line[0]:
            
            if len(sv_pos_list) == 0 and len(sv_end_list) == 0:
                del sv_pos_list[:]
                del sv_end_list[:]
                output.write(line)
                continue
            
            # count the locations and find most frquent one
            result_pos = Counter(sv_pos_list).most_common(1)
            begin = -1
            if len(result_pos) > 0 and result_pos[0][1] > 1:
                begin = result_pos[0][0]
            
            result_end = Counter(sv_end_list).most_common(1)
            end = -1
            if len(result_end) > 0 and result_end[0][1] > 1:
                end = result_end[0][0]
            
            output.write('\t'.join([sv_id, sv_chrom, str(begin), str(end)]) + '\n')  # write to file
            
            del sv_pos_list[:]
            del sv_end_list[:]
            output.write(line)
            continue

        line = line.split('\t')
        
        sv_id = line[0]
        sv_chrom = line[1]
        sv_loc = line[2]
        sv_soft_place = line[3]
        
        if '-1' in sv_soft_place:
            sv_pos_list.append(sv_loc)
        else:
            sv_end_list.append(sv_loc)
            
    f.close()
    
output.close()
