# MIT License
#
# Copyright (c) 2023 Akmuhammet Ashyralyyev
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import argparse
import random

__MİN_SV_LENGHT = 50

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str)
parser.add_argument("-c", "--chr", type=str, default="")
parser.add_argument("-l", "--length", type=int, default=__MİN_SV_LENGHT)
parser.add_argument("-o", "--output", type=str)
parser.add_argument("--tag", type=str)
parser.add_argument("--DEL", type=str)
parser.add_argument("--INS", type=str)
parser.add_argument("--INV", type=str)
args = parser.parse_args()

arg_vcf = args.input
arg_chrom = args.chr
arg_min_sv_length = args.length
arg_out = args.output
arg_tag = args.tag
arg_del = args.DEL
arg_ins = args.INS
arg_inv = args.INV


id_index = 1

if arg_vcf is None:
    print("Please enter vcf file name.")
    exit(-1)

if len(arg_vcf) < 4 or not arg_vcf[-4:] == ".vcf":
    print("Please enter valid vcf file.")
    exit(-1)

if arg_out is None:
    arg_out = arg_vcf[:-4] + ".sim.vcf"
    print("Output file is set '{}' as defualt".format(arg_out))

if arg_tag == "None":
    arg_tag = None

if arg_tag is not None and arg_del is None:
    arg_del = "DEL"
    print("Default tag for DELETION is set to '{}'".format(arg_del))

if arg_tag is not None and arg_ins is None:
    arg_ins = "INS"
    print("Default tag for INSERTION is set to '{}'".format(arg_ins))

if arg_tag is not None and arg_inv is None:
    arg_inv = "INV"
    print("Default tag for INVERSION is set to '{}'".format(arg_inv))

if arg_tag is not None:
    
    with open(arg_vcf, "r") as f:
        info_tag = "##INFO=<ID=" + arg_tag

        for line in f:
            if line.startswith(info_tag):
                break
            if line.startswith("#CHROM"):
                print("Tag name '{}' is not fount in meta field. Please provide valid tag name ofr SV type.".format(arg_tag))
                exit(-1)

out_f = open(arg_out, "w")

with open(arg_vcf, "r") as f:
    intro = True
    description = False

    for line in f:
        if intro:
            if line.startswith("##INFO"):
                if not description:
                    out_f.write('##INFO=<ID=SVELDT,Number=1,Type=String,Description="The SV is tagged by SVELDT program:SIMULATED=The SV is only simulated var varsim.py and not processed by sveldt yet, SUCCESS=SVELDT was able to refine all given intervals, PARTIAL=SVELDT was able to refine only one of the points, INCORRECT=SVELDT detected invalid SV."\n')
                    description = True
            if line.startswith("#CHROM"):
                if not description:
                    out_f.write('##INFO=<ID=SVELDT,Number=1,Type=String,Description="The SV is tagged by SVELDT program:SIMULATED=The SV is only simulated var varsim.py and not processed by sveldt yet, SUCCESS=SVELDT was able to refine all given intervals, PARTIAL=SVELDT was able to refine only one of the points, INCORRECT=SVELDT detected invalid SV."\n')
                    description = True
                intro = False
            
            out_f.write(line)
            continue
        
        splitted_line = line.split("\t")
                
        # check validity of chromosome
        if splitted_line[0].startswith("chr"):
            splitted_line[0] = splitted_line[0][3:]
        
        if arg_chrom != "" and splitted_line[0] != arg_chrom:
            continue
        
        # find sv type
        sv_type = ""
        sv_info_tag = "Invalid"
        sv_len = -1

        if arg_tag is not None:
            if arg_tag + "=" + arg_del in splitted_line[7]:
                sv_info_tag = arg_tag + "=" + arg_del
                sv_type = "DEL"
            elif arg_tag + "=" + arg_ins in splitted_line[7]:
                sv_info_tag = arg_tag + "=" + arg_ins
                sv_type = "INS"
            elif arg_tag + "=" + arg_inv in splitted_line[7]:
                sv_info_tag = arg_tag + "=" + arg_inv
                sv_type = "INV"
            else:
                # Leave other variation in vcf.
                out_f.write(line)
                continue
        else:
            if len(splitted_line[3]) > len(splitted_line[4]):
                sv_type = "DEL"
            elif len(splitted_line[3]) < len(splitted_line[4]):
                sv_type = "INS"
            else:
                # probably it is not dup either
                # mismathc, dup, etc.
                out_f.write(line)
                continue

        # set end
        end = str( int(splitted_line[1]) + 1 )

        if sv_type == "DEL":
            end = str( int(splitted_line[1]) + len(splitted_line[3]) - len(splitted_line[4]) + 1 )

        # find sv length
        if sv_type == "DEL" or sv_type == "INS":
            sv_len = len(splitted_line[4]) - len(splitted_line[3])
        
        # filter the sv according to its length
        if sv_len < arg_min_sv_length and sv_len > -arg_min_sv_length:
            out_f.write(line)
            continue

        # set confidence intervals
        outer_start = - abs( int( random.random() * sv_len * 0.06 + sv_len * 0.01 ) ) - 25
        inner_start = abs( int( random.random() * sv_len * 0.06 + sv_len * 0.01 ) ) + 25

        splitted_line[3] = splitted_line[3][0]
        splitted_line[4] = splitted_line[4][0]
       
        splitted_line[7] += ";CIPOS={},{}".format(outer_start, inner_start)
        
        if sv_type != "INS":
            inner_end = - abs( int( random.random() * sv_len * 0.06 + sv_len * 0.01 ) ) - 25
            outer_end = abs( int( random.random() * sv_len * 0.06 + sv_len * 0.01 ) ) + 25
            splitted_line[7] += ";CIEND={},{}".format(inner_end, outer_end)

        splitted_line[7] += ";END={}".format(end)
        splitted_line[7] += ";SVELDT=SIMULATED"
        
        if sv_info_tag in splitted_line[7]:
            splitted_line[7] = splitted_line[7].replace(sv_info_tag, "SVTYPE={}".format(sv_type))
        else:
            splitted_line[7] = splitted_line[7] + ";SVTYPE={}".format(sv_type)

        if splitted_line[2] == ".":
            splitted_line[2] = "GoldStandard{}".format(id_index)
            id_index += 1
        
        out_f.write("\t".join(splitted_line))
        
    print("Simualtion of vcf is successful.")

    out_f.close()
