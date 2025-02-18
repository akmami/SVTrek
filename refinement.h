#ifndef __REFINEMENT_H__
#define __REFINEMENT_H__

#include "params.h"

void deletion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, t_arg *params);

void insertion(int chrom, int outer_start, int inner_start, int imprecise_pos, t_arg *params);

void inversion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, t_arg *params);

#endif // __REFINEMENT_H__
