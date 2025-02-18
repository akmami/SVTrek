#ifndef __REFINEMENT_H__
#define __REFINEMENT_H__

#include "params.h"
#include <stdlib.h>

#define abs(a) (a < 0 ? -a : a)

void deletion(int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t inner_end, uint32_t outer_end, uint32_t imprecise_pos, uint32_t imprecise_end, t_arg *params, uint32_t *res_start, uint32_t *res_end);

void insertion(int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t imprecise_pos, t_arg *params, uint32_t *res_start);

void inversion(int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t inner_end, uint32_t outer_end, uint32_t imprecise_pos, uint32_t imprecise_end, t_arg *params, uint32_t *res_start, uint32_t *res_end);

#endif // __REFINEMENT_H__
