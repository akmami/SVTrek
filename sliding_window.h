#ifndef SLIDING_WINDOW_H
#define SLIDING_WINDOW_H

#include "params.h" 
#include "refinement.h" 

int sliding_window_ins(int chrom, interval inter, uint32_t imprecise_pos, t_arg *params, int windowSize, int slideSize);

#endif // SLIDING_WINDOW_H
