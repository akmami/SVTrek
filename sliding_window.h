#ifndef __SLIDING_WINDOW_H__
#define __SLIDING_WINDOW_H__

#include "params.h" 
#include "utils.h"
#include "refinement.h" 
#include <stdlib.h>
#include <stdio.h>
#include <htslib/sam.h>

int refine_ins_disc(int chrom, interval inter, t_arg *params, int window_size, int slide_size);

#endif // SLIDING_WINDOW_H
