#ifndef __REFINEMENT_H__
#define __REFINEMENT_H__

#include "params.h"
#include <algorithm>      // for std::sort
#include <limits>         // for INT_MAX
#include <cmath>          // for abs
#include <iostream>
#include <vector>

result deletion(int chrom,
                int outer_start,
                int inner_start,
                int inner_end,
                int outer_end,
                int imprecise_pos,
                int imprecise_end,
                struct params &_params,
                struct thread_data &_thread_data);

result insertion(int chrom,
                 int outer_start,
                 int inner_start,
                 int imprecise_pos,
                 struct params &_params,
                 struct thread_data &_thread_data);

result inversion(int chrom,
                 int outer_start,
                 int inner_start,
                 int inner_end,
                 int outer_end,
                 int imprecise_pos,
                 int imprecise_end,
                 struct params &_params,
                 struct thread_data &_thread_data);

#endif // __REFINEMENT_H__
