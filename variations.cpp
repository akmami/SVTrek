#ifndef VARIATIONS_CPP
#define VARIATIONS_CPP

#include "params.h"
#include "refinement.cpp"

struct result {
	int refined_start;
	int refined_end;
};

namespace sv {

	result deletion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, params &_params) {
	    
		result res;
	    res.refined_start = refine_start(chrom, outer_start, inner_start, imprecise_pos, _params);
	    res.refined_end = refine_end(chrom, inner_end, outer_end, imprecise_end, _params);
	    
	    return res;
	};

	result insertion(int chrom, int outer_start, int inner_start, int imprecise_pos,  params &_params) {

		result res;
	    int refined_start = refine_point(chrom, outer_start, inner_start, imprecise_pos, _params);

	    return res;
	};

	result inversion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, params &_params) {

	    result res;
	    res.refined_start = refine_point(chrom, outer_start, inner_start, imprecise_pos, _params);
	    res.refined_end = refine_point(chrom, inner_end, outer_end, imprecise_end, _params);

	    return res;
	};

};

#endif