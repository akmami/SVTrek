#ifndef VARIATIONS_CPP
#define VARIATIONS_CPP

#include "params.h"
#include "refinement.cpp"


namespace sv {

	result deletion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, params &_params) {
	    
		result res;
	    refine_start(chrom, outer_start, inner_start, imprecise_pos, res, _params, DELETION);
	    refine_end(chrom, inner_end, outer_end, imprecise_end, res, _params, DELETION);
	    
	    return res;
	};

	result insertion(int chrom, int outer_start, int inner_start, int imprecise_pos,  params &_params) {

		result res;
	    refine_point(chrom, outer_start, inner_start, imprecise_pos, res, _params, INSERTION);

	    return res;
	};

	result inversion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, params &_params) {

	    result res;
	   	refine_point(chrom, outer_start, inner_start, imprecise_pos, res, _params);
	    refine_point(chrom, inner_end, outer_end, imprecise_end, res, _params);

	    return res;
	};

};

#endif