#ifndef VARIATIONS_CPP
#define VARIATIONS_CPP

#include "params.h"
#include "refinement.cpp"


namespace sv {

	result deletion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, params &_params) {

		result res;
		res.refined_start = -1;
		res.refined_end = -1;
		if ( ( imprecise_end - imprecise_pos ) < __CI_MIN_FILTER_LENGTH || ( imprecise_end - imprecise_pos ) * _params.ci_max_length > inner_start - outer_start ) {
			refine_start(chrom, outer_start, inner_start, imprecise_pos, res, _params, DELETION);
		}
		if ( ( imprecise_end - imprecise_pos ) < __CI_MIN_FILTER_LENGTH || ( imprecise_end - imprecise_pos ) * _params.ci_max_length > outer_end - inner_end ) {
			refine_end(chrom, inner_end, outer_end, imprecise_end, res, _params, DELETION);
		}

		return res;
	};

	result insertion(int chrom, int outer_start, int inner_start, int imprecise_pos,  params &_params) {

		result res;
		res.refined_start = -1;
		res.svlen = -1;
		refine_point(chrom, outer_start, inner_start, imprecise_pos, res, _params, INSERTION);

		return res;
	};

	result inversion(int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, params &_params) {

		result res;
		res.refined_start = -1;
		res.refined_end = -1;
		if ( ( imprecise_end - imprecise_pos ) < __CI_MIN_FILTER_LENGTH || ( imprecise_end - imprecise_pos ) * _params.ci_max_length > inner_start - outer_start ) {
			refine_point(chrom, outer_start, inner_start, imprecise_pos, res, _params);
		}
		if ( ( imprecise_end - imprecise_pos ) < __CI_MIN_FILTER_LENGTH || ( imprecise_end - imprecise_pos ) * _params.ci_max_length > outer_end - inner_end ) {
			refine_point(chrom, inner_end, outer_end, imprecise_end, res, _params);
		}

		return res;
	};

};

#endif