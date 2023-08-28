#ifndef PARAMS_H
#define PARAMS_H

#include <fstream>

struct params {
    const char* bam_file;
    const char* vcf_file;
    const char* output_file;
    std::ofstream out_vcf;

    bool verbose;

    float ci_max_length;
    int wider_interval;
    int narrow_interval;
    int consensus_interval;
    int consensus_min_count;
    int thread_number;

    bool isInit;
};

void initialize_params(params &_params) {

	if (_params.isInit)
		return;

    _params.out_vcf.open(_params.output_file);

    _params.isInit = true;
};

void deinitialize_params(params &_params) {

	if (!_params.isInit)
		return;

    _params.out_vcf.close();
    
    _params.isInit = false;
};

#endif