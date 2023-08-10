#ifndef PARAMS_H
#define PARAMS_H

#include "htslib.h"
#include <fstream>

struct htslib_params {
    samFile *fp_in;
    bam_hdr_t *bam_hdr;
    hts_idx_t *bam_file_index;
};

struct params {
    const char* bam_file;
    const char* vcf_file;
    const char* output_file;
    std::ofstream out_vcf;

    bool verbose;
    htslib_params htslib;

    float ci_max_length;
    int wider_interval;
    int narrow_interval;
    int consensus_interval;
    int consensus_min_count;

    bool isInit;
};

void initialize_params(params &_params) {

	if (_params.isInit)
		return;

	_params.htslib.fp_in = hts_open(_params.bam_file, "r");            //open bam file
    _params.htslib.bam_hdr = sam_hdr_read(_params.htslib.fp_in);    //read header
    _params.htslib.bam_file_index = sam_index_load( _params.htslib.fp_in, _params.bam_file );

    _params.out_vcf.open(_params.output_file);

    _params.isInit = true;
};

void deinitialize_params(params &_params) {

	if (!_params.isInit)
		return;

    sam_close(_params.htslib.fp_in);
    hts_idx_destroy(_params.htslib.bam_file_index);
    bam_hdr_destroy(_params.htslib.bam_hdr);

    _params.htslib.fp_in = NULL;
    _params.htslib.bam_hdr = NULL;
    _params.htslib.bam_file_index = NULL;

    _params.out_vcf.close();
    
    _params.isInit = false;
};

#endif