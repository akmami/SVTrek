#ifndef HTSLIB_H
#define HTSLIB_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>

struct htslib_params {
    samFile *fp_in;
    bam_hdr_t *bam_hdr;
    hts_idx_t *bam_file_index;
};

#endif