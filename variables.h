#ifndef STATIC_H
#define STATIC_H

//                                Consumes query   reference   Op
#define CIGAR_ALIGNMENT_MATCH       0   //  yes     yes         M
#define CIGAR_INSERTION             1   //  yes     no          I
#define CIGAR_DELETION              2   //  no      yes         D
#define CIGAR_SKIPPED               3   //  no      yes         N
#define CIGAR_SOFT_CLIP             4   //  yes     no          S
#define CIGAR_HARD_CLIP             5   //  no      no          H
#define CIGAR_PADDING               6   //  no      no          P
#define CIGAR_SEQUENCE_MATCH        7   //  yes     yes         =
#define CIGAR_SEQUENCE_MISMATCH     8   //  yes     yes         X

#define FLAG_MULTIPLE_SEGMENTS          0x1
#define FLAG_SECONDARY_ALIGNMENT        0x100
#define FLAG_SUPPLEMENTARY_ALIGNMENT    0x800

#define CI_MAX_LENGTH           0.1
#define WIDER_INTERVAL          40000
#define NARROW_INTERVAL         2000
#define CONSENSUS_INTEVAL       10
#define CONSENSUS_MIN_COUNT     1
//#define CONSENSUS_COUNT_PERC    0.3

#include "htslib.h"


samFile *fp_in;
bam_hdr_t *bamHdr;
hts_idx_t *bam_file_index;
bam1_t *aln;

bool isInit;

void init_bam_var(char const *bam) {

	if (isInit)
		return;

	fp_in = hts_open(bam, "r");             //open bam file
    bamHdr = sam_hdr_read(fp_in);           //read header
    bam_file_index = sam_index_load( fp_in, bam );
    aln = bam_init1();                      //initialize an alignment

    isInit = true;
}

void deinitialize() {

	if (!isInit)
		return;

    bam_destroy1(aln);
    sam_close(fp_in);
    hts_idx_destroy(bam_file_index);
    bam_hdr_destroy(bamHdr);

    fp_in = NULL;
    bamHdr = NULL;
    bam_file_index = NULL;
    aln = NULL;

    isInit = false;
}

#endif