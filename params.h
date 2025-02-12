#ifndef PARAMS_H
#define PARAMS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <fstream>

//                                Consumes     query   reference   Op
#define __CIGAR_ALIGNMENT_MATCH       0     //  yes     yes         M
#define __CIGAR_INSERTION             1     //  yes     no          I
#define __CIGAR_DELETION              2     //  no      yes         D
#define __CIGAR_SKIPPED               3     //  no      yes         N
#define __CIGAR_SOFT_CLIP             4     //  yes     no          S
#define __CIGAR_HARD_CLIP             5     //  no      no          H
#define __CIGAR_PADDING               6     //  no      no          P
#define __CIGAR_SEQUENCE_MATCH        7     //  yes     yes         =
#define __CIGAR_SEQUENCE_MISMATCH     8     //  yes     yes         X

#define __FLAG_MULTIPLE_SEGMENTS          0x1
#define __FLAG_SECONDARY_ALIGNMENT        0x100
#define __FLAG_SUPPLEMENTARY_ALIGNMENT    0x800

#define __CI_MIN_FILTER_LENGTH    5000
#define __CI_MAX_LENGTH           0.1
#define __WIDER_INTERVAL          40000
#define __NARROW_INTERVAL         2000
#define __CONSENSUS_INTEVAL       10
#define __CONSENSUS_MIN_COUNT     3
#define __SV_MIN_LENGTH           50
//#define CONSENSUS_COUNT_PERC    0.3

#define __THREAD_NUMBER           32

typedef struct _args {
    // input arguments
    const char* bam_file;
    const char* vcf_file;
    const char* output_file;
    int thread_number;
    bool verbose;
    // program arguments
    float ci_max_length;
    int wider_interval;
    int narrow_interval;
    int consensus_interval;
    int consensus_min_count;
} args;

typedef struct _htslib_params {
    samFile *fp_in;
    bam_hdr_t *bam_hdr;
    hts_idx_t *bam_file_index;
} htslib_params;

enum svtype {
    DELETION,
    INSERTION,
    INVERSION
};

#endif // PARAMS_H
