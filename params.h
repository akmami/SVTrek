#ifndef PARAMS_H
#define PARAMS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>

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

#define __WIDER_INTERVAL           20000
#define __MEDIAN_INTERVAL          10000
#define __NARROW_INTERVAL          2000
#define __CONSENSUS_INTERVAL_RANGE 500
#define __CONSENSUS_INTEVAL        5
#define __CONSENSUS_MIN_COUNT      3
#define __SV_MIN_LENGTH            50
//#define CONSENSUS_COUNT_PERC     0.3

#define __THREAD_NUMBER           4
#define __THREAD_POOL_LOAD_FACTOR 2

typedef struct _args {
    // input arguments
    const char* bam_file;
    const char* vcf_file;
    const char* output_file;
    int thread_number;
    int verbose;
    int tload_factor;
    // program arguments
    int wider_interval;
    int median_interval;
    int narrow_interval;
    int consensus_interval_range;
    int consensus_interval;
    int consensus_min_count;
} args;

typedef struct _line_queue {
    char **lines;   /**< Queue to store lines extracted from VCF file for threads. */
    int size;       /**< The size of the queue. */
    int capacity;   /**< Capacity of the queue. */
    int front;      /**< The index for the pushing point. */
    int rear;       /**< The index for the popping point. */
} line_queue;

typedef struct _htslib_params {
    samFile *fp_in;
    bam_hdr_t *bam_hdr;
    hts_idx_t *bam_file_index;
} htslib_params;

typedef struct _t_arg {
    int wider_interval;
    int median_interval;
    int narrow_interval;
    int consensus_interval_range;
    int consensus_interval;
    int consensus_min_count;
    htslib_params hargs;
    line_queue *queue;
    pthread_mutex_t *queue_mutex;
    pthread_mutex_t *out_err_mutex;
    pthread_cond_t *cond_not_full;
    pthread_cond_t *cond_not_empty;
    int *exit_signal;
} t_arg;

typedef enum {
    SV_UNKNOWN = 0,
    SV_INS,    // Insertion
    SV_DEL,    // Deletion
    SV_INV,    // Inversion
    SV_DUP,    // Duplication
    SV_TRA,    // Translocation
    SV_BND     // Breakend
} sv_type_t;

typedef struct _interval {
    uint32_t start;
    uint32_t end;
} interval;

#endif // PARAMS_H
