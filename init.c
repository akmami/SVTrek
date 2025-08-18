#include "init.h"

void printUsage() {
    printf("Usage: ./svtrek [MODE] [OPTIONS]\n");
    printf("Mode:\n");
    printf("    disc    Variation discovery on graph alignment result.\n");
    printf("    audt    Audit the reported variations on VCF using BAM.\n");
}

void printDiscUsage() {
    printf("Usage: ./svtrek disc [-r|--gfa r/GFA] [-a|--gaf GAF] [-q|--fq FASTQ] [OPTIONS]\n");
    printf("Options:\n");
    printf("    [-o|--ouput] <filename>           Output filename [Default: svtrek.out]\n");
    printf("    -t <num>                          Thread number [Default: %d]\n", __THREAD_NUMBER);
    printf("    --verbose                         Verbose [Default: false]\n\n");
    printf("    --consensus-interval-range <num>  The interval to limit refinement range [DEFAULT: %d]\n", __CONSENSUS_INTERVAL_RANGE);
    printf("    --consensus-interval <num>        The interval that is considered into the same position [DEFAULT: %d]\n", __CONSENSUS_INTERVAL);
    printf("    --consensus-min-count <num>       Minimum number of elements needs for the consensus [Default: %d]\n\n", __CONSENSUS_MIN_COUNT);
}

void printAudtUsage() {
    printf("Usage: ./svtrek audt [-b|--bam BAM] [-v|--vcf VCF file] [OPTIONS]\n");
    printf("Options:\n");
    printf("    [-o|--ouput] <filename>           Output filename [Default: svtrek.out]\n");
    printf("    -t <num>                          Thread number [Default: %d]\n", __THREAD_NUMBER);
    printf("    --verbose                         Verbose [Default: false]\n\n");
    printf("    --wider-interval <num>            Interval for the offset of the reads to start [Default: %d]\n", __WIDER_INTERVAL);
    printf("    --median-interval <num>           Interval for the offset of the reads (for point) [Default: %d]\n", __MEDIAN_INTERVAL);
    printf("    --narrow-interval <num>           Interval for the offset of the reads to end [Default: %d]\n", __NARROW_INTERVAL);
    printf("    --consensus-interval-range <num>  The interval to limit refinement range [DEFAULT: %d]\n", __CONSENSUS_INTERVAL_RANGE);
    printf("    --consensus-interval <num>        The interval that is considered into the same position [DEFAULT: %d]\n", __CONSENSUS_INTERVAL);
    printf("    --consensus-min-count <num>       Minimum number of elements needs for the consensus [Default: %d]\n\n", __CONSENSUS_MIN_COUNT);
}

void validate_file(const char *filename, const char *message) {
    if (filename == NULL) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "[ERROR]: File couldn't be opened %s\n", filename);
    }

    fclose(file);
}

void init_audt(int argc, char *argv[], audt_args *params) {

    if (argc < 2) {
        printAudtUsage();
        exit(1);
    }

    params->bam_file = NULL;
    params->vcf_file = NULL;
    params->output_file = "svtrek.out";
    params->verbose = 0;
    params->tload_factor = __THREAD_POOL_LOAD_FACTOR;
    params->thread_number = __THREAD_NUMBER;
    params->wider_interval = __WIDER_INTERVAL;
    params->median_interval = __MEDIAN_INTERVAL;
    params->narrow_interval = __NARROW_INTERVAL;
    params->consensus_interval_range = __CONSENSUS_INTERVAL_RANGE;
    params->consensus_interval = __CONSENSUS_INTERVAL;
    params->consensus_min_count = __CONSENSUS_MIN_COUNT;
    params->mode = MODE_AUDT;

    struct option long_options[] = {
        {"bam", required_argument, NULL, 1},
        {"vcf", required_argument, NULL, 2},
        {"output", required_argument, NULL, 3},
        {"verbose", no_argument, NULL, 4},
        {"wider-interval", required_argument, NULL, 5},
        {"median-interval", required_argument, NULL, 6},
        {"narrow-interval", required_argument, NULL, 7},
        {"consensus-interval-range", required_argument, NULL, 8},
        {"consensus-interval", required_argument, NULL, 9},
        {"consensus-min-count", required_argument, NULL, 10},
        {"help", no_argument, NULL, 11},
        {NULL, 0, NULL, 0}
    };

    int opt;
    int long_index;

    while ((opt = getopt_long(argc, argv, "b:v:o:t:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'b':
                params->bam_file = optarg;
                break;
            case 1:
                params->bam_file = optarg;
                break;	
            case 'v':
                params->vcf_file = optarg;
                break;
            case 2:
                params->vcf_file = optarg;
                break;
            case 'o':
                params->output_file = optarg;
                break;
            case 3:
                params->output_file = optarg;
                break;
            case 4:
                params->verbose = 1;
                break;
            case 't':
                params->thread_number = atoi(optarg);
                break;
            case 5:
                params->wider_interval = atoi(optarg);
                break;
            case 6:
                params->median_interval = atoi(optarg);
                break;
            case 7:
                params->narrow_interval = atoi(optarg);
                break;
            case 8:
                params->consensus_interval_range = atoi(optarg);
                break;
            case 9:
                params->consensus_interval = atoi(optarg);
                break;
            case 10:
                params->consensus_min_count = atoi(optarg);
                break;
            case 'h':
                printAudtUsage();
                exit(EXIT_SUCCESS);
            case 11:
                printAudtUsage();
                exit(EXIT_SUCCESS);
            default:
                printf("[ERROR] Option %d is invalid.\n", opt);
                printAudtUsage();
                exit(EXIT_FAILURE);
        }
    }

    validate_file(params->bam_file, "[ERROR] BAM file is not provided.");
    validate_file(params->vcf_file, "[ERROR] VCF file is not provided.");
}

void init_disc(int argc, char *argv[], disc_args *params) {

    if (argc < 2) {
        printDiscUsage();
        exit(1);
    }

    params->gfa_file = NULL;
    params->gaf_file = NULL;
    params->fq_file = NULL;
    params->output_file = "svtrek.out";
    params->verbose = 0;
    params->tload_factor = __THREAD_POOL_LOAD_FACTOR;
    params->thread_number = __THREAD_NUMBER;
    params->consensus_interval_range = __CONSENSUS_INTERVAL_RANGE;
    params->consensus_interval = __CONSENSUS_INTERVAL;
    params->consensus_min_count = __CONSENSUS_MIN_COUNT;
    params->mode = MODE_DISC;

    struct option long_options[] = {
        {"gfa", required_argument, NULL, 1},
        {"gaf", required_argument, NULL, 2},
        {"fq", required_argument, NULL, 3},
        {"output", required_argument, NULL, 4},
        {"verbose", no_argument, NULL, 5},
        {"consensus-interval-range", required_argument, NULL, 6},
        {"consensus-interval", required_argument, NULL, 7},
        {"consensus-min-count", required_argument, NULL, 8},
        {"help", no_argument, NULL, 9},
        {NULL, 0, NULL, 0}
    };

    int opt;
    int long_index;

    while ((opt = getopt_long(argc, argv, "r:a:q:o:t:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'r':
                params->gfa_file = optarg;
                break;
            case 1:
                params->gfa_file = optarg;
                break;	
            case 'a':
                params->gaf_file = optarg;
                break;
            case 2:
                params->gaf_file = optarg;
                break;
            case 'q':
                params->fq_file = optarg;
                break;
            case 3:
                params->fq_file = optarg;
                break;
            case 'o':
                params->output_file = optarg;
                break;
            case 4:
                params->output_file = optarg;
                break;
            case 5:
                params->verbose = 1;
                break;
            case 't':
                params->thread_number = atoi(optarg);
                break;
            case 6:
                params->consensus_interval_range = atoi(optarg);
                break;
            case 7:
                params->consensus_interval = atoi(optarg);
                break;
            case 8:
                params->consensus_min_count = atoi(optarg);
                break;
            case 'h':
                printDiscUsage();
                exit(EXIT_SUCCESS);
            case 9:
                printDiscUsage();
                exit(EXIT_SUCCESS);
            default:
                printf("[ERROR] Option %d is invalid.\n", opt);
                printDiscUsage();
                exit(EXIT_FAILURE);
        }
    }

    validate_file(params->gfa_file, "[ERROR] r/GFA file is not provided.");
    validate_file(params->gaf_file, "[ERROR] GAF file is not provided.");
    validate_file(params->fq_file, "[ERROR] FASTQ file is not provided.");
}

