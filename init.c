#include "init.h"

void printUsage() {
    printf("Usage: ./svtrek [-b|--bam BAM] [-v|--vcf VCF file] [OPTIONS]\n\n");
    printf("Options:\n");
    printf("\t[-o|--ouput] <filename>           Output filename [Default: svtrek.out]\n");
    printf("\t-t <num>                          Thread number [Default: %d]\n", __THREAD_NUMBER);
    printf("\t--verbose                         Verbose [Default: false]\n\n");
    printf("\t--wider-interval <num>            Interval for the offset of the reads to start [Default: %d]\n", __WIDER_INTERVAL);
    printf("\t--median-interval <num>           Interval for the offset of the reads (for point) [Default: %d]\n", __MEDIAN_INTERVAL);
    printf("\t--narrow-interval <num>           Interval for the offset of the reads to end [Default: %d]\n", __NARROW_INTERVAL);
    printf("\t--consensus-interval-range <num>  The interval to limit refinement range [DEFAULT: %d]\n", __CONSENSUS_INTERVAL_RANGE);
    printf("\t--consensus-interval <num>        The interval that is considered into the same position [DEFAULT: %d]\n", __CONSENSUS_INTEVAL);
    printf("\t--consensus-min-count <num>       Minimum number of elements needs for the consensus [Default: %d]\n\n", __CONSENSUS_MIN_COUNT);
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

void init(int argc, char *argv[], args *params) {

    if (argc < 2) {
        printUsage();
        exit(1);
    }

    params->bam_file = NULL;
    params->vcf_file = NULL;
    params->output_file = "svtrek.out"; // for now
    params->verbose = 0;
    params->tload_factor = __THREAD_POOL_LOAD_FACTOR;
    params->thread_number = __THREAD_NUMBER;
    params->wider_interval = __WIDER_INTERVAL;
    params->median_interval = __MEDIAN_INTERVAL;
    params->narrow_interval = __NARROW_INTERVAL;
    params->consensus_interval_range = __CONSENSUS_INTERVAL_RANGE;
    params->consensus_interval = __CONSENSUS_INTEVAL;
    params->consensus_min_count = __CONSENSUS_MIN_COUNT;

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
        {NULL, 0, NULL, 0}
    };

    int opt;
    int long_index;

    while ((opt = getopt_long(argc, argv, "b:v:o:t:", long_options, &long_index)) != -1) {
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
            default:
                exit(EXIT_FAILURE);
        }
    }

    validate_file(params->bam_file, "[ERROR] BAM file is not provided.");
    validate_file(params->vcf_file, "[ERROR] VCF file is not provided.");
}