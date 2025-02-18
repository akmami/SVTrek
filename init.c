#include "init.h"

void printUsage() {
    printf("Usage: ./svtrek [-b|--bam BAM] [-v|--vcf VCF file] [OPTIONS]\n\n");
    printf("Options:\n");
    printf("\t[-o|--ouput] [filename]       \n\n");
    printf("\t-t [num]                      [Default: %d]\n\n", __THREAD_NUMBER);
    printf("\t--ci-max-length [num]         [Default: %f]\n\n", __CI_MAX_LENGTH);
    printf("\t--wider-interval [num]        [Default: %d]\n\n", __WIDER_INTERVAL);
    printf("\t--narrow-interval [num]       [Default: %d]\n\n", __NARROW_INTERVAL);
    printf("\t--consensus-interval [num]    [DEFAULT: %d]\n\n", __CONSENSUS_INTEVAL);
    printf("\t--consensus-min-count [num]   [Default: %d]\n\n", __CONSENSUS_MIN_COUNT);
    printf("\t--verbose                     \n\n");
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
    params->ci_max_length = __CI_MAX_LENGTH;
    params->narrow_interval = __NARROW_INTERVAL;
    params->consensus_interval = __CONSENSUS_INTEVAL;
    params->consensus_min_count = __CONSENSUS_MIN_COUNT;

    struct option long_options[] = {
        {"bam", required_argument, NULL, 1},
        {"vcf", required_argument, NULL, 2},
        {"output", required_argument, NULL, 3},
        {"verbose", required_argument, NULL, 4},
        {"ci-max-length", required_argument, NULL, 5},
        {"wider-interval", required_argument, NULL, 6},
        {"narrow-interval", required_argument, NULL, 7},
        {"consensus-interval", required_argument, NULL, 8},
        {"consensus-min-count", required_argument, NULL, 9},
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
                params->ci_max_length = atoi(optarg);
                break;
            case 6:
                params->wider_interval = atoi(optarg);
                break;
            case 7:
                params->narrow_interval = atoi(optarg);
                break;
            case 8:
                params->consensus_interval = atoi(optarg);
                break;
            case 9: 
                params->consensus_min_count = atoi(optarg);
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    validate_file(params->bam_file, "[ERROR] BAM file is not provided.");
    validate_file(params->vcf_file, "[ERROR] VCF file is not provided.");
}