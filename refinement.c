#include "refinement.h"

void quicksort(int *array, int low, int high) {
    if (low < high) {
        int pivot = array[high];
        int i = low - 1;

        for (int j = low; j < high; j++) {
            if (array[j] < pivot) {
                i++;
                int temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }

        int temp = array[i+1];
        array[i+1] = array[high];
        array[high] = temp;

        quicksort(array, low, i);
        quicksort(array, i + 2, high);
    }
}

/**
 * consensus() picks the most common (or "peak") value in a range of integers
 */
int consensus(int *arr, int size, int consensus_min_count, int consensus_interval) {
    int consensus_val = -1;
    int max_count = consensus_min_count - 1;

    quicksort(arr, 0, size-1);

    for (int i = 0; i < size; i++) {
        int count = 1;
        // count how many values fall within consensus_interval of lengths[i]
        for (int j = i + 1; j < size && arr[j] <= arr[i] + consensus_interval; j++) {
            count++;
        }
        if (count > max_count) {
            max_count = count;
            consensus_val = arr[i];
        }
    }
    return consensus_val;
}

/**
 * consensus_pos() is similar, but picks the best cluster *and*
 * tries to be closest to imprecise_pos
 */
int consensus_pos(int *locations, int size, int imprecise_pos, int consensus_min_count, int consensus_interval) {
    int consensus_val = -1;
    int max_count = consensus_min_count - 1;
    int distance = 0x7fffffff;

    quicksort(locations, 0, size-1);

    for (int i = 0; i < size; i++) {
        int count = 1;
        for (int j = i + 1; j < size && locations[j] <= locations[i] + consensus_interval; j++) {
            count++;
        }
        int candidate = locations[i];

        if (count > max_count || (count == max_count && abs(imprecise_pos - candidate) < distance)) {
            max_count = count;
            consensus_val = candidate;
            distance = abs(imprecise_pos - candidate);
        }
    }
    return consensus_val;
}

/**
 * refine_start()
 */
int refine_start(sv_type_t sv_type, int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *start_locations = (int *)malloc(sizeof(int)*capacity);
    if(start_locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for start positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    // query the region from [outer_start - wider_interval] to [inner_start + narrow_interval]
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, outer_start-1, inner_start-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);

            // n_cigar is unsigned, so use uint32_t in the loop
            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (sv_type == SV_DEL && bam_cigar_op(cigar[i]) == __CIGAR_DELETION && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        int *temp = (int *)realloc(start_locations, sizeof(int)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate end locations array.\n");
                            return -1;
                        }
                        start_locations = temp;
                    }
                    start_locations[size++] = reference_pos + 1;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inner_start) {
                    break;
                }
            }

            if (bam_cigar_op(cigar[0]) == __CIGAR_SOFT_CLIP && outer_start <= bam_cigar_oplen(cigar[0]) && bam_cigar_oplen(cigar[0]) <= inner_start) {
                if (capacity == size) {
                    capacity = capacity * 1.5;
                    int *temp = (int *)realloc(start_locations, sizeof(int)*capacity);
                    if (temp == NULL) {
                        fprintf(stderr, "[ERROR] Couldn't reallocate start locations array.\n");
                        return -1;
                    }
                    start_locations = temp;
                }
                // +1 because positions are 1-based
                start_locations[size++] = reference_pos + 1;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    // pick the best refined start
    return consensus_pos(start_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval);
}

/**
 * refine_end()
 */
int refine_end(sv_type_t sv_type, int chrom, uint32_t inner_end, uint32_t outer_end, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *end_locations = (int *)malloc(sizeof(int)*capacity);
    if(end_locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for end positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, inner_end-1, outer_end-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            uint32_t i = 0;

            for (; i < aln->core.n_cigar; i++) {
                if (sv_type == SV_INS && bam_cigar_op(cigar[i]) == __CIGAR_INSERTION && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        int *temp = (int *)realloc(end_locations, sizeof(int)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate end locations array.\n");
                            return -1;
                        }
                        end_locations = temp;
                    }
                    end_locations[size++] = reference_pos + 1;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > outer_end) {
                    break;
                }
            }

            if (bam_cigar_op(cigar[i]) == __CIGAR_SOFT_CLIP && inner_end <= reference_pos && reference_pos <= outer_end) {
                if (capacity == size) {
                    capacity = capacity * 1.5;
                    int *temp = (int *)realloc(end_locations, sizeof(int)*capacity);
                    if (temp == NULL) {
                        fprintf(stderr, "[ERROR] Couldn't reallocate end locations array.\n");
                        return -1;
                    }
                    end_locations = temp;
                }
                end_locations[size++] = reference_pos + 1;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    // pick the best refined end
    return consensus_pos(end_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval);
}

/**
 * refine_point() - used for insertions or single-point refinement
 */
int refine_point(sv_type_t sv_type, int chrom, uint32_t start, uint32_t end, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *locations = (int *)malloc(sizeof(int)*capacity);
    if(locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, start-1, end-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar   = bam_get_cigar(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (sv_type == SV_INS && bam_cigar_op(cigar[i]) == __CIGAR_DELETION && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        int *temp = (int *)realloc(locations, sizeof(int)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate locations array.\n");
                            return -1;
                        }
                        locations = temp;
                    }
                    locations[size++] = reference_pos + 1;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > end) {
                    break;
                }
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    return consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval);
}

void deletion(int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t inner_end, uint32_t outer_end, uint32_t imprecise_pos, uint32_t imprecise_end, t_arg *params, uint32_t *res_start, uint32_t *res_end) {
    (*res_start) = refine_start(SV_DEL, chrom, outer_start, inner_start, imprecise_pos, params);
    (*res_end) = refine_end(SV_DEL, chrom, inner_end, outer_end, imprecise_end, params);
}

void insertion(int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t imprecise_pos, t_arg *params, uint32_t *res_start) {
    (*res_start) = refine_point(SV_INS, chrom, outer_start, inner_start, imprecise_pos, params);
}

void inversion(int chrom, uint32_t outer_start, uint32_t inner_start, uint32_t inner_end, uint32_t outer_end, uint32_t imprecise_pos, uint32_t imprecise_end, t_arg *params, uint32_t *res_start, uint32_t *res_end) {
    (*res_start) = refine_point(SV_INV, chrom, outer_start, inner_start, imprecise_pos, params);
    (*res_end) = refine_point(SV_INV, chrom, inner_end, outer_end, imprecise_end, params);
}
