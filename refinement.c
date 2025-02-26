#include "refinement.h"

int lower_bound(int *arr, int size, int location) {
    for (int i = 0; i < size; i++) {
        if (arr[i] > location) {
            return i == 0 ? 0 : i-1;
        }
    }
    return size - 1;
}

int upper_bound(int *arr, int size, int location) {
    for (int i = 0; i < size; i++) {
        if (arr[i] < location) {
            return i;
        }
    }
    return size - 1;
}

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

int consensus_pos(int *locations, int size, int pos, int consensus_min_count, int consensus_interval, int consensus_interval_range) {
    
    if (size < consensus_min_count) {
        return -1;
    }

    int consensus_val_left = -1;
    int max_count_left = consensus_min_count - 1;
    int distance_left = 0x7fffffff;
    int consensus_val_right = -1;
    int max_count_right = consensus_min_count - 1;
    int distance_right = 0x7fffffff;
    
    quicksort(locations, 0, size-1);

    int point = lower_bound(locations, size, pos + __SV_MIN_LENGTH / 2);

    for (int i = point; 0 <= i && abs(pos-locations[i]) < consensus_interval_range; i--) {
        int count = 1;
        uint64_t total = locations[i];
        for (int j = i - 1; 0 <= j && locations[i] <= locations[j] + consensus_interval; j--) {
            count++;
            total += locations[j];
        }
        int candidate = ( total + count / 2 ) / count;

        if (count > max_count_left) {
            if (abs(pos - candidate) < consensus_interval) {
                return candidate;
            } else if (abs(pos - candidate) < distance_left) {
                max_count_left = count;
                consensus_val_left = candidate;
                distance_left = abs(pos - candidate);
            }
        }
    }

    point = upper_bound(locations, size, pos - __SV_MIN_LENGTH / 2);

    for (int i = point; i < size && abs(pos-locations[i]) < consensus_interval_range; i++) {
        int count = 1;
        uint64_t total = locations[i];
        for (int j = i + 1; j < size && locations[j] <= locations[i] + consensus_interval; j++) {
            count++;
            total += locations[j];
        }
        int candidate = ( total + count / 2 ) / count;

        if (count > max_count_right) {
            if (abs(pos - candidate) < consensus_interval) {
                return candidate;
            } else if (abs(pos - candidate) < distance_right) {
                max_count_right = count;
                consensus_val_right = candidate;
                distance_right = abs(pos - candidate);
            }
        }
    }

    return distance_left < distance_right ? consensus_val_left : consensus_val_right;
}

int refine_start(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *start_locations = (int *)malloc(sizeof(int)*capacity);
    if(start_locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for start positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, inter.start-1, inter.end-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            int check_soft_clip = bam_cigar_op(cigar[aln->core.n_cigar-1]) == __CIGAR_SOFT_CLIP;

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
                    start_locations[size++] = reference_pos;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    check_soft_clip = 0;
                    break;
                }
            }
            
            if (check_soft_clip && inter.start <= reference_pos && reference_pos <= inter.end) {
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
                start_locations[size++] = reference_pos;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    // pick the best refined start
    return consensus_pos(start_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
}

int refine_end(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *end_locations = (int *)malloc(sizeof(int)*capacity);
    if(end_locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for end positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, inter.start-1, inter.end-1);

    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (sv_type == SV_DEL && bam_cigar_op(cigar[i]) == __CIGAR_DELETION && __SV_MIN_LENGTH < bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        int *temp = (int *)realloc(end_locations, sizeof(int)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate end locations array.\n");
                            return -1;
                        }
                        end_locations = temp;
                    }
                    end_locations[size++] = reference_pos + bam_cigar_oplen(cigar[i]) + 1;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    break;
                }
            }

            if (bam_cigar_op(cigar[0]) == __CIGAR_SOFT_CLIP && inter.start <= aln->core.pos && aln->core.pos <= inter.end) {
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
    return consensus_pos(end_locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
}

int refine_point(sv_type_t sv_type, int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *locations = (int *)malloc(sizeof(int)*capacity);
    if(locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, inter.start-1, inter.end-1);
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
                    locations[size++] = reference_pos;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    break;
                }
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);
    return consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
}

int refine_ins(int chrom, interval inter, uint32_t imprecise_pos, t_arg *params) {

    int capacity = 100;
    int size = 0;
    int *locations = (int *)malloc(sizeof(int)*capacity);
    if(locations == NULL) {
        fprintf(stderr, "Couldn't allocate array for positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();

    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom-1, inter.start-1, inter.end-1);
    int count = 0;
    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            count++;
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar   = bam_get_cigar(aln);

            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (bam_cigar_op(cigar[i]) == __CIGAR_INSERTION && __SV_MIN_LENGTH <= bam_cigar_oplen(cigar[i])) {
                    if (capacity == size) {
                        capacity = capacity * 1.5;
                        int *temp = (int *)realloc(locations, sizeof(int)*capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate locations array.\n");
                            return -1;
                        }
                        locations = temp;
                    }
                    locations[size++] = reference_pos;
                }

                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION && bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }

                if (reference_pos > inter.end) {
                    break;
                }
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);
    return consensus_pos(locations, size, imprecise_pos, params->consensus_min_count, params->consensus_interval, params->consensus_interval_range);
}

void deletion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter) {
    res_inter->start = refine_start(SV_DEL, chrom, begin, sv_inter.start, params);
    res_inter->end = refine_end(SV_DEL, chrom, end, sv_inter.end, params);
}

void insertion(int chrom, interval begin, uint32_t pos, t_arg *params, uint32_t *res_start) {
    (*res_start) = refine_ins(chrom, begin, pos, params);
}

void inversion(int chrom, interval begin, interval end, interval sv_inter, t_arg *params, interval *res_inter) {
    res_inter->start = refine_point(SV_INV, chrom, begin, sv_inter.start, params);
    res_inter->end = refine_point(SV_INV, chrom, end, sv_inter.end, params);
}
