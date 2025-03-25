#include <stdlib.h>
#include <stdio.h>
#include <htslib/sam.h>
#include "sliding_window.h"
#include "refinement.h" 

void quicksort(int *array, int low, int high);

int sliding_window_ins(int chrom, interval inter, uint32_t imprecise_pos, t_arg *params, int window_size, int slide_Size) {
    int capacity = 100;
    int size = 0;
    int *locations = (int *)malloc(sizeof(int) * capacity);
    if (locations == NULL) {
        fprintf(stderr, "[ERROR] Couldn't allocate array for positions.\n");
        return -1;
    }

    bam1_t *aln = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(params->hargs.bam_file_index, chrom - 1, inter.start - 1, inter.end - 1);
    if (iter) {
        while (sam_itr_next(params->hargs.fp_in, iter, aln) > 0) {
            uint32_t reference_pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
                if (bam_cigar_op(cigar[i]) == __CIGAR_INSERTION &&
                    bam_cigar_oplen(cigar[i]) >= __SV_MIN_LENGTH) {
                    if (size == capacity) {
                        capacity = (int)(capacity * 1.5);
                        int *temp = (int *)realloc(locations, sizeof(int) * capacity);
                        if (temp == NULL) {
                            fprintf(stderr, "[ERROR] Couldn't reallocate locations array.\n");
                            free(locations);
                            bam_destroy1(aln);
                            return -1;
                        }
                        locations = temp;
                    }
                    locations[size++] = reference_pos;
                }
                if (bam_cigar_op(cigar[i]) != __CIGAR_INSERTION &&
                    bam_cigar_op(cigar[i]) != __CIGAR_SOFT_CLIP) {
                    reference_pos += bam_cigar_oplen(cigar[i]);
                }
                if (reference_pos > inter.end)
                    break;
            }
        }
        sam_itr_destroy(iter);
    }
    bam_destroy1(aln);

    //sliding window principle
    if (size == 0) {
        free(locations);
        return -1;
    }
    
    quicksort(locations, 0, size - 1);
    
    int best_cons = -1;
    int highest_sup = 0;
    
    for (int i = 0; i < size; i += slide_size) {
        int end = i;
        while (end < size && (locations[end] - locations[i]) <= window_size) {
            end++;
        }
        int compare = end - i;
        if (compare >= params->consensus_min_count && compare > highest_sup) {
            highest_sup = compare;
            int sum = 0;
            for (int j = i; j < end; j++) {
                sum += locations[j];
            }
            best_cons = (sum + compare / 2) / compare; 
        
        }
    }
    //printing once
    if (best_cons != -1) {
        printf("INS Discovery at position %d\n", best_cons);
    }
    return best_cons;
    
}