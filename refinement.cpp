#ifndef REFINEMENT_CPP
#define REFINEMENT_CPP

#include <iostream>
#include "htslib.h"
#include "static_variables.h"
#include "params.h"
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

// MARK: CONSENSUS

int consensus(std::vector<int> locations, int imprecise_pos, params &_params) {

    int length = locations.size(), i;
    int consensus = -1, max_count = _params.consensus_min_count - 1, distance = INT_MAX;
    std::vector<int>::iterator ptr_i;
    std::vector<int>::iterator ptr_j;

    sort(locations.begin(), locations.end());

    for (ptr_i = locations.begin(), i = 0; ptr_i < locations.end(); ptr_i++, i++) {
        int item_i = *ptr_i;
        int sum = 0, count = 1;
        for (ptr_j = locations.begin() + i; ptr_j < locations.end() && *ptr_j <= *ptr_i + _params.consensus_interval; ptr_j++) {
            sum += (*ptr_j) - (*ptr_i);
            count++;
        }
        
        sum = (*ptr_i) + round(sum / count);
        
        if (count > max_count) {
            max_count = count;
            consensus = sum;
            distance = std::abs(imprecise_pos-sum);
        } else if (count == max_count && distance > std::abs(imprecise_pos-sum)) {
            max_count = count;
            consensus = sum;
            distance = std::abs(imprecise_pos-sum);
        }
    }
    
    return consensus;
};

// MARK: REFINEMENT

int refine_start(int chrom, int outer_start, int inner_start, int imprecise_pos, params &_params) {

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    bam1_t *aln = bam_init1();
    hts_itr_t *iter;

    iter = sam_itr_queryi( _params.htslib.bam_file_index, chrom - 1, outer_start - _params.wider_interval, inner_start + _params.narrow_interval);
    std::vector<int> start_positions;
    
    if (iter) {
        while (sam_itr_next( _params.htslib.fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read ends in between (outer_start, inner_start)
            if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == __CIGAR_SOFT_CLIP) {
                int reference_pos = pos;
                for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                    if ( bam_cigar_op( cigar[i] ) != __CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != __CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_PADDING ) {
                        reference_pos += bam_cigar_oplen( cigar[i] );
                    }
                }
                
                if ( outer_start <= reference_pos && reference_pos <= inner_start ) {
                    start_positions.push_back(reference_pos+1);
                }
            }
        }
    } else {
        std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }

    bam_destroy1(aln);
    sam_itr_destroy(iter);

    return consensus(start_positions, imprecise_pos, _params);
};


int refine_end(int chrom, int inner_end, int outer_end, int imprecise_pos, params &_params) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    bam1_t *aln = bam_init1();
    hts_itr_t *iter;

    std::vector<int> end_positions;
    iter = sam_itr_queryi( _params.htslib.bam_file_index, chrom - 1, inner_end - _params.wider_interval, outer_end + _params.narrow_interval);
    
    if (iter) {
        while (sam_itr_next( _params.htslib.fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read starts in between (inner_end, outer_end)
            if ( bam_cigar_op( cigar[0] ) == __CIGAR_SOFT_CLIP && inner_end <= pos && pos <= outer_end  ) {
                end_positions.push_back(pos+1);
            }
        }
    } else {
        std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }
    
    bam_destroy1(aln);
    sam_itr_destroy(iter);
    
    return consensus(end_positions, imprecise_pos, _params);
};


int refine_point(int chrom, int start, int end, int imprecise_pos, params &_params) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    bam1_t *aln = bam_init1();
    hts_itr_t *iter;

    iter = sam_itr_queryi( _params.htslib.bam_file_index, chrom - 1, start - _params.wider_interval, end + _params.narrow_interval);
    std::vector<int> positions;
    
    if (iter) {
        while (sam_itr_next( _params.htslib.fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read ends in between (start, end)
            if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == __CIGAR_SOFT_CLIP) {
                int reference_pos = pos;
                for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                    if ( bam_cigar_op( cigar[i] ) != __CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != __CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_PADDING ) {
                        reference_pos += bam_cigar_oplen( cigar[i] );
                    }
                }
                if ( start <= reference_pos && reference_pos <= end ) {
                    positions.push_back(reference_pos+1);
                }
            }
            // If read starts in between (start, end)
            if ( bam_cigar_op( cigar[0] ) == __CIGAR_SOFT_CLIP && start <= pos && pos <= end ) {
                positions.push_back(pos);
            }
        }
    } else {
        std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }

    bam_destroy1(aln);
    sam_itr_destroy(iter);
    
    return consensus(positions, imprecise_pos, _params);
};

#endif