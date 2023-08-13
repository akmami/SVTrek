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

enum svtype {
    DELETION,
    INSERTION,
    INVERSION,
};

struct result {
    int refined_start;
    int refined_end;
    int svlen;
};


// MARK: CONSENSUS

int consensus(std::vector<int> lengths, params &_params) {

    int i;
    int consensus = -1, max_count = _params.consensus_min_count - 1, distance = INT_MAX;
    std::vector<int>::iterator ptr_i;
    std::vector<int>::iterator ptr_j;

    sort(lengths.begin(), lengths.end());

    for (ptr_i = lengths.begin(), i = 0; ptr_i < lengths.end(); ptr_i++, i++) {
        int sum = 0, count = 1;
        for (ptr_j = lengths.begin() + i; ptr_j < lengths.end() && *ptr_j <= *ptr_i + _params.consensus_interval; ptr_j++) {
            sum += (*ptr_j) - (*ptr_i);
            count++;
        }
        
        sum = (*ptr_i) + round(sum / count);
        
        if (count > max_count) {
            max_count = count;
            consensus = sum;
        }
    }
    
    return consensus;
};

int consensus_pos(std::vector<int> locations, int imprecise_pos, params &_params) {

    int i;
    int consensus = -1, max_count = _params.consensus_min_count - 1, distance = INT_MAX;
    std::vector<int>::iterator ptr_i;
    std::vector<int>::iterator ptr_j;

    sort(locations.begin(), locations.end());

    for (ptr_i = locations.begin(), i = 0; ptr_i < locations.end(); ptr_i++, i++) {
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

void refine_start(int chrom, int outer_start, int inner_start, int imprecise_pos, result &res, params &_params, svtype type=DELETION) {

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
            int reference_pos = pos;
            for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                // If deletion is large enough to be considered as sv
                // and it is covered in cigar string as the read is long enough.
                if ( outer_start <= reference_pos && reference_pos <= inner_start && bam_cigar_op( cigar[i]) == __CIGAR_DELETION && bam_cigar_oplen( cigar[i]) >= __SV_MIN_LENGTH ) {
                    start_positions.push_back(reference_pos+1);
                }

                if ( bam_cigar_op( cigar[i] ) != __CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != __CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_PADDING ) {
                    reference_pos += bam_cigar_oplen( cigar[i] );
                }
            }
            
            // Refinement base on Breakpoints
            if ( type == DELETION && bam_cigar_op( cigar[aln->core.n_cigar-1] ) == __CIGAR_SOFT_CLIP && outer_start <= reference_pos && reference_pos <= inner_start ) {
                start_positions.push_back(reference_pos+1);
            }

            if ( inner_start < reference_pos )
                    break;
        }
    } else {
        _params.verbose && std::cout << " Invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }

    bam_destroy1(aln);
    sam_itr_destroy(iter);

    res.refined_start = consensus_pos(start_positions, imprecise_pos, _params);
};


void refine_end(int chrom, int inner_end, int outer_end, int imprecise_pos, result &res, params &_params, svtype type=DELETION) {
    
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

            int reference_pos = pos;
            for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                if ( bam_cigar_op( cigar[i] ) != __CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != __CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_PADDING ) {
                    reference_pos += bam_cigar_oplen( cigar[i] );
                }

                // If deletion is large enough to be considered as sv
                // and it is covered in cigar string as the read is long enough.
                if ( type == DELETION && inner_end <= reference_pos && reference_pos <= outer_end && bam_cigar_op( cigar[i]) == __CIGAR_DELETION && bam_cigar_oplen( cigar[i]) >= __SV_MIN_LENGTH ) {
                    end_positions.push_back(reference_pos+1);
                }

                if ( outer_end < reference_pos )
                    break;
            }
        }
    } else {
        _params.verbose && std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }
    
    bam_destroy1(aln);
    sam_itr_destroy(iter);
    
    res.refined_end = consensus_pos(end_positions, imprecise_pos, _params);
};


void refine_point(int chrom, int start, int end, int imprecise_pos, result &res, params &_params, svtype type=INSERTION) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    bam1_t *aln = bam_init1();
    hts_itr_t *iter;

    iter = sam_itr_queryi( _params.htslib.bam_file_index, chrom - 1, start - _params.wider_interval, end + _params.narrow_interval);
    std::vector<int> positions;
    std::vector<int> lengths;
    
    if (iter) {
        while (sam_itr_next( _params.htslib.fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);

            // If read starts in between (start, end)
            if ( bam_cigar_op( cigar[0] ) == __CIGAR_SOFT_CLIP && start <= pos && pos <= end  ) {
                positions.push_back(pos+1);
            }
            
            // If read ends in between (start, end)
            int reference_pos = pos;
            for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                if ( bam_cigar_op( cigar[i] ) != __CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != __CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != __CIGAR_PADDING ) {
                    reference_pos += bam_cigar_oplen( cigar[i] );
                }

                if ( type == INSERTION && start <= reference_pos && reference_pos <= end && bam_cigar_op( cigar[i] ) == __CIGAR_INSERTION && bam_cigar_oplen( cigar[i] ) > __SV_MIN_LENGTH ) {
                    positions.push_back(reference_pos+1);
                    lengths.push_back( bam_cigar_oplen( cigar[i] ) );
                }
            }
            if ( start <= reference_pos && reference_pos <= end ) {
                positions.push_back(reference_pos+1);
            }
        }
    } else {
        _params.verbose && std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }

    bam_destroy1(aln);
    sam_itr_destroy(iter);
    
    res.refined_start = consensus_pos(positions, imprecise_pos, _params);
    if ( type == INSERTION ) {
        res.svlen = consensus(lengths, _params);
    }
};

#endif