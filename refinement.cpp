#ifndef REFINEMENT_CPP
#define REFINEMENT_CPP

#include <iostream>
#include "variables.h"
#include "htslib.h"
#include "consensus.cpp"

int find_start(int chrom, int outer_start, int inner_start, int imprecise_pos) {

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
    hts_itr_t *iter;
    iter = sam_itr_queryi( bam_file_index, chrom - 1, outer_start - WIDER_INTERVAL, inner_start + NARROW_INTERVAL);
    std::vector<int> start_positions;
    
    if (iter) {
        while (sam_itr_next( fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read ends in between (outer_start, inner_start)
            if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == CIGAR_SOFT_CLIP) {
                int reference_pos = pos;
                for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                    if ( bam_cigar_op( cigar[i] ) != CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_PADDING ) {
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
    sam_itr_destroy(iter);

    return consensus(start_positions, imprecise_pos);
};


int find_end(int chrom, int inner_end, int outer_end, int imprecise_pos) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
    hts_itr_t *iter;
    std::vector<int> end_positions;
    iter = sam_itr_queryi( bam_file_index, chrom - 1, inner_end - NARROW_INTERVAL, outer_end + NARROW_INTERVAL);
    
    if (iter) {
        while (sam_itr_next( fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read starts in between (inner_end, outer_end)
            if ( bam_cigar_op( cigar[0] ) == CIGAR_SOFT_CLIP && inner_end <= pos && pos <= outer_end  ) {
                end_positions.push_back(pos+1);
            }
        }
    } else {
        std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }
    
    sam_itr_destroy(iter);
    
    return consensus(end_positions, imprecise_pos);
};


int find_start_or_end(int chrom, int start, int end, int imprecise_pos) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    hts_itr_t *iter;
    iter = sam_itr_queryi( bam_file_index, chrom - 1, start - NARROW_INTERVAL, end + NARROW_INTERVAL);
    std::vector<int> positions;
    
    if (iter) {
        while (sam_itr_next( fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read ends in between (start, end)
            if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == CIGAR_SOFT_CLIP) {
                int reference_pos = pos;
                for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                    if ( bam_cigar_op( cigar[i] ) != CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_PADDING ) {
                        reference_pos += bam_cigar_oplen( cigar[i] );
                    }
                }
                if ( start <= reference_pos && reference_pos <= end ) {
                    positions.push_back(reference_pos+1);
                }
            }
            // If read starts in between (start, end)
            if ( bam_cigar_op( cigar[0] ) == CIGAR_SOFT_CLIP && start <= pos && pos <= end ) {
                positions.push_back(pos);
            }
        }
    } else {
        std::cout << "# invalid interval, iter is null. pos: " << imprecise_pos << std::endl;
    }
    sam_itr_destroy(iter);
    
    return consensus(positions, imprecise_pos);
};


#endif