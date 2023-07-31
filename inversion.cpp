#ifndef INVERSION_CPP
#define INVERSION_CPP

#include <iostream>
#include <string>
#include "refinement.cpp"

bool inversion(std::string id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end, int imprecise_pos, int imprecise_end, bool verbose=false) {

    // sam_read1 variables
    int refined_start = find_start_or_end(chrom, outer_start, inner_start, imprecise_pos);
    int refined_end = find_start_or_end(chrom, inner_end, outer_end, imprecise_end);
    
    std::cout << chrom << "\t" << id << "\t" << "inv" << "\t" << refined_start << "\t" <<  refined_end << std::endl;
    
    return true;
}

#endif