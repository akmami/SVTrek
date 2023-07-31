#ifndef INSERTION_CPP
#define INSERTION_CPP

#include <iostream>
#include <string>
#include "refinement.cpp"

bool insertion(std::string id, int chrom, int outer_start, int inner_start, int imprecise_pos, bool verbose=false) {

    int refined_start = find_start_or_end(chrom, outer_start, inner_start, imprecise_pos);
    
    std::cout << chrom << "\t" << id << "\t" << "ins" << "\t" << refined_start << std::endl;
    
    return true;
}

#endif