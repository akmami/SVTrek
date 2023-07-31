#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include "argh.h"
#include "variables.h"
#include "process_line.cpp"


int run(char const *bam, char const *vcf, bool verbose);

// MARK: - main program

int main(int argc, char *argv[]){

    std::cout << "Program begins..." << std::endl;
    
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    std::string bam_str, vcf_str, output_str;
    bool verbose;
    
    if ( !( cmdl({"--bam", "b"}) >> bam_str ) ) {
        std::cout << "BAM file is missing." << std::endl;
        exit(-1);
    }
    
    if ( !( cmdl({"--vcf", "f"}) >> vcf_str ) ) {
        std::cout << "VCF file is missing." << std::endl;
        exit(-1);
    }

    if ( cmdl[{ "-v", "--verbose" }] ) {
        verbose = true;
    } else {
        verbose = false;
    }
    
    
    const char* bam = bam_str.c_str();
    const char* vcf = vcf_str.c_str();
    
    struct stat buffer_bam;
    if ( stat ( bam, &buffer_bam) != 0) {
        std::cout << "BAM file doesn't exist." << std::endl;
        return -1;
    }
    
    struct stat buffer_vcf;
    if ( stat ( vcf, &buffer_vcf) != 0) {
        std::cout << "VCF file doesn't exist." << std::endl;
        return -1;
    }
    
    run(bam, vcf, verbose);

    std::cout << "End of the program execution" << std::endl;
    
    return 0;
}

// MARK: - implementation of helper functions

int run(char const *bam, char const *vcf, bool verbose) {

    init_bam_var(bam);
    
    std::ifstream vcf_file(vcf);                   // open sv file
    std::string line;
    
    while (getline(vcf_file, line)) {
        if (line.at(0) == '#')
            continue;
        process_line(line, verbose);
    }
    
    vcf_file.close();
    deinitialize();
    
    return 0;
}
