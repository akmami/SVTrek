#include <string>
#include <iostream>
#include <fstream>
#include "argp/argp.h"
#include "static_variables.h"
#include "htslib.h"
#include "params.h"
#include "process_line.cpp"


int process_vcf(params &_params);


// MARK: - main program

int main(int argc, char *argv[]){

    argp::parser parser;
    parser.parse(argc, argv);
    params _params;

    std::string bam_file, vcf_file, output_file;

    parser({"--bam", "-b"}) >> bam_file;
    parser({"--vcf", "-f"}) >> vcf_file;
    parser({"--output", "-o"}) >> output_file;

    if ( bam_file == "" ) {
        std::cout << "BAM file is missing." << std::endl;
        exit(-1);
    }
    
    if ( vcf_file == "" ) {
        std::cout << "VCF file is missing." << std::endl;
        exit(-1);
    }

    if ( output_file == "" ) {
        output_file = "output.vcf";
    }

    std::ifstream bam(vcf_file);
    if( !bam.good() ) {
        std::cout << "Unable to open bam file." << std::endl;
        exit(-1);
    }
    
    std::ifstream vcf(vcf_file);
    if( !vcf.good() ) {
        std::cout << "Unable to open vcf file." << std::endl;
        exit(-1);
    }

    _params.bam_file = bam_file.c_str();
    _params.vcf_file = vcf_file.c_str();
    _params.output_file = output_file.c_str();
    _params.verbose = false;

    if ( parser[{ "-v", "--verbose" }] ) {
        _params.verbose = true;
    }

    std::string _ci_max_length, _wider_interval, _narrow_interval, _consensus_interval, _consensus_min_count;
    float ci_max_length = __CI_MAX_LENGTH;
    int wider_interval = __WIDER_INTERVAL;
    int narrow_interval = __NARROW_INTERVAL;
    int consensus_interval = __CONSENSUS_INTEVAL;
    int consensus_min_count = __CONSENSUS_MIN_COUNT;


    parser("--ci-max-length") >> _ci_max_length;
    parser("--wider-interval") >> _wider_interval;
    parser("--narrow-interval") >> _narrow_interval;
    parser("--consensus-interval") >> _consensus_interval;
    parser("--consensus-min-count") >> _consensus_min_count;

    try {
        _params.ci_max_length = stof(_ci_max_length);
    } catch (...) {}

    try {
        _params.wider_interval = stoi(_wider_interval);
    } catch (...) {}

    try {
        _params.narrow_interval = stoi(_narrow_interval);
    } catch (...) {}

    try {
        _params.consensus_interval = stoi(_consensus_interval);
    } catch (...) {}

    try {
        _params.consensus_min_count = stoi(_consensus_min_count);
    } catch (...) {}

    _params.verbose && std::cout << "Program begins..." << std::endl;
    
    process_vcf(_params);

    _params.verbose && std::cout << "End of the program execution" << std::endl;
    
    return 0;
}

// MARK: - implementation of helper functions

int process_vcf(params &_params) {

    initialize_params(_params);
    
    std::ifstream vcf(_params.vcf_file);                // open sv file
    std::string line;
    
    while (getline(vcf, line)) {
        if (line.at(0) == '#') {
            _params.out_vcf << line << std::endl;
            continue;
        }
        process_line(_params, line);
    }
    
    vcf.close();
    deinitialize_params(_params);
    
    return 0;
}
