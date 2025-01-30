#include "argp/argp.h"
#include "params.h"
#include "process_vcf.h"
#include <string>
#include <iostream>

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

    std::string _ci_max_length, _wider_interval, _narrow_interval, _consensus_interval, _consensus_min_count, _thread_number;

    parser("--ci-max-length") >> _ci_max_length;
    parser("--wider-interval") >> _wider_interval;
    parser("--narrow-interval") >> _narrow_interval;
    parser("--consensus-interval") >> _consensus_interval;
    parser("--consensus-min-count") >> _consensus_min_count;
    parser({"-t", "--thread_number"}) >> _thread_number;

    try {
        _params.ci_max_length = stof(_ci_max_length);
    } catch (...) {
        _params.ci_max_length = __CI_MAX_LENGTH;
    }

    try {
        _params.wider_interval = stoi(_wider_interval);
    } catch (...) {
        _params.wider_interval = __WIDER_INTERVAL;
    }

    try {
        _params.narrow_interval = stoi(_narrow_interval);
    } catch (...) {
        _params.narrow_interval = __NARROW_INTERVAL;
    }

    try {
        _params.consensus_interval = stoi(_consensus_interval);
    } catch (...) {
        _params.consensus_interval = __CONSENSUS_INTEVAL;
    }

    try {
        _params.consensus_min_count = stoi(_consensus_min_count);
    } catch (...) {
        _params.consensus_min_count = __CONSENSUS_MIN_COUNT;
    }

    try {
        _params.thread_number = stoi(_thread_number);
    } catch (...) {
        _params.thread_number = __THREAD_NUMBER;
    }

    _params.verbose && std::cout << "Program begins..." << std::endl;
    
    process_vcf(_params);

    _params.verbose && std::cout << "End of the program execution" << std::endl;
    
    return 0;
};
