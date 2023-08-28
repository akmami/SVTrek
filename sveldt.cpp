#include <string>
#include <iostream>
#include "argp/argp.h"
#include "static_variables.h"
#include "htslib.h"
#include "params.h"
#include "thread_data.h"
#include "threads.cpp"


// MARK: - Helper functions' declaration

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


// MARK: - implementation of helper functions

int process_vcf(params &_params) {

    initialize_params(_params);
    
    std::ifstream vcf(_params.vcf_file);                // open sv file
    std::string line;

    // Multithreading
    _params.verbose && std::cout << "Initializing threads..." << std::endl;
    
    // Create all the threads
    thread_data threads[_params.thread_number];
    int thread_index;
    
    // Creating mutex for signaling from child to parent
    std::mutex _mutex;
    std::condition_variable _condition_variable_parent;
    std::unique_lock<std::mutex> _unique_lock(_mutex);
    
    for ( thread_index = 0; thread_index < _params.thread_number; thread_index++ ) {
        threads[thread_index].thread_id = thread_index;
        threads[thread_index]._thread = std::thread(thread_func, &_params, &(threads[thread_index]), &_condition_variable_parent);
    }

    _params.verbose && std::cout << "Total number of threads being initialized is " << _params.thread_number << std::endl;
    
    thread_index = 0;
    int line_index = 0;

    while (getline(vcf, line)) {
        
        line_index++;
        
        if (line.at(0) == '#') {
            _params.out_vcf << line << std::endl;
            continue;
        }

        // check if previous thread is available
        if (threads[thread_index]._run ) {
            _condition_variable_parent.wait(_unique_lock, [&] { return threads[thread_index]._run ? false : true; });
        }
        // Print processed value
        threads[thread_index]._available && _params.out_vcf << threads[thread_index].processed << std::endl;
        threads[thread_index]._available = false;

        //_params.verbose && std::cout << "Line is processed.. id: " << thread_index << std::endl;

        threads[thread_index].line = line;
        threads[thread_index]._run = true;
        threads[thread_index]._condition_variable_thread.notify_one();

        thread_index = ( thread_index + 1 ) % _params.thread_number;
    }

     _params.verbose && std::cout << "Ending all threads..." << std::endl;

    // Wait for all the tasks to be completed...
    for ( int temp = thread_index; temp < thread_index + _params.thread_number; temp++ ) { 
        // check if previous thread is available
        if (threads[temp % _params.thread_number]._run ) {
            _condition_variable_parent.wait(_unique_lock, [&] { return threads[temp % _params.thread_number]._run ? false : true; } );
        }
        threads[temp % _params.thread_number]._available && _params.out_vcf << threads[temp % _params.thread_number].processed << std::endl;
        threads[temp % _params.thread_number]._available = false;
    }

    _params.verbose && std::cout << "Stoping all threads..." << std::endl;
    
    // Send stop signal to all threads and join them...
    for ( thread_index = 0; thread_index < _params.thread_number; thread_index++ ) {
        threads[thread_index]._stop = true;
        threads[thread_index]._condition_variable_thread.notify_one();
    }

    _params.verbose && std::cout << "Joining all threads..." << std::endl;

    // Join all the threads
    for ( thread_index = 0; thread_index < _params.thread_number; thread_index++ ) {
        threads[thread_index]._thread.join();
    }

    _params.verbose && std::cout << "All threads are detached." << std::endl;
    
    vcf.close();
    deinitialize_params(_params);

    _params.verbose && std::cout << "Total number of lines being processed: " << line_index << std::endl;
    
    return 0;
};



