#include "process_vcf.h"

// Example of splitting a line on a delimiter (e.g., '\t')
std::vector<std::string> split_line(const std::string &line, char delimiter) {
    std::vector<std::string> fields;
    size_t start = 0, endPos;
    while ((endPos = line.find(delimiter, start)) != std::string::npos) {
        fields.push_back(line.substr(start, endPos - start));
        start = endPos + 1;
    }
    // Last field
    fields.push_back(line.substr(start));
    return fields;
}

std::string &process_line(params &_params, std::string &line, thread_data &_thread_data) {
    // Only proceed if line has CIPOS and SVTYPE (for instance).
    // You can adapt this condition to your logic.
    if (line.find("CIPOS") != std::string::npos && line.find("SVTYPE=") != std::string::npos) {
        std::istringstream iss(line);
        std::string _chrom, _pos, _id, _ref, _alt, _qual, _filter, _info, _format;

        // Read each tab-delimited column
        std::getline(iss, _chrom, '\t');
        std::getline(iss, _pos, '\t');
        std::getline(iss, _id, '\t');
        std::getline(iss, _ref, '\t');
        std::getline(iss, _alt, '\t');
        std::getline(iss, _qual, '\t');
        std::getline(iss, _filter, '\t');
        std::getline(iss, _info, '\t');
        iss >> _format;  // read whatever remains

        // Variables to hold numeric conversions
        int chrom, pos, end;
        // CIPOS and CIEND produce "outer" and "inner" coordinates
        int outer_start = 0, inner_start = 0;
        int outer_end = 0,   inner_end = 0;
        std::string sv_type;

        // 1. Parse CHROM
        try {
            chrom = std::stoi(_chrom);
        } catch (...) {
            // Fallback if e.g. "chr5"
            try {
                chrom = std::stoi(_chrom.substr(3));  // skip "chr"
            } catch (...) {
                return line; // can't parse => just return unmodified
            }
        }

        // 2. Parse POS
        try {
            pos = std::stoi(_pos);
        } catch (...) {
            return line; // can't parse => just return unmodified
        }

        // 3. Parse SVTYPE=
        size_t index_1 = _info.find("SVTYPE=");
        if (index_1 == std::string::npos) {
            return line;
        }
        size_t index_2 = _info.find(";", index_1);
        sv_type = _info.substr(index_1 + 7, index_2 - (index_1 + 7));

        // 4. Parse END=
        index_1 = _info.find("END=");
        if (index_1 == std::string::npos) {
            return line;
        }
        index_2 = _info.find(";", index_1);
        try {
            end = std::stoi(_info.substr(index_1 + 4, index_2 - (index_1 + 4)));
        } catch (...) {
            return line;
        }

        // 5. Parse CIPOS= => outer_start, inner_start
        index_1 = _info.find("CIPOS=");
        if (index_1 != std::string::npos) {
            index_2 = _info.find(";", index_1 + 6);
            if (index_2 == std::string::npos)
                index_2 = _info.size();

            std::string cipos_sub = _info.substr(index_1 + 6,
                                                 index_2 - (index_1 + 6));
            std::istringstream iss_cipos(cipos_sub);
            std::string left_str, right_str;
            std::getline(iss_cipos, left_str, ',');
            std::getline(iss_cipos, right_str, ',');

            try {
                outer_start = std::stoi(left_str) + pos;
                inner_start = std::stoi(right_str) + pos;
            } catch (...) {
                outer_start = pos;
                inner_start = pos;
            }
        } else {
            // Fallback if CIPOS not found
            outer_start = pos;
            inner_start = pos;
        }

        // 6. Parse CIEND= => outer_end, inner_end
        index_1 = _info.find("CIEND=");
        if (index_1 != std::string::npos) {
            index_2 = _info.find(";", index_1 + 6);
            if (index_2 == std::string::npos)
                index_2 = _info.size();

            std::string ciend_sub = _info.substr(index_1 + 6,
                                                 index_2 - (index_1 + 6));
            std::istringstream iss_ciend(ciend_sub);
            std::string left_str, right_str;
            std::getline(iss_ciend, left_str, ',');
            std::getline(iss_ciend, right_str, ',');

            try {
                outer_end = std::stoi(left_str) + end;
                inner_end = std::stoi(right_str) + end;
            } catch (...) {
                outer_end = end;
                inner_end = end;
            }
        } else {
            // Fallback if CIEND not found
            outer_end = end;
            inner_end = end;
        }

        // 7. Call the appropriate variation function from 'variations.h'
        if (sv_type == "INS" || sv_type == "INS:ME") {
            result res = insertion(chrom, outer_start,
                                           inner_start,
                                           pos,  // imprecise_pos
                                           _params,
                                           _thread_data);
            if (res.refined_start != -1) {
                _pos = std::to_string(res.refined_start);
                _info += ";SVELDT=SUCCESS";
            } else {
                _info += ";SVELDT=INCORRECT";
            }
        }
        else if (sv_type == "DEL" || sv_type == "DEL:ME") {
            result res = deletion(chrom,
                                          outer_start,
                                          inner_start,
                                          inner_end,
                                          outer_end,
                                          pos,  // imprecise_pos
                                          end,  // imprecise_end
                                          _params,
                                          _thread_data);
            if (res.refined_start != -1) {
                _pos = std::to_string(res.refined_start);
                _info += ";SVELDT=SUCCESS";
            } else {
                _info += ";SVELDT=INCORRECT";
            }
        }
        else if (sv_type == "INV") {
            result res = inversion(chrom,
                                           outer_start,
                                           inner_start,
                                           inner_end,
                                           outer_end,
                                           pos,  // imprecise_pos
                                           end,  // imprecise_end
                                           _params,
                                           _thread_data);
            if (res.refined_start != -1) {
                _pos = std::to_string(res.refined_start);
                _info += ";SVELDT=SUCCESS";
            } else {
                _info += ";SVELDT=INCORRECT";
            }
        }

        // Reconstruct the line with possibly updated POS and INFO
        line = _chrom + "\t" + _pos + "\t" + _id + "\t"
               + _ref + "\t" + _alt + "\t" + _qual + "\t"
               + _filter + "\t" + _info + "\t" + _format;
    }

    // Return the (possibly) modified line
    return line;
}

void thread_func(struct params* _params, struct thread_data* _thread_data, std::condition_variable* _condition_variable_parent) {
    _thread_data->htslib.fp_in    = hts_open(_params->bam_file, "r");
    _thread_data->htslib.bam_hdr  = sam_hdr_read(_thread_data->htslib.fp_in);
    _thread_data->htslib.bam_file_index = sam_index_load(_thread_data->htslib.fp_in, _params->bam_file);

    std::unique_lock<std::mutex> _unique_lock(_thread_data->_mutex, std::defer_lock);
    while (true) {
        _thread_data->_condition_variable_thread.wait(_unique_lock, [_thread_data]() {
            return _thread_data->_stop || _thread_data->_run;
        });

        if (!_thread_data->_run && _thread_data->_stop) {
            break;
        }
        _thread_data->processed = process_line(*_params, _thread_data->line, *_thread_data);

        _thread_data->_run = false;
        _thread_data->_available = true;
        _condition_variable_parent->notify_one();
    }

    sam_close(_thread_data->htslib.fp_in);
    hts_idx_destroy(_thread_data->htslib.bam_file_index);
    bam_hdr_destroy(_thread_data->htslib.bam_hdr);

    _thread_data->htslib.fp_in        = nullptr;
    _thread_data->htslib.bam_hdr      = nullptr;
    _thread_data->htslib.bam_file_index = nullptr;
}

int process_vcf(params &_params) {

    initialize_params(_params);

    std::ifstream vcf(_params.vcf_file);                // open sv file
    std::string line;

    // Multithreading
    _params.verbose && std::cout << "Initializing threads..." << std::endl;

    thread_data threads[_params.thread_number];
    int thread_index;

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
    bool description = false;

    while (getline(vcf, line)) {

        line_index++;

        if (line.at(0) == '#') {
            if ( !strncmp(line.c_str(), "##INFO", 6 ) ) {
                if ( !description ) {
                    _params.out_vcf << "##INFO=<ID=SVELDT,Number=1,Type=String,Description=\"The SV is tagged by SVELDT program:SIMULATED=The SV is only simulated var varsim.py and not processed by sveldt yet, SUCCESS=SVELDT was able to refine all given intervals, PARTIAL=SVELDT was able to refine only one of the points, INCORRECT=SVELDT detected invalid SV.\"" << std::endl;
                    description = true;
                }
            }

            if ( !strncmp(line.c_str(), "#CHROM", 6 ) ) {
                if ( !description ) {
                    _params.out_vcf << "##INFO=<ID=SVELDT,Number=1,Type=String,Description=\"The SV is tagged by SVELDT program:SIMULATED=The SV is only simulated var varsim.py and not processed by sveldt yet, SUCCESS=SVELDT was able to refine all given intervals, PARTIAL=SVELDT was able to refine only one of the points, INCORRECT=SVELDT detected invalid SV.\"" << std::endl;
                    description = true;
                }
            }

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