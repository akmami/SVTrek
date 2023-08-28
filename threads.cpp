#ifndef THREADS_CPP
#define THREADS_CPP

#include "htslib.h"
#include "params.h"
#include "thread_data.h"
#include "process_line.cpp"


void thread_func(params* _params, thread_data* _thread_data, std::condition_variable* _condition_variable_parent) {
    
    _thread_data->htslib.fp_in = hts_open(_params->bam_file, "r");         //open bam file
    _thread_data->htslib.bam_hdr = sam_hdr_read(_thread_data->htslib.fp_in);    //read header
    _thread_data->htslib.bam_file_index = sam_index_load( _thread_data->htslib.fp_in, _params->bam_file );

    std::unique_lock<std::mutex> _unique_lock(_thread_data->_mutex, std::defer_lock);
    while (true) {

        // Wait signal (either stop or new task)
        _thread_data->_condition_variable_thread.wait(_unique_lock, [_thread_data] () { return _thread_data->_stop || _thread_data->_run; });
        if ( !(_thread_data->_run) && _thread_data->_stop) { 
            goto exit_thread; 
        }

        // Execute the task!
        _thread_data->processed = process_line(*_params, _thread_data->line, *_thread_data);

        // Make thread available to next task
        _thread_data->_run = false;
        _thread_data->_available = true;
        _condition_variable_parent->notify_one();
    }

    exit_thread:

    sam_close(_thread_data->htslib.fp_in);
    hts_idx_destroy(_thread_data->htslib.bam_file_index);
    bam_hdr_destroy(_thread_data->htslib.bam_hdr);

    _thread_data->htslib.fp_in = NULL;
    _thread_data->htslib.bam_hdr = NULL;
    _thread_data->htslib.bam_file_index = NULL;
};


#endif