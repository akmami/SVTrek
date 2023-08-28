#ifndef THREAD_DATA_CPP
#define THREAD_DATA_CPP

#include <thread>
#include <condition_variable>
#include <mutex>


struct thread_data {
    int thread_id;
    std::thread _thread;
    //std::packaged_task<std::string()> _task;
    std::condition_variable _condition_variable_thread;
    std::mutex _mutex; // Mutex used for avoiding data races
    htslib_params htslib;
    bool _stop = false;
    bool _run = false;
    bool _available = false;
    std::string line;
    std::string processed;
};


#endif