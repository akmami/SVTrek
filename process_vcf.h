#ifndef __PROCESS_VCF_H__
#define __PROCESS_VCF_H__

#include "params.h"
#include "tpool.h"
#include "refinement.h"
#include <string.h>
#include <pthread.h>

/**
 * @brief Processes a VCF file using multithreading.
 *
 * This function reads a Variant Call Format (VCF) file and processes its entries 
 * using a thread pool. It creates multiple threads that retrieve and process 
 * VCF lines in parallel, improving efficiency for large datasets.
 *
 * @param params Pointer to the argument structure containing necessary parameters 
 *               for VCF processing.
 *
 * @return 0 on successful processing, or an error code if processing fails.
 */
int process_vcf(args *params);

#endif
