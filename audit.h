#ifndef __AUDIT_H__
#define __AUDIT_H__

#include "params.h"
#include "init.h"
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
 */
int audit(int argc, char *argv[]);

#endif
