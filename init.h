#ifndef __INIT_C__
#define __INIT_C__

#include "params.h"
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

/**
 * @brief Parses command-line arguments and initializes parameters.
 *
 * This function processes command-line arguments and initializes the `args` structure 
 * with the necessary parameters for execution.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param params Pointer to the structure where parsed arguments will be stored.
 */
void init(int argc, char *argv[], args *params);

#endif