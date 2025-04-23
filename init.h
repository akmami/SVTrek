#ifndef __INIT_C__
#define __INIT_C__

#include "params.h"
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

/**
 * @brief Print usage manual to the console.asm
 */
void printUsage();

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
void init_audt(int argc, char *argv[], audt_args *params);


void init_disc(int argc, char *argv[], disc_args *params);

#endif