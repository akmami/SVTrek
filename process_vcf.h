#ifndef __PROCESS_VCF_H__
#define __PROCESS_VCF_H__

#include "params.h"
#include "refinement.h"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cstring>  // for strncmp

int process_vcf(struct params &_params);

#endif // __PROCESS_VCF_H__
