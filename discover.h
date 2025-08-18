#ifndef __DISCOVER_H__
#define __DISCOVER_H__

#include "params.h"
#include "init.h"
#include "utils.h"
#include "khashl.h"
#include <stdint.h>
#include <htslib/kseq.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>

#ifdef DEBUG
#define READ_MODE "rb"
#define WRITE_MODE "wb"
#else
#define READ_MODE "r"
#define WRITE_MODE "w"
#endif

int discover(int argc, char *argv[]);

#endif