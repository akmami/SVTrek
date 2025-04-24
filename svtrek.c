#include "discover.h"
#include "audit.h"


int main(int argc, char *argv[]){

    if (argc < 2) {
        printUsage();
        exit(1);
    }

    if (strcmp(argv[1], "disc") == 0) {
        discover(argc, argv);
    } else if (strcmp(argv[1], "audt") == 0) {
        audit(argc, argv);
    } else {
        printUsage();
        exit(1);
    }

    return 0;
};
