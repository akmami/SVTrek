#include "init.h"
#include "params.h"
#include "process_vcf.h"

int main(int argc, char *argv[]){

    args params;

    init(argc, argv, &params);

    if (params.verbose) {
        printf("Program begins...\n");
    }
    
    process_vcf(&params);

    if (params.verbose) {
        printf("End of the program execution\n");
    }

    return 0;
};
