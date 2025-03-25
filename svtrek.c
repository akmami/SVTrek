#include "init.h"
#include "params.h"
#include "process_vcf.h"

int main(int argc, char *argv[]){

    args params;

    init(argc, argv, &params);

    if (params.verbose) {
        printf("Program begins...\n");
        if (params.mode == MODE_EVAL)
            printf("Mode: Evaluation\n");
        else if (params.mode == MODE_DISC)
            printf("Mode: Discovery\n");
    }
    
    process_vcf(&params);

    if (params.verbose) {
        printf("End of the program execution\n");
    }

    return 0;
};
