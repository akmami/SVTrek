#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "params.h"      
#include "refinement.h"  
#include "utils.h"
#define ROW_SIZE 1024
#define TYPE_IND 0xE0000000U   
#define INDEX_IND  0x1FFFFFFFU  


int process_row(FILE *fp, t_arg *params, uint32_t imprecise_pos) {
    uint32_t arr[ROW_SIZE];
    size_t read_result = fread(arr, sizeof(arr), 1, fp);
    if (read_result == 0) {
        if (feof(fp)) {
            return -1; 
        } else {
            fprintf(stderr, "Error reading file");
            return -1;
        }
    }

    quicksort((int*)arr, 0, ROW_SIZE-1);

    int ins_cap = 128, ins_cnt = 0;
    int del_cap = 128, del_cnt = 0;
    int mis_cap = 128, mis_cnt = 0;
    
    uint32_t *ins_temp = malloc(ins_cap * sizeof(uint32_t));
    uint32_t *del_temp = malloc(del_cap * sizeof(uint32_t));
    uint32_t *mis_temp = malloc(mis_cap * sizeof(uint32_t));
    
    if (!ins_temp || !del_temp || !mis_temp) {
        fprintf(stderr, "Memory allocation failed\n");
        free(ins_temp);
        free(del_temp);
        free(mis_temp);
        return -1;
    }

    for (int i = 0; i < ROW_SIZE; i++) {
        uint32_t t = (arr[i] &TYPE_IND) >> 29;
        uint32_t idx = arr[i] &INDEX_IND;
        
        if (t == 1) {  
            if (ins_cnt == ins_cap) {
                ins_cap *= 2;
                uint32_t *new_temp = realloc(ins_temp, ins_cap * sizeof(uint32_t));
                if (!new_temp) {
                    fprintf(stderr, "Memory reallocation failed for insertion\n");
                    free(ins_temp);
                    free(del_temp);
                    free(mis_temp);
                    return -1;
                }
                ins_temp = new_temp;
            }
            ins_temp[ins_cnt++] = idx;
        }
        else if (t == 2) {  // deletion (002) will be implemented later
            
        }
        else if (t == 3) {  // mismatch (003) will be implemented later
            
        }
    }
    //memory allocation controls may be deleted in the future developments as they maybe unnecessary   
    int cons = -1;
    
    if (ins_cnt > 0) {
        cons = consensus_pos((int*)ins_temp, ins_cnt, imprecise_pos, params->consensus_min_count,
                            params->consensus_interval, params->consensus_interval_range);
        
        if (cons < 0) {
            fprintf(stderr, "Insertion consensus failed\n");
        } else {
            printf("Row consensus INS: %d\n", cons);
        }
    }
    
    if (del_cnt > 0) {
        //Deletion consensus 
    }
    
    if (mis_cnt > 0) {
        //Mismatch consensus
    }

    free(ins_temp);
    free(del_temp);
    free(mis_temp);
    
    return cons;
}
