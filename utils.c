#include "utils.h"
#include "refinement.h"  

void quicksort(int *array, int low, int high) {
    if (low < high) {
        int pivot = array[high];
        int i = low - 1;

        for (int j = low; j < high; j++) {
            if (array[j] < pivot) {
                i++;
                int temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }

        int temp = array[i+1];
        array[i+1] = array[high];
        array[high] = temp;

        quicksort(array, low, i);
        quicksort(array, i + 2, high);
    }
}

void reverse(uint64_t *arr, int size) {
    int i = 0;
    while (i < size / 2) {
        uint64_t temp = arr[i];
        arr[i] = arr[size-i-1];
        arr[size-i-1] = temp;
        i++;
    }
}

void fix_indices(int start, int end, int len, int *new_start, int *new_end) {
    int preLen = start;
    int postLen = len - end;
    int total_len = len - preLen - postLen;
    *new_start = postLen;
    *new_end = postLen + total_len;
}

int parse_cigar(char *cigar, char *ops, int max_ops, int rev) {
    int i = 0, j = 0, num = 0;
    while (cigar[i]) {
        if ('0' <= cigar[i] && cigar[i] <= '9') {
            num = num * 10 + (cigar[i] - '0');
        } else {
            if (j + num >= max_ops) return -1; // out of range
            for (int k = 0; k < num; k++) {
                ops[j++] = cigar[i];
            }
            num = 0;
        }
        i++;
    }
    
    if (rev) {
        int k=0;
        while (k < j/2) {
            char temp = ops[k];
            ops[k] = ops[j-k-1];
            ops[j-k-1] = temp;
            k++;
        }
    }
    return j;  // number of ops
}

int next_node(const char *path, uint64_t *id, char *strand) {
    if (*path != '\0') {
        *strand = *path;
        path++;
        int index = 1; uint64_t i = 0;
        while (*path >= '0' && *path <= '9') {
            i = i * 10 + (*path - '0');
            path++; index++;
        }
        *id = i;
        return index;
    }
    return 0;
}

FILE *io_open(const char* file_path, char **line, int cap) {

    FILE* file = fopen(file_path, "r");
    if (!file) { fprintf(stderr, "[ERROR] Failed to open file %s\n", file_path); exit(EXIT_FAILURE); }
    
    *line = (char *)malloc(cap);

    if (!(*line)) { fprintf(stderr, "[ERROR] Failed to allocate memory."); exit(EXIT_FAILURE); }

    return file;
}

void io_close(FILE *file, char **str) {
    if (file)
        fclose(file);
    if (*str) {
        free(*str);
    }
    file = NULL;
}

int io_read(FILE *file, char **str, size_t *cap) {

    char *line = *str;
    size_t line_len = *cap;

    if (fgets(line, line_len, file)) {
        size_t len = strlen(line);

        // read the line and fit it to `line`
        while (len == line_len - 1 && line[len - 1] != '\n') {
            line_len *= 2;
            *cap = line_len;
            char *temp_line = (char *)realloc(line, line_len);
            if (!temp_line) { fprintf(stderr, "[ERROR] Memory reallocation failed.\n"); return 0; }
            line = temp_line;
            *str = temp_line;

            if (fgets(line + len, line_len - len, file) == NULL) { break; }
            len = strlen(line);
        }
        if (len && line[len - 1] == '\n') { line[len - 1] = '\0'; len--; }

        return len;
    }

    return 0;
}

void free_segments(segment **segments, int segment_size) {
    segment *temp = *segments;
    for (int i = 0; i < segment_size; i++) {
        free(temp[i].seq);
    }
    free(*segments);
}

int consensus2(FILE *fp, t_arg *params, uint32_t imprecise_pos) {
    if (!fp) {
        fprintf(stderr, "Invalid file pointer\n");
        return -1;
    }

    uint32_t arr[ROW_SIZE];
    size_t items_read = fread(arr, sizeof(uint32_t), ROW_SIZE, fp);
    if (items_read != ROW_SIZE) {
        if (feof(fp)) {
            return -1; 
        } else {
            fprintf(stderr, "Error reading file: %s\n", strerror(errno));
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
