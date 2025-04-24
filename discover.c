#include "discover.h"


// KHASHL_MAP_INIT(SCOPE, HType, prefix, khkey_t, kh_val_t, __hash_fn, __hash_eq)
KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)


int parse_nodes(const char *readName, char *readPath, uint64_t **nodes, segment *segments, map32_t *h1) {

    uint64_t *temp_nodes = (uint64_t *)malloc(sizeof(uint64_t) * MAX_CIGAR);
    int node_size = 0, path_index, fwd_strand_count = 0, rev_strand_count = 0;
    const char *path = readPath; uint64_t id = 0; char strand = '>';
    int invalid_path = 0;
    
    while ((path_index = next_node(path, &id, &strand))) {
        path += path_index;
        khint_t k = map32_get(h1, id);
        if (k == kh_end(h1)) {
            fprintf(stderr, "[ERROR] Segment %lu in read %s not found.\n", id, readName); 
            invalid_path = 1;
            break;
        }
        temp_nodes[node_size++] = id;

        if (strand == '>') fwd_strand_count++;
        else rev_strand_count++;
        // validation 1
        segment *temp = segments + kh_val(h1, k);
        if (temp->rank > 1) {
            fprintf(stderr, "[ERROR] Read %s contains invalid ranks.\n", readName);
            invalid_path = 1;
            break;
        }
        // validation 2
        if (fwd_strand_count && rev_strand_count) {
            fprintf(stderr, "[ERROR] Read %s aligned in mixed strands.\n", readName);
            invalid_path = 1;
            break;
        }
    }
    // Continue if the wile loop terminated because of missing segment
    if (invalid_path) {
        free(temp_nodes);
        return 0;
    }

    *nodes = temp_nodes;
    return node_size;
}

int parse_gaf(const char* file_path, segment *segments, map32_t *h1) {

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;

    while ((line_len = io_read(file, &line, &line_cap))) {

        alignment aln;
        memset(&aln, 0, sizeof(alignment));

        char *saveptr, *token;
        int field = 0;

        token = strtok_r(line, "\t", &saveptr);
        while (token != NULL) {
            switch (field) {
                case 0: aln.readName = token; break;
                case 1: aln.readLen = atoi(token); break;
                case 2: aln.readStart = atoi(token); break;
                case 3: aln.readEnd = atoi(token); break;
                case 4: aln.strand = token[0]; break;
                case 5: aln.path = strdup(token); break;
                case 6: aln.pathLen = atoi(token); break;
                case 7: aln.pathStart = atoi(token); break;
                case 8: aln.pathEnd = atoi(token); break;
                case 9: aln.matches = atoi(token); break;
                case 10: aln.blockLen = atoi(token); break;
                case 11: aln.qual = atoi(token); break;
                default:
                    if (strncmp(token, "cg:Z:", 5) == 0) aln.cigar = strdup(token + 5);
                    break;
            }
            field++;
            token = strtok_r(NULL, "\t", &saveptr);
        }

        // Discard reads that has no quality
        if (aln.qual == 0) {
            free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // Get nodes and store in array. It is needed in case if it is mapped as reverse complement
        uint64_t *nodes;
        int node_size = parse_nodes(aln.readName, aln.path, &nodes, segments, h1);
        // Continue if the wile loop terminated because of missing segment
        if (!node_size) {
            free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // If read mapped as rc, reverse the order of the nodes
        if (aln.path[0] == '<') {
            reverse(nodes, node_size);
            {
                int new_start = 0, new_end = 0;
                fix_indices(aln.pathStart, aln.pathEnd, aln.pathLen, &new_start, &new_end);
                aln.pathStart = new_start;
                aln.pathEnd = new_end;
            }
            {
                int new_start = 0, new_end = 0;
                fix_indices(aln.readStart, aln.readEnd, aln.readLen, &new_start, &new_end);
                aln.pathStart = new_start;
                aln.pathEnd = new_end;
            }
        }

        // Initialize the cigar strings
        char ops[MAX_CIGAR]; int op_index = 0;

        // Soft clip for the prefix part of the read
        for (int i = 0; i < aln.readStart; i++) ops[op_index++] = __CIGAR_SOFT_CLIP;

        char cigar_ops[MAX_CIGAR];
        int n_ops = parse_cigar(aln.cigar, cigar_ops, MAX_CIGAR, aln.path[0] == '<');
        if (n_ops < 0) {
            fprintf(stderr, "[ERROR] Unable to parse CIGAR %s in read %s\n", aln.cigar, aln.readName);
            free(aln.path); free(nodes);
            if (aln.cigar) free(aln.cigar);
            continue;
        }

        // Initialize the segments
        int node_index = 0;
        khint_t k = map32_get(h1, nodes[node_index]);
        segment *seg = segments + kh_val(h1, k);
        segment *temp_prev_seg = (seg->rank == 0 ? seg : NULL);
        int p_length = strlen(seg->seq)-aln.pathStart;

        int reference_start = (seg->rank == 0 ? seg->start+aln.pathStart+1 : -1);
        int cigar_op_index = 0, is_reference_start_set = 0;

        while (cigar_op_index < n_ops) {
            char op = cigar_ops[cigar_op_index++];

            if (!is_reference_start_set && seg->rank == 0 && CIGAR_REF(op)) {
                if (op == __CIGAR_SEQUENCE_MATCH) is_reference_start_set = 1;
                else reference_start++;
            }

            if (seg->rank == 0) ops[op_index++] = op;
            else if (CIGAR_QUE(op)) ops[op_index++] = __CIGAR_INSERTION;

            if (CIGAR_REF(op)) {
                p_length--;

                if (p_length) continue;

                node_index++;
                if (node_index == node_size) break;
                khint_t k = map32_get(h1, nodes[node_index]);
                seg = segments + kh_val(h1, k);
                p_length = strlen(seg->seq);
                if (seg->rank == 0) {
                    if(!is_reference_start_set) reference_start = seg->start;

                    if (temp_prev_seg != NULL) {
                        for (int index = temp_prev_seg->end; index < seg->start; index++)
                            ops[op_index++] = __CIGAR_DELETION;
                    }

                    temp_prev_seg = seg;
                }
            }
        }

        // Soft clip for the end part of the read
        for (int i = aln.readEnd; i < aln.readLen; i++) ops[op_index++] = __CIGAR_SOFT_CLIP;

        // Validate cigar
        int cigar_query_count = 0;
        for (int i=0; i<op_index; i++) {
            if (CIGAR_QUE(ops[i])) cigar_query_count++;
        }
        if (cigar_query_count != aln.readLen) 
            fprintf(stderr, "[ERROR] CIGAR (query) mismatch in length. %d %d\n", cigar_query_count, aln.readLen);
        
        // Cleanup
        free(aln.path); free(nodes);
        if (aln.cigar) free(aln.cigar);
    }    
    fclose(file);
    free(line);

    return 0;
}

int parse_gfa(const char* file_path, segment **segments, int *size, map32_t *h1) {

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;
    
    int segment_size = 0, segment_capacity = 1000000;
    segment *temp_segments = (segment *)malloc(segment_capacity * sizeof(segment));

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'S') {
            if (segment_size == segment_capacity) {
                segment_capacity = segment_capacity * 2;
                segment *temp = realloc(temp_segments, segment_capacity * sizeof(segment));
                if (temp == NULL) {
                    fprintf(stderr, "[ERROR] Realloc failed.\n");
                    free(line); free_segments(&temp_segments, segment_size); map32_destroy(h1);
                    return 1;
                }
                temp_segments = temp;
            }
            
            char *saveptr;
            char *token = strtok_r(line, "\t", &saveptr); // Skip 'S'
            // ID
            token = strtok_r(NULL, "\t", &saveptr); // ID
            temp_segments[segment_size].id = strtoull(token, NULL, 10);
            
            // Sequence
            token = strtok_r(NULL, "\t", &saveptr); // Sequence string
            temp_segments[segment_size].seq = malloc(strlen(token) + 1);
            if (!temp_segments[segment_size].seq) {
                fprintf(stderr, "[ERROR] Memory allocation failed\n");
                return 1;
            }
            strcpy(temp_segments[segment_size].seq, token);
            temp_segments[segment_size].rank = 1;

            // INFO
            token = strtok_r(NULL, "\t", &saveptr); // process rest 
            while (token) {
                // process SN:Z:chr22	SO:i:10510000	SR:i:0
                char *info_saveptr;
                char *tag = strtok_r(token, ":", &info_saveptr);  // get TAG
                char *type = strtok_r(NULL, ":", &info_saveptr);  // get TYPE
                char *value = strtok_r(NULL, ":", &info_saveptr); // get VALUE

                if (tag && type && value) {
                    if (strcmp(tag, "SN") == 0 && strcmp(type, "Z") == 0) {
                        // temp_segments[segment_size].seq_name = strdup(value);
                    } else if (strcmp(tag, "SO") == 0 && strcmp(type, "i") == 0) {
                        temp_segments[segment_size].start = atoi(value);
                        temp_segments[segment_size].end = atoi(value) + strlen(temp_segments[segment_size].seq);
                    } else if (strcmp(tag, "SR") == 0 && strcmp(type, "i") == 0) {
                        temp_segments[segment_size].rank = atoi(value);
                    } 
                }

                token = strtok_r(NULL, "\t", &saveptr); // move to the next field
            }
            
            khint_t k; int absent;
            k = map32_put(h1, temp_segments[segment_size].id, &absent);
            kh_val(h1, k) = segment_size; // set value as index in segments
            segment_size++;
        } else if (line[0] == 'L') {
            uint64_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);

            if (overlap) {
                fprintf(stderr, "[ERROR] Overlaps are not zero, cannot make conversion.\n");
                return 1;
            }
        } else if (line[0] == 'P') {
            char *token = strtok(line, "\t");
            token = strtok(NULL, "\t"); // skip path ID
            
            token = strtok(NULL, "\t"); // process segment IDs

            char *segment_token = strtok(token, ",");
            segment *prev_segment = NULL;
            int ref_pos = 0;

            while (segment_token) {
                size_t len = strlen(segment_token);
                if (segment_token[len - 1] == '+') segment_token[len - 1] = '\0'; // remove the trailing '+'

                uint64_t seg_id = strtoull(segment_token, NULL, 10);
                khint_t k = map32_get(h1, seg_id);
                segment *current_segment = temp_segments + kh_val(h1, k);

                if (prev_segment) prev_segment->rank = 0;

                current_segment->start = ref_pos;
                ref_pos += strlen(current_segment->seq);
                current_segment->end = ref_pos; // for now, not needed.
                prev_segment = current_segment; // update previous
                segment_token = strtok(NULL, ","); // next segment
            }
            if (prev_segment) prev_segment->rank = 0;
        }
    }
    io_close(file, &line);

    *segments = temp_segments;
    *size = segment_size;

    return 0;
}

int discover(int argc, char *argv[]) {
    
    disc_args params;
    init_disc(argc, argv, &params);

    map32_t *h1 = map32_init();
    segment *segments;
    int segment_size;
    if (parse_gfa(params.gfa_file, &segments, &segment_size, h1)) {
        fprintf(stderr, "[ERROR] GFA file parsing resulted in null.\n");
        return 0;
    }

    if(parse_gaf(params.gaf_file, segments, h1)) {
        fprintf(stderr, "[ERROR] GAF file parsing failed.\n");
        return 0;
    }

    for (int i = 0; i < segment_size; i++) {
        if (segments[i].seq) free(segments[i].seq);
    }
    free(segments);

    map32_destroy(h1);


    return 1;
}