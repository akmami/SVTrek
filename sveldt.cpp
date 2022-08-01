#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <iostream>

using namespace std;
/**

    run ./sveldt /mnt/storage1/projects/giab/pacbio/HG002/HG002_PacBio_GRCh38.bam HG002_sv_summary.txt del > out.txt
 
    run ./sveldt --bam bam_file --vcf vcf_file --output output_file_name(default_output.txt)
*/

int main(int argc, char *argv[]){

    printf("# Program begins...\n");
    
    if (argc != 7 && argc != 5) {
        cout<<"Missing parameters. Please run the program as: "<<endl;
        cout<<"./sveldt --bam bam_file --vcf vcf_file --output output_file_name(default_output.txt)"<<endl;
        return -1;
    }

    if ( strcmp(argv[1], "--bam") != 0 || strcmp(argv[3], "--vcf") != 0) {
        cout << "Wrong input format. Please run the program as: " << endl;
        cout << "./sveldt --bam bam_file --vcf vcf_file --output output_file_name(default_output.txt)" << endl;
        return -1;
    }
    
    struct stat buffer_bam;
    if ( stat ( argv[2], &buffer_bam) != 0) {
        cout << "BAM file doesn't exist." << endl;
        return -1;
    }
    
    struct stat buffer_vcf;
    if ( stat ( argv[4], &buffer_vcf) != 0) {
        cout << "VCF file doesn't exist." << endl;
        return -1;
    }
    
    samFile *fp_in = hts_open(argv[1],"r");     //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);     //read header
    hts_idx_t *bam_file_index = sam_index_load( fp_in, argv[1] );
    bam1_t *aln = bam_init1();             //initialize an alignment

    FILE *sv_fp_in = fopen(argv[2], "r"); // open sv file
    char line[255];
    size_t len = 0;
    ssize_t read;
    if (sv_fp_in == NULL) {
        printf("# There is no such sv file.\n");
        return -1;
    }

    int comp;
    int chrom;
    char *id;
    int sv_pos, sv_end;
    char *alt;
    int outer_start, inner_start;
    int inner_end, outer_end;

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    // skip first line
    char *first_line = fgets(line, sizeof(line), sv_fp_in);
    printf("# ID\tCHROM\tPOS\tEND\tALT\tOUTER_START\tINNER_START\tINNER_END\tOUTER_END\n");

    while (fgets(line, sizeof(line), sv_fp_in) != NULL) {
        id = strtok(line, "\t");
        chrom = atoi(strtok(NULL, "\t")) - 1;
        sv_pos = atoi(strtok(NULL, "\t")) - 1;
        sv_end = atoi(strtok(NULL, "\t")) - 1;
        alt = strtok(NULL, "\t");
        outer_start = atoi(strtok(NULL, "\t")) + sv_pos - 1;
        inner_start = atoi(strtok(NULL, "\t")) + sv_pos - 1 ;
        inner_end = atoi(strtok(NULL, "\t")) + sv_end - 1;
        outer_end = atoi(strtok(NULL, "\t")) + sv_end - 1;
        
        printf("# %s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", id, chrom + 1, alt, sv_pos + 1, sv_end + 1, outer_start + 1, inner_start + 1, inner_end + 1, outer_end + 1);
        
        // ----------------------------------------------------------------------------------------------------------------------
        // Deletion
        // ----------------------------------------------------------------------------------------------------------------------
        
        if ( strcmp(alt, "del") == 0 ) {
            hts_itr_t *iter;
            iter = sam_itr_queryi( bam_file_index, chrom, outer_start - 40000, inner_start + 2000);
            if (iter == NULL) {
                printf("# invalid interval, iter is null\n");
            } else {
                while (sam_itr_next( fp_in, iter, aln ) > 0) {
                    pos = aln->core.pos;
                    uint32_t *cigar = bam_get_cigar(aln);

                    if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == 4) {
                        int reference_pos = pos;
                        for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                            if ( bam_cigar_op( cigar[i] ) != 1 && bam_cigar_op( cigar[i] ) != 4 && bam_cigar_op( cigar[i] ) != 5 && bam_cigar_op( cigar[i] ) != 6 ) {
                                reference_pos += bam_cigar_oplen( cigar[i] );
                            }
                        }
                        if ( outer_start <= reference_pos && reference_pos <= inner_start ) {
                            printf("%s\t%d\t%d\t-1\n", id, chrom + 1, reference_pos + 1);
                        }
                    }
                }
            }
            sam_itr_destroy(iter);
            
            iter = sam_itr_queryi( bam_file_index, chrom, inner_end - 2000, outer_end + 2000);
            if (iter == NULL) {
                printf("# invalid interval, iter is null\n");
            } else {
                while (sam_itr_next( fp_in, iter, aln ) > 0) {
                    pos = aln->core.pos;
                    uint32_t *cigar = bam_get_cigar(aln);

                    // Even if the read starts with soft clip, the pos is matched with aligned part
                    // Hence, some reads can be without soft clip and directly matched with the reference genome.
                    if ( bam_cigar_op( cigar[0] ) == 4 && inner_end <= pos && pos <= outer_end ) {
                        printf("%s\t%d\t%d\t1\n", id, chrom + 1, pos + 1);
                    }
                }
            }
            sam_itr_destroy(iter);
        }
        
        // ----------------------------------------------------------------------------------------------------------------------
        // Insertion
        // ----------------------------------------------------------------------------------------------------------------------
        
        if ( strcmp(alt, "ins") == 0 ) {
            hts_itr_t *iter;
            iter = sam_itr_queryi( bam_file_index, chrom, outer_start - 40000, inner_start + 2000);
            if (iter == NULL) {
                printf("# invalid interval, iter is null\n");
            } else {
                while (sam_itr_next( fp_in, iter, aln ) > 0) {
                    pos = aln->core.pos;
                    uint32_t *cigar = bam_get_cigar(aln);
                    
                    // If read ends in between (outer_start, inner_start)
                    if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == 4) {
                        int reference_pos = pos;
                        for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                            if ( bam_cigar_op( cigar[i] ) != 1 && bam_cigar_op( cigar[i] ) != 4 && bam_cigar_op( cigar[i] ) != 5 && bam_cigar_op( cigar[i] ) != 6 ) {
                                reference_pos += bam_cigar_oplen( cigar[i] );
                            }
                        }
                        if ( outer_start <= reference_pos && reference_pos <= inner_start ) {
                            printf("%s\t%d\t%d\t-1\n", id, chrom + 1, reference_pos + 1);
                        }
                    }
                }
            }
            sam_itr_destroy(iter);
            
            iter = sam_itr_queryi( bam_file_index, chrom, inner_end - 2000, outer_end + 2000);
            if (iter == NULL) {
                printf("# invalid interval, iter is null\n");
            } else {
                while (sam_itr_next( fp_in, iter, aln ) > 0) {
                    pos = aln->core.pos;
                    uint32_t *cigar = bam_get_cigar(aln);
                    
                    // If read starts in between (inner_end, outer_end)
                    if ( bam_cigar_op( cigar[0] ) == 4 && inner_end <= pos && pos <= outer_end  ) {
                        printf("%s\t%d\t%d\t1\n", id, chrom + 1, pos + 1);
                    }
                }
            }
            sam_itr_destroy(iter);
        }
        
        // ----------------------------------------------------------------------------------------------------------------------
        // Inversion
        // ----------------------------------------------------------------------------------------------------------------------
        
        if ( strcmp(alt, "inv") == 0 ) {
            hts_itr_t *iter;
            iter = sam_itr_queryi( bam_file_index, chrom - 1, outer_start - 40000, inner_start + 2000);
            if (iter == NULL) {
                printf("# invalid interval, iter is null\n");
            } else {
                while (sam_itr_next( fp_in, iter, aln ) > 0) {
                    pos = aln->core.pos;
                    uint32_t *cigar = bam_get_cigar(aln);
                    
                    // If read ends in between (outer_start, inner_start)
                    if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == 4) {
                        int reference_pos = pos;
                        for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                            if ( bam_cigar_op( cigar[i] ) != 1 && bam_cigar_op( cigar[i] ) != 4 && bam_cigar_op( cigar[i] ) != 5 && bam_cigar_op( cigar[i] ) != 6 ) {
                                reference_pos += bam_cigar_oplen( cigar[i] );
                            }
                        }
                        if ( outer_start <= reference_pos && reference_pos <= inner_start ) {
                            printf("%s\t%d\t%d\t-1\n", id, chrom + 1, reference_pos + 1);
                        }
                    }

                    // If read starts in between (outer_start, inner_start)
                    if ( bam_cigar_op( cigar[0] ) == 4 && outer_start <= pos && pos <= inner_start ) {
                        printf("%s\t%d\t%d\t-1\n", id, chrom + 1, pos + 1);
                    }
                }
            }
            sam_itr_destroy(iter);
            
            iter = sam_itr_queryi( bam_file_index, chrom, inner_end - 40000, outer_end + 2000);
            if (iter == NULL) {
                printf("# invalid interval, iter is null\n");
            } else {
                while (sam_itr_next( fp_in, iter, aln ) > 0) {
                    pos = aln->core.pos;
                    uint32_t *cigar = bam_get_cigar(aln);
                    
                    // If read ends in between (inner_end, outter_end)
                    if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == 4) {
                        int reference_pos = pos;
                        for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                            if ( bam_cigar_op( cigar[i] ) != 1 && bam_cigar_op( cigar[i] ) != 4 && bam_cigar_op( cigar[i] ) != 5 && bam_cigar_op( cigar[i] ) != 6 ) {
                                reference_pos += bam_cigar_oplen( cigar[i] );
                            }
                        }
                        if ( inner_end <= reference_pos && reference_pos <= outer_end ) {
                            printf("%s\t%d\t%d\t1\n", id, chrom + 1, reference_pos + 1);
                        }
                    }
                    
                    // If read starts in between (inner_end, outter_end)
                    if ( bam_cigar_op( cigar[0] ) == 4 && inner_end <= pos && pos <= outer_end ) {
                        printf("%s\t%d\t%d\t1\n", id, chrom + 1, pos + 1);
                    }
                }
            }
            sam_itr_destroy(iter);
        }
    }
    printf("# End of the program execution\n");

    fclose(sv_fp_in);
    bam_destroy1(aln);
    sam_close(fp_in);
    hts_idx_destroy(bam_file_index);
    bam_hdr_destroy(bamHdr);
    return 0;
}
