//
//  main.c
//  Bioinformatics
//
//  Created by Akmuhammet Ashyralyyev on 17.03.2022.
//

#include <stdio.h>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include <zlib.h>
#include "htslib/faidx.h"

int main(int argc, const char * argv[]) {
    // Compulsory data structures & variables
    htsFile* bam_file; /* file pointer to the SAM/BAM/CRAM file */
    hts_idx_t* bam_file_index; /* pointer to the BAM/CRAM index */
    hts_itr_t *iter; /* iterator to read and parse the lines in SAM/BAM/CRAM file*/
    bam_hdr_t* bam_header; /* pointer to the header information of the BAM/CRAM file */
    bam1_core_t bam_alignment_core; /*  Structure for core alignment information */
    bam1_t* bam_alignment; /* Structure for one alignment */
    
    // Needs initialization
    char *path = "/home/akmuhammet/HG002.GRCh38.300x.bam"; /* input SAM/BAM/CRAM file name */
    char *fai_file = "/mnt/compgen/inhouse/share/rg_annot/b38/human_v38.fasta.fai"; /* file name for the reference genome FAI index as generated using,
                       e.g., samtools faidx ref.fa */
    
    // Initialization
    bam_file = hts_open( path, "r"); /* open SAM/BAM/CRAM file. htslib automatically detects the format of the file. Returns NULL if there is an error. */
    bam_header = sam_hdr_read(bam_file); /* Read in BAM header information */
    bam_file_index = sam_index_load( bam_file, path); /* load BAM/CRAM index */
    hts_set_fai_filename( bam_file, fai_file); /* pass the reference index to enable reading CRAM files */
    
    bam_alignment = bam_init1();
    
    int return_value = sam_read1( bam_file, bam_header, bam_alignment);
        /* return >= 0 on successfully reading a new record, -1 on end of stream, < -1 on error */
    
    int chrom_id = 1; /* ID of the chromosome as defined in the BAM/CRAM header. First chromosome ID is 0.*/
    int start = 100; /* starting coordinate: 0-based */
    int end = 1000; /* ending coordinate: 0-based */
    iter = sam_itr_queryi( bam_file_index, chrom_id, start, end);
        /* returns NULL if it fails */
    if (iter == NULL) {
        printf("Iter is NULL.\n");
    }

    while( sam_itr_next( bam_file, iter, bam_alignment) > 0) {
        bam_alignment_core = bam_alignment->core;
	printf(bam_alignment_core.pos);
	printf(" alignment \n");
    }

    if ( !(sam_itr_next( bam_file, iter, bam_alignment) > 0) ) {
        printf("Not valid iteration.\n");
    }
    // Cleanup
    sam_itr_destroy( iter);
    
    bam_destroy1( bam_alignment); /* free alignment pointer */
    bam_hdr_destroy( bam_header); /* free BAM/CRAM header pointer */
    hts_idx_destroy( bam_file_index); /* free BAM/CRAM index pointer */
    return_value = hts_close( bam_file); /* close the BAM/CRAM file */
    
    printf("End of the program...\n");
    return 0;
}
