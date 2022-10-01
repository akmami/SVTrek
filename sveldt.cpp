#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>

using namespace std;

/**

    run ./sveldt /mnt/storage1/projects/giab/pacbio/HG002/HG002_PacBio_GRCh38.bam HG002_sv_summary.txt del > out.txt
 
    run ./sveldt --bam bam_file --vcf vcf_file --output output_file_name(default_output.txt)
*/

// MARK: - helper function declaration

int run(char argv2[], char argv4[]);
bool process_line(string line);
bool insertion(char *id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end);
bool deletion(char *id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end);
bool inversion(char *id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end);

// MARK: - global variables

samFile *fp_in;
bam_hdr_t *bamHdr;
hts_idx_t *bam_file_index;
bam1_t *aln;

// MARK: - main program

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
    
    run(argv[2], argv[4]);
    
    return 0;
}

// MARK: - implementation of helper functions

int run(char argv2[], char argv4[]) {
    // init global variables
    fp_in = hts_open(argv2,"r");            //open bam file
    bamHdr = sam_hdr_read(fp_in);           //read header
    bam_file_index = sam_index_load( fp_in, argv2 );
    aln = bam_init1();                      //initialize an alignment
    
    ifstream vcf_file(argv4);                   // open sv file
    string line;
    
    while (getline(vcf_file, line)) {
        if (line.at(0) == '#')
            continue;
        process_line(line);
    }
    
    cout << "# End of the program execution" << endl;
    
    vcf_file.close();
    
    // de-init global variables, prevent memory leaks
    bam_destroy1(aln);
    sam_close(fp_in);
    hts_idx_destroy(bam_file_index);
    bam_hdr_destroy(bamHdr);
    
    return 0;
}

bool process_line(string line) {
    
    // extract data
    if (line.find("IMPRECISE") == string::npos || line.find("SVLEN") != string::npos || line.find("CIPOS") != string::npos) {
        
        // local variables
        int comp, chrom, sv_pos;
        string id, alt;
        
        vector<string> splitted_line;
        boost::split(splitted_line, line, boost::is_any_of("\t"));
        try {
            chrom = stoi(splitted_line[0]);
        } catch (...) {
            try {
                chrom = stoi(splitted_line[0].substr(3));
            } catch (...) {
                // unlucky :(
                if ( (splitted_line[0].length() && (splitted_line[0][0] == 'X' || splitted_line[0][0] == 'x') ) ||
                        splitted_line[0].length() == 4 && (splitted_line[0][3] == 'X' || splitted_line[0][3] == 'x') ) {
                    // chrom is X
                } else if ( ( splitted_line[0].length() && (splitted_line[0][0] == 'Y' || splitted_line[0][0] == 'y') ) ||
                           ( splitted_line[0].length() == 4 && (splitted_line[0][3] == 'Y' || splitted_line[0][3] == 'y') ) ){
                    // chrom is Y
                } else {
                    cout << "Chromosome could not be cinverted into integer. CHROM tag : " << splitted_line[0] << endl;
                }
                return false;
            }
        }
        try {
            sv_pos = stoi(splitted_line[1]);      // it is set but never used
        } catch (...) {
            sv_pos = -1;
            return false;
        }
        
        id = splitted_line[2];
        alt = splitted_line[4];
        
        // CIPOS and CIEND info retrieval part
        int outer_start, inner_start, inner_end, outer_end;
        
        
        
        // call function according to sv type
        int type_index = splitted_line[7].find("SVTYPE");
        
        if (type_index == -1) {
            return false;
        }
        
        if ( splitted_line[7].find("INS", type_index + 7) != -1 || splitted_line[7].find("ins", type_index + 7) != -1 ) {
            insertion(id, chrom, outer_start, inner_start, inner_end, outer_end);
        } else if ( splitted_line[7].find("DEL", type_index + 7) != -1 || splitted_line[7].find("del", type_index + 7) != -1 ) {
            deletion(id, chrom, outer_start, inner_start, inner_end, outer_end);
        } else if ( splitted_line[7].find("INV", type_index + 7) != -1 || splitted_line[7].find("inv", type_index + 7) != -1 ) {
            inversion(id, chrom, outer_start, inner_start, inner_end, outer_end);
        } else {
            // Other type of sv
            return false;
        }
    }
    
    return true;
}

bool insertion(char *id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end) {

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
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
    
    return true;
}


bool deletion(char *id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end) {

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
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
    
    return true;
}

bool inversion(char *id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end) {

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
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
    
    return true;
}
