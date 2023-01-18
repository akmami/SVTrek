#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <cmath>
#include "argh.h"

//                                Consumes query   reference   Op
#define CIGAR_ALIGNMENT_MATCH       0   //  yes     yes         M
#define CIGAR_INSERTION             1   //  yes     no          I
#define CIGAR_DELETION              2   //  no      yes         D
#define CIGAR_SKIPPED               3   //  no      yes         N
#define CIGAR_SOFT_CLIP             4   //  yes     no          S
#define CIGAR_HARD_CLIP             5   //  no      no          H
#define CIGAR_PADDING               6   //  no      no          P
#define CIGAR_SEQUENCE_MATCH        7   //  yes     yes         =
#define CIGAR_SEQUENCE_MISMATCH     8   //  yes     yes         X

#define FLAG_MULTIPLE_SEGMENTS          0x1
#define FLAG_SECONDARY_ALIGNMENT        0x100
#define FLAG_SUPPLEMENTARY_ALIGNMENT    0x800


#define WIDER_INTERVAL          40000
#define NARROW_INTERVAL         2000
#define CONSENSUS_INTEVAL       3
#define CONSENSUS_MIN_COUNT     2
//#define CONSENSUS_COUNT_PERC    0.3

using namespace std;

// MARK: - helper function declaration

int run(char const *bam, char const *vcf);
bool process_line(string line);
bool insertion(string id, int chrom, int outer_start, int inner_start);
bool deletion(string id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end);
bool inversion(string id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end);
int find_start(int chrom, int outer_start, int inner_start);
int find_end(int chrom, int inner_end, int outer_end);
int find_start_or_end(int chrom, int start, int end);
int consensus(vector<int> locations);

// MARK: - global variables

samFile *fp_in;
bam_hdr_t *bamHdr;
hts_idx_t *bam_file_index;
bam1_t *aln;

// MARK: - main program

int main(int argc, char *argv[]){

    printf("# Program begins...\n");
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);
    string bam_str, vcf_str, output_str;
    
    if ( !( cmdl({"--bam", "b"}) >> bam_str ) ) {
        cout << "BAM file is missing." << endl;
        exit(-1);
    }
    
    if ( !( cmdl({"--vcf", "v"}) >> vcf_str ) ) {
        cout << "VCF file is missing." << endl;
        exit(-1);
    }
    
    
    const char* bam = bam_str.c_str();
    const char* vcf = vcf_str.c_str();
    
    struct stat buffer_bam;
    if ( stat ( bam, &buffer_bam) != 0) {
        cout << "BAM file doesn't exist." << endl;
        return -1;
    }
    
    struct stat buffer_vcf;
    if ( stat ( vcf, &buffer_vcf) != 0) {
        cout << "VCF file doesn't exist." << endl;
        return -1;
    }
    
    run(bam, vcf);
    
    return 0;
}

// MARK: - implementation of helper functions

int run(char const *bam, char const *vcf) {
    // init global variables
    fp_in = hts_open(bam, "r");             //open bam file
    bamHdr = sam_hdr_read(fp_in);           //read header
    bam_file_index = sam_index_load( fp_in, bam );
    aln = bam_init1();                      //initialize an alignment
    
    ifstream vcf_file(vcf);                   // open sv file
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
        int comp, chrom, sv_pos, sv_end;
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
                        splitted_line[0].length() == CIGAR_SOFT_CLIP && (splitted_line[0][3] == 'X' || splitted_line[0][3] == 'x') ) {
                    // chrom is X
                } else if ( ( splitted_line[0].length() && (splitted_line[0][0] == 'Y' || splitted_line[0][0] == 'y') ) ||
                           ( splitted_line[0].length() == CIGAR_SOFT_CLIP && (splitted_line[0][3] == 'Y' || splitted_line[0][3] == 'y') ) ){
                    // chrom is Y
                } else {
                    cout << "# Chromosome could not be converted into integer. CHROM tag : " << splitted_line[0] << endl;
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
        int type_index = splitted_line[7].find("SVTYPE=");
        
        if (type_index == -1) {
            return false;
        }
        
        // Retrieve END
        int sv_end_start = splitted_line[7].find("END=");
        if (sv_end_start != 0) {
            sv_end_start = splitted_line[7].find(";END=");
        }
        int sv_end_end = splitted_line[7].find(";", sv_end_start+1);
        
        if (sv_end_start == -1) {
            return false;
        }
        if (sv_end_end == -1) {
            sv_end_end = splitted_line[7].length();
        }
        
        try {
            if (sv_end_start == 0) {
                sv_end = stoi( splitted_line[7].substr(sv_end_start+4, sv_end_end) );
            } else {
                sv_end = stoi( splitted_line[7].substr(sv_end_start+5, sv_end_end) );
            }
        } catch(...) {
            cout << "# Unable to convert into int" << splitted_line[7].substr(sv_end_start+4, sv_end_end) << endl;
        }
        
        
        // Retrieve CIPOS and CIEND index
        // CIPOS
        int CIPOS_start = splitted_line[7].find("CIPOS=");
        int CIPOS_comma = splitted_line[7].find(",", CIPOS_start + 1);
        int CIPOS_end = splitted_line[7].find(";", CIPOS_start + 1);
        
        if (CIPOS_end == -1)
            CIPOS_end = splitted_line[7].length();
        
        try {
            outer_start = stoi( splitted_line[7].substr(CIPOS_start+6, CIPOS_comma) );
        } catch(...) {
            cout << "# Unable to convert into int" << splitted_line[7].substr(CIPOS_start+6, CIPOS_comma) << endl;
        }
        
        try {
            inner_start = stoi( splitted_line[7].substr(CIPOS_comma+1, CIPOS_end) );
        } catch(...) {
            cout << "# Unable to convert into int" << splitted_line[7].substr(CIPOS_comma+1, CIPOS_end) << endl;
        }

        if ( splitted_line[7].find("INS", type_index + 7) != -1 || splitted_line[7].find("ins", type_index + 7) != -1 ) {
            insertion(id, chrom, sv_pos+outer_start, sv_pos+inner_start);
            return true;
        }
        
        // CIEND
        int CIEND_start = splitted_line[7].find("CIEND=");
        int CIEND_comma = splitted_line[7].find(",", CIEND_start + 1);
        int CIEND_end = splitted_line[7].find(";", CIEND_start + 1);
        
        if (CIEND_end == -1)
            CIEND_end = splitted_line[7].length();
        
        try {
            inner_end = stoi( splitted_line[7].substr(CIEND_start+6, CIEND_comma) );
        } catch(...) {
            cout << "# Unable to convert into int" << splitted_line[7].substr(CIEND_start+6, CIEND_comma) << endl;
        }
        
        try {
            outer_end = stoi( splitted_line[7].substr(CIEND_comma+1, CIEND_end) );
        } catch(...) {
            cout << "# Unable to convert into int" << splitted_line[7].substr(CIEND_comma+1, CIEND_end) << endl;
        }
        
        // Call relevant function to refine SV
        if ( splitted_line[7].find("DEL", type_index + 7) != -1 || splitted_line[7].find("del", type_index + 7) != -1 ) {
            deletion(id, chrom, sv_pos+outer_start, sv_pos+inner_start, sv_end+inner_end, sv_end+outer_end);
        } else if ( splitted_line[7].find("INV", type_index + 7) != -1 || splitted_line[7].find("inv", type_index + 7) != -1 ) {
            inversion(id, chrom, sv_pos+outer_start, sv_pos+inner_start, sv_end+inner_end, sv_end+outer_end);
        } else {
            // Other type of sv
            return false;
        }
    }
    
    return true;
}

int find_start(int chrom, int outer_start, int inner_start) {

    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
    hts_itr_t *iter;
    iter = sam_itr_queryi( bam_file_index, chrom - 1, outer_start - WIDER_INTERVAL, inner_start + NARROW_INTERVAL);
    vector<int> start_positions;
    
    if (iter == NULL) {
        printf("# invalid interval, iter is null\n");
    } else {
        while (sam_itr_next( fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read ends in between (outer_start, inner_start)
            if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == CIGAR_SOFT_CLIP) {
                int reference_pos = pos;
                for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                    if ( bam_cigar_op( cigar[i] ) != CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_PADDING ) {
                        reference_pos += bam_cigar_oplen( cigar[i] );
                    }
                }
                
                if ( outer_start <= reference_pos && reference_pos <= inner_start ) {
                    start_positions.push_back(reference_pos+1);
                }
            }
        }
    }
    sam_itr_destroy(iter);

    return consensus(start_positions);
}

int find_end(int chrom, int inner_end, int outer_end) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;
    
    hts_itr_t *iter;
    vector<int> end_positions;
    iter = sam_itr_queryi( bam_file_index, chrom - 1, inner_end - NARROW_INTERVAL, outer_end + NARROW_INTERVAL);
    
    if (iter == NULL) {
        printf("# invalid interval, iter is null\n");
    } else {
        while (sam_itr_next( fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read starts in between (inner_end, outer_end)
            if ( bam_cigar_op( cigar[0] ) == CIGAR_SOFT_CLIP && inner_end <= pos && pos <= outer_end  ) {
                end_positions.push_back(pos+1);
            }
        }
    }
    
    sam_itr_destroy(iter);
    
    return consensus(end_positions);
}

int find_start_or_end(int chrom, int start, int end) {
    
    // sam_read1 variables
    int32_t pos;
    char *chr;
    size_t read_len;
    uint32_t flag;

    hts_itr_t *iter;
    iter = sam_itr_queryi( bam_file_index, chrom - 1, start - NARROW_INTERVAL, end + NARROW_INTERVAL);
    vector<int> positions;
    
    if (iter == NULL) {
        printf("# invalid interval, iter is null\n");
    } else {
        while (sam_itr_next( fp_in, iter, aln ) > 0) {
            pos = aln->core.pos;
            uint32_t *cigar = bam_get_cigar(aln);
            
            // If read ends in between (start, end)
            if ( bam_cigar_op( cigar[aln->core.n_cigar-1] ) == CIGAR_SOFT_CLIP) {
                int reference_pos = pos;
                for ( int i = 0; i < aln->core.n_cigar; i++ ) {
                    if ( bam_cigar_op( cigar[i] ) != CIGAR_INSERTION && bam_cigar_op( cigar[i] ) != CIGAR_SOFT_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_HARD_CLIP && bam_cigar_op( cigar[i] ) != CIGAR_PADDING ) {
                        reference_pos += bam_cigar_oplen( cigar[i] );
                    }
                }
                if ( start <= reference_pos && reference_pos <= end ) {
                    positions.push_back(reference_pos+1);
                }
            }
            // If read starts in between (start, end)
            if ( bam_cigar_op( cigar[0] ) == CIGAR_SOFT_CLIP && start <= pos && pos <= end ) {
                positions.push_back(pos);
            }
        }
    }
    sam_itr_destroy(iter);
    
    return consensus(positions);
}

int consensus(vector<int> locations) {
    int length = locations.size(), i;
    int consensus = -1, max_count = CONSENSUS_MIN_COUNT - 1;
    vector<int>::iterator ptr_i;
    vector<int>::iterator ptr_j;

    sort(locations.begin(), locations.end());
    
    for (ptr_i = locations.begin(), i = 0; ptr_i<locations.end(); ptr_i++, i++) {
        int item_i = *ptr_i;
        int sum = 0, count = 1;
        for (ptr_j = locations.begin()+i; ptr_j<locations.end() && *ptr_j<*ptr_i+CONSENSUS_INTEVAL; ptr_j++) {
            sum += (*ptr_j) - (*ptr_i);
            count++;
        }
        
        sum = (*ptr_i) + round(sum / count);
        
        if (count > max_count) {
            max_count = count;
            consensus = sum;
        }
    }
    
    return consensus;
}

bool insertion(string id, int chrom, int outer_start, int inner_start) {

    int refined_start = find_start_or_end(chrom, outer_start, inner_start);
    
    cout << chrom << "\t" << id << "\t" << "ins" << "\t" << refined_start << endl;
    
    return true;
}

bool deletion(string id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end) {
    
    int refined_start = find_start(chrom, outer_start, inner_start);
    int refined_end = find_end(chrom, inner_end, outer_end);
    
    cout << chrom << "\t" <<  id << "\t" << "del" << "\t" << refined_start << "\t" <<  refined_end << endl;

    return true;
}

bool inversion(string id, int chrom, int outer_start, int inner_start, int inner_end, int outer_end) {

    // sam_read1 variables
    int refined_start = find_start_or_end(chrom, outer_start, inner_start);
    int refined_end = find_start_or_end(chrom, inner_end, outer_end);
    
    cout << chrom << "\t" << id << "\t" << "inv" << "\t" << refined_start << "\t" <<  refined_end << endl;
    
    return true;
}
