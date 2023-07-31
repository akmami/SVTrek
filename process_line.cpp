#ifndef PROCESS_LINE_CPP
#define PROCESS_LINE_CPP

#include "variables.h"
#include "deletion.cpp"
#include "insertion.cpp"
#include "inversion.cpp"

bool process_line(std::string &line, bool verbose=false) {
    
    // extract data
    if (line.find("IMPRECISE") == std::string::npos || line.find("SVLEN") != std::string::npos || line.find("CIPOS") != std::string::npos) {
        
        // local variables
        int comp, chrom, sv_pos, sv_end;
        std::string id, alt, token;
        
        std::vector<std::string> splitted_line;
        std::istringstream iss(line);

        while(std::getline(iss, token, '\t'))
            splitted_line.push_back(token);
        
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
                    std::cout << "# Chromosome could not be converted into integer. CHROM tag : " << splitted_line[0] << std::endl;
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
            std::cout << "# Unable to convert into int" << splitted_line[7].substr(sv_end_start+4, sv_end_end) << std::endl;
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
            std::cout << "# Unable to convert into int" << splitted_line[7].substr(CIPOS_start+6, CIPOS_comma) << std::endl;
        }
        
        try {
            inner_start = stoi( splitted_line[7].substr(CIPOS_comma+1, CIPOS_end) );
        } catch(...) {
            std::cout << "# Unable to convert into int" << splitted_line[7].substr(CIPOS_comma+1, CIPOS_end) << std::endl;
        }

        if ( splitted_line[7].find("INS", type_index + 7) != -1 || splitted_line[7].find("ins", type_index + 7) != -1 ) {
            insertion(id, chrom, sv_pos+outer_start, sv_pos+inner_start, sv_pos, verbose);
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
            std::cout << "# Unable to convert into int" << splitted_line[7].substr(CIEND_start+6, CIEND_comma) << std::endl;
        }
        
        try {
            outer_end = stoi( splitted_line[7].substr(CIEND_comma+1, CIEND_end) );
        } catch(...) {
            std::cout << "# Unable to convert into int" << splitted_line[7].substr(CIEND_comma+1, CIEND_end) << std::endl;
        }
        
        // Call relevant function to refine SV
        if ( splitted_line[7].find("DEL", type_index + 7) != -1 || splitted_line[7].find("del", type_index + 7) != -1 ) {
            deletion(id, chrom, sv_pos+outer_start, sv_pos+inner_start, sv_end+inner_end, sv_end+outer_end, sv_pos, sv_end, verbose);
        } else if ( splitted_line[7].find("INV", type_index + 7) != -1 || splitted_line[7].find("inv", type_index + 7) != -1 ) {
            inversion(id, chrom, sv_pos+outer_start, sv_pos+inner_start, sv_end+inner_end, sv_end+outer_end, sv_pos, sv_end, verbose);
        } else {
            // Other type of sv
            return false;
        }
    }
    
    return true;
}

#endif