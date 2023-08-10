#ifndef PROCESS_LINE_CPP
#define PROCESS_LINE_CPP

#include "static_variables.h"
#include "params.h"
#include "variations.cpp"


bool process_line(params &_params, std::string &line) {
    
    // extract data
    if (line.find("IMPRECISE") == std::string::npos || line.find("SVLEN") != std::string::npos || line.find("CIPOS") != std::string::npos) {
        
        // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
        std::string _chrom, _pos, _id, _ref, _alt, _qual, _filter, _info, _format;
        std::istringstream iss(line);

        std::getline(iss, _chrom, '\t');
        std::getline(iss, _pos, '\t');
        std::getline(iss, _id, '\t');
        std::getline(iss, _ref, '\t');
        std::getline(iss, _alt, '\t');
        std::getline(iss, _qual, '\t');
        std::getline(iss, _filter, '\t');
        std::getline(iss, _info, '\t');
        iss >> _format;
        
        int chrom, pos;
        std::string sv_type;

        // MARK: Processing CHROM field
        try {
            chrom = stoi(_chrom);
        } catch (...) {
            try {
                chrom = stoi(_chrom.substr(3));
            } catch (...) {
                return false;
            }
        }

        // MARK: Processing POS field
        try {
            pos = stoi(_pos);
        } catch (...) {
            return false;
        }

        result res;

        // MARK: Processing INFO field
        int index_1, index_2, end, outer_start, inner_start, inner_end, outer_end;
        std::string temp;

        // END
        index_1 = _info.find("SVTYPE=");
        if (index_1 == std::string::npos) return false;
        index_2 = _info.find(";", index_1);
        sv_type = _info.substr( index_1 + 7, index_2 );

        index_1 = _info.find("END=") != 0 ? _info.find(";END=") : 0;
        if (index_1 == std::string::npos) return false;
        index_2 = _info.find(";", index_1);

        try {
            end = stoi( _info.substr( index_1 + (index_1 == 0 ? 4 : 5), index_2 ) );
        } catch(...) {
            _params.verbose && std::cout << "# Unable to convert into int" << _info.substr( index_1 + (index_1 == 0 ? 4 : 5), index_2 ) << std::endl;
        }

        // CIPOS
        index_1 = _info.find("CIPOS=");
        if (index_1 == std::string::npos) return false;
        index_2 = _info.find(";", index_1);
        iss.clear();
        iss.str( _info.substr( index_1 + 6, index_2 ) );
        
        std::getline(iss, temp, ',');
        
        try {
            outer_start = stoi( temp ) + pos;
        } catch(...) {
            _params.verbose && std::cout << "# Unable to convert into int " << temp << std::endl;
        }
        iss >> temp;
        
        try {
            inner_start = stoi( temp ) + pos;
        } catch(...) {
            _params.verbose && std::cout << "# Unable to convert into int " << temp << std::endl;
        }


        if ( sv_type == "INS" || sv_type == "INS:ME")
            goto refining_ins;

        // CIEND
        index_1 = _info.find("CIEND=");    
        if (index_1 == std::string::npos) return false;
        index_2 = _info.find(";", index_1);
        iss.clear();
        iss.str( _info.substr( index_1 + 6, index_2 ) );
        
        std::getline(iss, temp, ',');
        try {
            inner_end = stoi( temp ) + end;
        } catch(...) {
            _params.verbose && std::cout << "# Unable to convert into int " << temp << std::endl;
        }
        
        iss >> temp;
        try {
            outer_end = stoi( temp ) + end;
        } catch(...) {
            _params.verbose && std::cout << "# Unable to convert into int " << temp << std::endl;
        }

        goto refining_other;
        
        refining_ins:
        res = sv::insertion(chrom, outer_start, inner_start, pos, _params);
        return true;
        
        
        refining_other:
        if ( sv_type == "DEL" || sv_type == "DEL:ME" ) {
            
            res = sv::deletion(chrom, outer_start, inner_start, inner_end, outer_end, pos, end, _params);

            if ( res.refined_start != -1 ) {
                index_1 = _info.find("CIPOS=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 );
                _pos = std::to_string(res.refined_start); 
            }
            
            if ( res.refined_end != -1 ) {
                index_1 = _info.find("CIEND=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 );
            }
            
            _params.out_vcf << _chrom << '\t' << _pos << _id << _ref << _alt << _qual << _filter << _info << _format << std::endl;
            
            return true;
        } 

        if ( sv_type == "INV" ) {
            
            res = sv::inversion(chrom, outer_start, inner_start, inner_end, outer_end, pos, end, _params);
            
            if ( res.refined_start != -1 ) {
                index_1 = _info.find("CIPOS=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 );
                _pos = std::to_string(res.refined_start); 
            }
            
            if ( res.refined_end != -1 ) {
                index_1 = _info.find("CIEND=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 );
            }
            
            _params.out_vcf << _chrom << '\t' << _pos << _id << _ref << _alt << _qual << _filter << _info << _format << std::endl;

            return true;
        } 

        return false;
    }
    
    _params.out_vcf << line;
    return true;
}

#endif