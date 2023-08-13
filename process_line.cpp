#ifndef PROCESS_LINE_CPP
#define PROCESS_LINE_CPP

#include "static_variables.h"
#include "params.h"
#include "variations.cpp"


bool process_line(params &_params, std::string &line) {
    
    // extract data
    if ( line.find("CIPOS") != std::string::npos && line.find("SVTYPE=") != std::string::npos ) {
        
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
                _params.verbose && std::cout << " Unable to convert chrom into int for sv id " << _id << std::endl;
                return false;
            }
        }

        // MARK: Processing POS field
        try {
            pos = stoi(_pos);
        } catch (...) {
            _params.verbose && std::cout << " Unable to convert pos into int for sv id " << _id << std::endl;
            return false;
        }

        result res;

        // MARK: Processing INFO field
        int index_1, index_2, end, outer_start, inner_start, inner_end, outer_end;
        std::string temp;

        // SV Type
        index_1 = _info.find("SVTYPE=");
        index_2 = _info.find(";", index_1);
        sv_type = _info.substr( index_1 + 7, index_2 - index_1 - 7 );

        // END
        index_1 = _info.find("END=") != 0 ? _info.find(";END=") + 1 : 0;
        index_2 = _info.find(";", index_1);

        try {
            end = stoi( _info.substr( index_1 + 4, index_2 - index_1 - 4 ) );
        } catch(...) {
            _params.verbose && std::cout << " Unable to convert pos into int for sv id" << _id << std::endl;
            return false;
        }

        // CIPOS
        index_1 = _info.find("CIPOS=");
        index_2 = _info.find(";", index_1);
        iss.clear();
        iss.str( _info.substr( index_1 + 6, index_2 - index_1 - 6 ) );
        
        std::getline(iss, temp, ',');
        
        try {
            outer_start = stoi( temp ) + pos;
        } catch(...) {
            _params.verbose && std::cout << " Unable to convert into int " << temp << " for cipos" << std::endl;
            return false;
        }
        iss >> temp;
        
        try {
            inner_start = stoi( temp ) + pos;
        } catch(...) {
            _params.verbose && std::cout << " Unable to convert into int " << temp << " for cipos" << std::endl;
            return false;
        }


        if ( sv_type == "INS" || sv_type == "INS:ME")
            goto refining_ins;

        if ( sv_type != "DEL" && sv_type != "DEL:ME" && sv_type != "INV" ) {
            _params.verbose && std::cout << " SV Type is " << sv_type << ". Program is unable to handle it." << std::endl;
            return false;
        }

        if ( line.find("CIEND") == std::string::npos ) {
            _params.verbose && std::cout << " Missing CIEND in info field for DEL and INV." << std::endl;
            return false;
        }

        // CIEND
        index_1 = _info.find("CIEND=");    
        index_2 = _info.find(";", index_1);
        iss.clear();
        iss.str( _info.substr( index_1 + 6, index_2 - index_1 - 6 ) );
        
        std::getline(iss, temp, ',');
        try {
            inner_end = stoi( temp ) + end;
        } catch(...) {
            _params.verbose && std::cout << " Unable to convert into int " << temp << " for ciend" << std::endl;
            return false;
        }
        
        iss >> temp;
        try {
            outer_end = stoi( temp ) + end;
        } catch(...) {
            _params.verbose && std::cout << "Unable to convert into int " << temp << " for ciend" << std::endl;
            return false;
        }

        goto refining_other;
        
        refining_ins:
        res = sv::insertion(chrom, outer_start, inner_start, pos, _params);
        if ( res.refined_start != -1 ) {
            index_1 = _info.find("CIPOS=");
            index_2 = _info.find(";", index_1);
            _info.erase( index_1, index_2 - index_1 + 1 );
            _pos = std::to_string(res.refined_start);

            index_1 = _info.find("END=") != 0 ? _info.find(";END=") + 1 : 0;
            index_2 = _info.find(";", index_1);
            _info.erase( index_1, index_2 - index_1 + 1 );
            _info = _info + std::string( ";END=" ) + std::to_string(res.refined_start+1);

            // If there is no consensus on sv length, then it will place SVLEN
            if ( res.svlen != -1 ) {
                if ( _info.find("SVLEN=") != std::string::npos ) {
                    index_1 = _info.find("SVLEN=");
                    index_2 = _info.find(";", index_1);
                    _info.erase( index_1, index_2 - index_1 + 1);
                }
                
                _info = _info + std::string( ";SVLEN=" ) + std::to_string( res.svlen );
            }
        }
        
        _params.out_vcf << _chrom << '\t' << _pos << '\t' <<_id << '\t' << _ref << '\t' << _alt << '\t' << _qual << '\t' << _filter << '\t' << _info << '\t' << _format << std::endl;
        return true;
        
        
        refining_other:
        if ( sv_type == "DEL" || sv_type == "DEL:ME" ) {

            res = sv::deletion(chrom, outer_start, inner_start, inner_end, outer_end, pos, end, _params);

            if ( res.refined_start != -1 ) {
                index_1 = _info.find("CIPOS=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 + 1 );
                _pos = std::to_string(res.refined_start); 
            }
            
            if ( res.refined_end != -1 ) {
                index_1 = _info.find("CIEND=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 + 1 );
            }

            if ( res.refined_start != -1 && res.refined_end != -1 ) {
                if ( _info.find("SVLEN=") != std::string::npos ) {
                    index_1 = _info.find("SVLEN=");
                    index_2 = _info.find(";", index_1);
                    _info.erase( index_1, index_2 - index_1 + 1);
                }
                
                _info = _info + std::string( ";SVLEN=" ) + std::to_string( abs(end - pos + 1) );
            }
            
            _params.out_vcf << _chrom << '\t' << _pos << '\t' <<_id << '\t' << _ref << '\t' << _alt << '\t' << _qual << '\t' << _filter << '\t' << _info << '\t' << _format << std::endl;
            return true;
        } 

        if ( sv_type == "INV" ) {
            
            res = sv::inversion(chrom, outer_start, inner_start, inner_end, outer_end, pos, end, _params);
            
            if ( res.refined_start != -1 ) {
                index_1 = _info.find("CIPOS=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 + 1 );
                _pos = std::to_string(res.refined_start); 
            }
            
            if ( res.refined_end != -1 ) {
                index_1 = _info.find("CIEND=");
                index_2 = _info.find(";", index_1);
                _info.erase( index_1, index_2 - index_1 + 1);
            }

            if ( res.refined_start != -1 && res.refined_end != -1 ) {
                if ( _info.find("SVLEN=") != std::string::npos ) {
                    index_1 = _info.find("SVLEN=");
                    index_2 = _info.find(";", index_1);
                    _info.erase( index_1, index_2 - index_1 + 1);
                }
                
                _info = _info + std::string( ";SVLEN=" ) + std::to_string( abs(end - pos + 1) );
            }
            
            _params.out_vcf << _chrom << '\t' << _pos << '\t' <<_id << '\t' << _ref << '\t' << _alt << '\t' << _qual << '\t' << _filter << '\t' << _info << '\t' << _format << std::endl;
            return true;
        }

        _params.verbose && std::cout << "No SV refiner is called." << std::endl;
        return false;
    }
    
    _params.out_vcf << line;
    return false;
}

#endif