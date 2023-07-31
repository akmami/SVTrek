#ifndef CONSENSUS_CPP
#define CONSENSUS_CPP

#include <vector>
#include "variables.h"

int consensus(std::vector<int> locations, int imprecise_pos) {

    int length = locations.size(), i;
    int consensus = -1, max_count = CONSENSUS_MIN_COUNT - 1, distance = INT_MAX;
    std::vector<int>::iterator ptr_i;
    std::vector<int>::iterator ptr_j;

    sort(locations.begin(), locations.end());

    for (ptr_i = locations.begin(), i = 0; ptr_i<locations.end(); ptr_i++, i++) {
        int item_i = *ptr_i;
        int sum = 0, count = 1;
        for (ptr_j = locations.begin()+i; ptr_j<locations.end() && *ptr_j<=*ptr_i+CONSENSUS_INTEVAL; ptr_j++) {
            sum += (*ptr_j) - (*ptr_i);
            count++;
        }
        
        sum = (*ptr_i) + round(sum / count);
        
        if (count > max_count) {
            max_count = count;
            consensus = sum;
            distance = std::abs(imprecise_pos-sum);
        } else if (count == max_count && distance > std::abs(imprecise_pos-sum)) {
            max_count = count;
            consensus = sum;
            distance = std::abs(imprecise_pos-sum);
        }
    }
    
    return consensus;
}

#endif