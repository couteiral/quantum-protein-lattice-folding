#include <string>

#include "solution.h"

bool cmp_solution(Solution s1, Solution s2) {

    if (s1.energy > s2.energy) {
        return true;
    } else {
        return false;
    }

}

