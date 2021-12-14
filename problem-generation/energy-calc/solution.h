#ifndef SOLUTION_H
#define SOLUTION_H

struct Solution {
    double      energy;
    std::string seq;
};

bool cmp_solution(Solution s1, Solution s2);

#endif
