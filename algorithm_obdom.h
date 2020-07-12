#ifndef ALGORITHM_OBDOM_H
#define ALGORITHM_OBDOM_H

#include "algorithm.h"
#include "problem.h"
#include "algorithm_greedydom.h"

class algorithm_obdom : public algorithm
{
private:
    refer current_permutation[MAX_VERTICES];
    refer new_permutation[MAX_VERTICES];
    refer possibilities[MAX_VERTICES];
    refer possibilities_count;
    bool dominated[MAX_VERTICES];
    problem *problemdom;
    algorithm_greedydom *greedydom;
public:
    algorithm_obdom();
    ~algorithm_obdom();
    void greedy_obdom(graph G, refer *permutation, bool *result, bool connected=false);
    void obdom(graph G, bool *result, long long time_limit, refer lower_bound, long long *output_t, bool multi_start=false, bool init_greedy=true, bool use_t_stag_max=false, bool connected=false);
};

#endif // ALGORITHM_OBDOM_H
