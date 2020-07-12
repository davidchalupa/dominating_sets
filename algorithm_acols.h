#ifndef ALGORITHM_ACOLS_H
#define ALGORITHM_ACOLS_H

#include "algorithm.h"

class algorithm_acols : public algorithm
{
public:
    algorithm_acols();
    ~algorithm_acols();
    void acols(graph G, bool *result, long long time_limit, unsigned long max_iter, refer max_ants, bool greedy_init, bool preprocessing, refer lower_bound, bool sparse, double update_param1, double update_param2, long long *output_t);
};

#endif // ALGORITHM_ACOLS_H
