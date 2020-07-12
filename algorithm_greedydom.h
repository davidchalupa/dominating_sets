#ifndef ALGORITHM_GREEDYDOM_H
#define ALGORITHM_GREEDYDOM_H

#include "algorithm.h"

class algorithm_greedydom : public algorithm
{
private:
    double priorities[MAX_VERTICES];
    refer max_priority_vertices[MAX_VERTICES];
    refer max_priority_vertices_count;
    bool dominated[MAX_VERTICES];
public:
    algorithm_greedydom();
    ~algorithm_greedydom();
    void greedydom(graph G, bool *result);
};

#endif // ALGORITHM_GREEDYDOM_H
