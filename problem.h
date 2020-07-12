#ifndef PROBLEM_H
#define PROBLEM_H
#include "random_generator.h"
#include "graphs.h"

class problem
{
protected:
    random_generator generator;
    bool dominated[MAX_VERTICES];
    refer nondominated_count;
    refer i,j;
public:
    problem();
    refer undominated_vertices_count(graph G, bool *result);
    refer undominated_vertices_count_component(graph G, bool *result, refer current_component, refer *component_labels);
    refer undominated_vertices_count_component_alternative(graph G, refer *candidate_set, refer candidate_set_size, refer component_size);
    refer dominating_set_size(graph G, bool *result);
};

#endif // PROBLEM_H
