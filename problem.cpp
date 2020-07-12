#include "problem.h"

problem::problem()
{
    refer i;
    for (i=0;i<MAX_VERTICES;i++)
    {
        dominated[i] = 0;
    }
}

refer problem::undominated_vertices_count(graph G, bool *result)
{
    refer count, v, j;
    bool dominated;

    count = 0;
    for (v=0; v<G->n; v++)
    {
        if (result[v] == 0)
        {
            dominated = false;
            for (j=0; j<G->V[v].edgecount; j++)
            {
                if (result[G->V[v].sibl[j]] == 1)
                {
                    dominated = true;
                    break;
                }
            }
            if (! dominated)
            {
                count++;
            }
        }
    }

    return count;
}

refer problem::undominated_vertices_count_component(graph G, bool *result, refer current_component, refer *component_labels)
{
    refer count, v, j;
    bool dominated;

    count = 0;
    for (v=0; v<G->n; v++)
    {
        if (component_labels[v] != current_component)
        {
            continue;
        }
        if (result[v] == 0)
        {
            dominated = false;
            for (j=0; j<G->V[v].edgecount; j++)
            {
                if (component_labels[G->V[v].sibl[j]] != current_component)
                {
                    continue;
                }
                if (result[G->V[v].sibl[j]] == 1)
                {
                    dominated = true;
                    break;
                }
            }
            if (! dominated)
            {
                count++;
            }
        }
    }

    return count;
}

refer problem::undominated_vertices_count_component_alternative(graph G, refer *candidate_set, refer candidate_set_size, refer component_size)
{
    nondominated_count = component_size;

    for (i=0;i<candidate_set_size;i++)
    {
        if (! dominated[candidate_set[i]])
        {
            nondominated_count--;
            if (0 == nondominated_count)
            {
                break;
            }
            dominated[candidate_set[i]] = 1;
        }
        for (j=0;j<G->V[candidate_set[i]].edgecount;j++)
        {
            if (! dominated[G->V[candidate_set[i]].sibl[j]])
            {
                nondominated_count--;
                if (0 == nondominated_count)
                {
                    break;
                }
                dominated[G->V[candidate_set[i]].sibl[j]] = 1;
            }
        }
    }
    for (i=0;i<candidate_set_size;i++)
    {
        dominated[candidate_set[i]] = 0;
        for (j=0;j<G->V[candidate_set[i]].edgecount;j++)
        {
            dominated[G->V[candidate_set[i]].sibl[j]] = 0;
        }
    }

    return nondominated_count;
}

refer problem::dominating_set_size(graph G, bool *result)
{
    refer count, v;

    count = 0;

    // MWDS
    if (G->weighted)
    {
        for (v=0; v<G->n; v++)
        {
            if (result[v] == 1)
            {
                count += G->V[v].weight;
            }
        }
    }
    // MDS
    else
    {
        for (v=0; v<G->n; v++)
        {
            if (result[v] == 1)
            {
                count++;
            }
        }
    }

    return count;
}
