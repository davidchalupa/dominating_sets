#include "algorithm_greedydom.h"

algorithm_greedydom::algorithm_greedydom()
{
}

algorithm_greedydom::~algorithm_greedydom()
{
}

void algorithm_greedydom::greedydom(graph G, bool *result)
{
    refer j, k, v, chosen_vertex;
    double max_priority;
    bool ok;

    for (v=0; v<G->n; v++)
    {
        if (G->weighted)
        {
            priorities[v] = (double) (G->V[v].edgecount) / (double) (G->V[v].weight);
        }
        else
        {
            priorities[v] = (double) (G->V[v].edgecount + 1);
        }
        dominated[v] = 0;
        result[v] = 0;
    }

    while (1)
    {
        ok = true;
        for (v=0; v<G->n; v++)
        {
            if (! dominated[v])
            {
                ok = false;
                break;
            }
        }
        if (ok)
        {
            //printf("count %ld\n", dominated_count);
            break;
        }
        // maximum priority
        max_priority = 0;
        for (v=0; v<G->n; v++)
        {
            if (result[v] == 0 && max_priority < priorities[v])
            {
                max_priority = priorities[v];
            }
        }
        // all maximum priority vertices
        max_priority_vertices_count = 0;
        for (v=0; v<G->n; v++)
        {
            if (result[v] == 0 && max_priority == priorities[v])
            {
                max_priority_vertices[max_priority_vertices_count] = v;
                max_priority_vertices_count++;
            }
        }
        // random choice from the maximum priority vertices
        chosen_vertex = max_priority_vertices[generator.random(0, max_priority_vertices_count-1)];
        //chosen_vertex = max_priority_vertices[0];

        // priorities update
        if (G->weighted)
        {
            for (j=0; j<G->V[chosen_vertex].edgecount; j++)
            {
                v = G->V[chosen_vertex].sibl[j];
                // if the chosen vertex was not dominated, then its neighbors now have lower priorities
                if (! dominated[chosen_vertex])
                {
                    priorities[v] -= 1.0 / (double) G->V[v].weight;
                }
                for (k=0; k<G->V[v].edgecount; k++)
                {
                    // if the neighbor was not dominated but now is, then each other neighbor has a lower priority now
                    if (! dominated[v])
                    {
                        priorities[G->V[v].sibl[k]] -= 1.0 / (double) G->V[G->V[v].sibl[k]].weight;
                    }
                }
            }
        }
        else
        {
            if (! dominated[chosen_vertex])
            {
                priorities[chosen_vertex]--;
            }
            for (j=0; j<G->V[chosen_vertex].edgecount; j++)
            {
                v = G->V[chosen_vertex].sibl[j];
                // if the chosen vertex was not dominated, then its neighbors now have lower priorities
                if (! dominated[chosen_vertex])
                {
                    priorities[v]--;
                }
                if (! dominated[v])
                {
                    priorities[v]--;
                }
                for (k=0; k<G->V[v].edgecount; k++)
                {
                    // if the neighbor was not dominated but now is, then each other neighbor has a lower priority now
                    if (! dominated[v])
                    {
                        priorities[G->V[v].sibl[k]]--;
                    }
                }
            }
        }

        // now dominated vertices
        result[chosen_vertex] = 1;
        dominated[chosen_vertex] = 1;
        for (j=0; j<G->V[chosen_vertex].edgecount; j++)
        {
            v = G->V[chosen_vertex].sibl[j];
            dominated[v] = 1;
        }

    }
}
