#include "algorithm_acols.h"
#include "algorithm_greedydom.h"
#include <QTime>
#define MAX_ANTS 30

double pheromone[MAX_VERTICES] = {0.0};
bool ants[MAX_ANTS][MAX_VERTICES];
long roulette_aco[MAX_VERTICES] = {0};
refer redundant[MAX_VERTICES] = {0};
bool dominated_vertices[MAX_VERTICES] = {0};
refer dominated_ways[MAX_VERTICES] = {0};
bool initial_greedy[MAX_VERTICES] = {0};

algorithm_acols::algorithm_acols()
{
}

algorithm_acols::~algorithm_acols()
{
}

void algorithm_acols::acols(graph G, bool *result, long long time_limit, unsigned long max_iter, refer max_ants, bool greedy_init, bool preprocessing, refer lower_bound, bool sparse, double update_param1, double update_param2, long long *output_t)
{
    refer i,j,k,m,a,r,v;
    refer F = G->n;
    double pheromone_initial = 10.0;
    double rho = 0.985;
    double pheromone_sum, pheromone_partial;
    long p_r_percent = 60;
    int M = 100;
    bool initial = true;
    refer t_stag, t_stag_max;
    bool ok;

    if (sparse)
    {
        //t_stag_max = 50000;
        t_stag_max = 1000000;
    }
    else
    {
        //t_stag_max = 1000;
        t_stag_max = 20000;
    }

    random_generator generator;
    long q, roulette_current;
    refer possibilities_count, vertex, vertex_index;
    refer best_ant_index, best_ant_cardinality, best_ant_count, cardinality;
    bool is_redundant;
    refer dominated_count;

    QTime qTimer;
    qTimer.start();

    refer redundant_count, lowest_degree, lowest_degree_count;

    for (j=0;j<G->n;j++)
    {
        result[j] = 1;
    }

    if (greedy_init)
    {
        for (j=0;j<G->n;j++)
        {
            pheromone[j] = pheromone_initial;
        }
        algorithm_greedydom *greedydom = new algorithm_greedydom();
        greedydom->greedydom(G, initial_greedy);
        for (j=0;j<G->n;j++)
        {
            if (initial_greedy[j])
            {
                pheromone[j] = pheromone_initial*100.0;
            }
            else
            {
                pheromone[j] = pheromone_initial;
            }
        }
        delete(greedydom);
    }
    else
    {
        for (j=0;j<G->n;j++)
        {
            pheromone[j] = pheromone_initial;
        }
    }

    // preprocessing using maximal independent sets
    if (preprocessing)
    {
        bool *used_in_independent_set = new bool[G->n];
        bool *already_covered = new bool[G->n];
        refer *independent_set_candidates = new refer[G->n];
        refer *independent_set_candidates_positioning = new refer[G->n];
        refer independent_set_candidates_count;

        for (i=0;i<M;i++)
        {
            // construct a maximal independent set using a greedy algorithm
            for (j=0;j<G->n;j++)
            {
                used_in_independent_set[j] = false;
                independent_set_candidates[j] = j;
                independent_set_candidates_positioning[j] = j;
            }
            independent_set_candidates_count = G->n;

            // choose the vertex
            while (independent_set_candidates_count > 0)
            {
                r = generator.random(0, independent_set_candidates_count-1);
                v = independent_set_candidates[r];
                used_in_independent_set[v] = true;
                already_covered[v] = true;

                // exclude the vertex and neighbours from the possibilities
                independent_set_candidates[independent_set_candidates_positioning[v]] = independent_set_candidates[independent_set_candidates_count-1];
                independent_set_candidates_positioning[independent_set_candidates[independent_set_candidates_count-1]] = independent_set_candidates_positioning[v];
                independent_set_candidates_count--;
                for (k=0;k<G->V[v].edgecount;k++)
                {
                    if (! already_covered[G->V[v].sibl[k]])
                    {
                        independent_set_candidates[independent_set_candidates_positioning[G->V[v].sibl[k]]] = independent_set_candidates[independent_set_candidates_count-1];
                        independent_set_candidates_positioning[independent_set_candidates[independent_set_candidates_count-1]] = independent_set_candidates_positioning[G->V[v].sibl[k]];
                        independent_set_candidates_count--;
                        already_covered[G->V[v].sibl[k]] = true;
                    }
                }
            }

            // update the pheromone using the maximal independent set
            for (j=0;j<G->n;j++)
            {
                if (used_in_independent_set[j])
                {
                    pheromone[j] += update_param1 / (update_param2 + (double) best_ant_cardinality - (double) F);
                }
                if (pheromone[j] < 0.08)
                {
                    pheromone[j] = 0.08;
                }
            }
        }

        delete[](used_in_independent_set);
        delete[](independent_set_candidates);
        delete[](independent_set_candidates_positioning);
        delete[](already_covered);
    }

    t_stag = 0;
    for (i=0;i<max_iter;i++)
    {
        if (qTimer.elapsed() > time_limit /*|| t_stag >= t_stag_max*/ || F == lower_bound)
        {
            break;
        }
        // empty ants in the beginning
        for (a=0;a<max_ants;a++)
        {
            for (j=0;j<G->n;j++)
            {
                ants[a][j] = 0;
            }
        }
        for (a=0;a<max_ants;a++)
        {
            dominated_count = 0;
            for (j=0;j<G->n;j++)
            {
                dominated_vertices[j] = 0;
            }
            initial = true;
            while (1)
            {
                // test if the current ant represents a dominating set
                if (dominated_count >= G->n)
                {
                    break;
                }

                // all possibilities
                if (! initial)
                {
                    possibilities_count = 0;
                    pheromone_sum = 0.0;
                    for (k=0;k<G->V[vertex].edgecount;k++)
                    {
                        j = G->V[vertex].sibl[k];
                        if (0 == ants[a][j])
                        {
                            pheromone_sum += pheromone[j];
                            possibilities_count++;
                        }
                    }
                    if (0 != possibilities_count)
                    {
                        // computing the borders for the roulette
                        roulette_current = 0;
                        pheromone_partial = 0.0;
                        q = 0;
                        for (k=0;k<G->V[vertex].edgecount;k++)
                        {
                            j = G->V[vertex].sibl[k];
                            if (0 == ants[a][j])
                            {
                                roulette_current = (long) (pheromone_partial * (double) (RAND_MAX) / pheromone_sum);
                                roulette_aco[q] = roulette_current;
                                pheromone_partial += pheromone[j];
                                q++;
                            }
                        }
                    }
                    else
                    {
                        initial = true;
                    }
                }
                if (initial)
                {
                    possibilities_count = 0;
                    pheromone_sum = 0.0;
                    for (j=0;j<G->n;j++)
                    {
                        if (0 == ants[a][j])
                        {
                            pheromone_sum += pheromone[j];
                            possibilities_count++;
                        }
                    }

                    // computing the borders for the roulette
                    roulette_current = 0;
                    pheromone_partial = 0.0;
                    q = 0;
                    for (j=0;j<G->n;j++)
                    {
                        if (0 == ants[a][j])
                        {
                            roulette_current = (long) (pheromone_partial * (double) (RAND_MAX) / pheromone_sum);
                            roulette_aco[q] = roulette_current;
                            pheromone_partial += pheromone[j];
                            q++;
                        }
                    }

                    // in the case of sparse version of the algorithm,
                    // we now switch to the condensed roulette
                    if (sparse)
                    {
                        initial = false;
                    }
                }

                // choosing the index of 0-bits to choose
                q = generator.random(0,RAND_MAX);
                vertex_index = 0;
                vertex = 0;
                for (j=0;j<possibilities_count;j++)
                {
                    if (roulette_aco[j] <= q && (j == possibilities_count-1 || q < roulette_aco[j+1]))
                    {
                        vertex_index = j;
                        break;
                    }
                }

                // finding to which vertex this corresponds
                for (j=0;j<G->n;j++)
                {
                    if (0 == ants[a][j])
                    {
                        if (0 == vertex_index)
                        {
                            vertex = j;
                            break;
                        }
                        else
                        {
                            vertex_index--;
                        }
                    }
                }

                // putting the vertex to the set
                if (0 == ants[a][vertex])
                {
                    // TODO EXPERIMENT: but only if it is currently not dominated
                    /*ok = false;
                    if (0 == dominated_vertices[vertex])
                    {
                        ok = true;
                    }
                    for (k=0;k<G->V[vertex].edgecount;k++)
                    {
                        if (0 == dominated_vertices[G->V[vertex].sibl[k]])
                        {
                            ok = true;
                        }
                    }
                    if (! ok)
                    {
                        continue;
                    }*/
                    ants[a][vertex] = 1;
                    if (0 == dominated_vertices[vertex])
                    {
                        dominated_vertices[vertex] = 1;
                        dominated_count++;
                    }
                    for (k=0;k<G->V[vertex].edgecount;k++)
                    {
                        if (0 == dominated_vertices[G->V[vertex].sibl[k]])
                        {
                            dominated_count++;
                            dominated_vertices[G->V[vertex].sibl[k]] = 1;
                        }
                    }
                }
            }

            // LS - iterative improvement of the ant
            redundant_count = 0;
            for (j=0;j<G->n;j++)
            {
                dominated_ways[j] = 0;
            }
            for (j=0;j<G->n;j++)
            {
                if (1 == ants[a][j])
                {
                    dominated_ways[j]++;
                }
                for (k=0;k<G->V[j].edgecount;k++)
                {
                    if (1 == ants[a][G->V[j].sibl[k]])
                    {
                        dominated_ways[j]++;
                    }
                }
            }
            for (j=0;j<G->n;j++)
            {
                if (1 == ants[a][j])
                {
                    if (1 >= dominated_ways[j])
                    {
                        is_redundant = false;
                    }
                    else
                    {
                        is_redundant = true;
                        for (k=0;k<G->V[j].edgecount;k++)
                        {
                            if (1 >= dominated_ways[G->V[j].sibl[k]])
                            {
                                is_redundant = false;
                                break;
                            }
                        }
                    }
                    if (is_redundant)
                    {
                        redundant[redundant_count] = j;
                        redundant_count++;
                    }
                }
            }

            while (1)
            {
                // getting the information on the redundant nodes
                if (0 == redundant_count)
                {
                    break;
                }

                if (generator.random(0,99) < p_r_percent)
                {
                    // random choice
                    q = generator.random(0,redundant_count-1);
                }
                else
                {
                    // lowest degree choice
                    q = 0;
                    lowest_degree = G->n+1;
                    lowest_degree_count = 1;
                    for (j=0;j<redundant_count;j++)
                    {
                        if (lowest_degree > G->V[redundant[j]].edgecount)
                        {
                            lowest_degree = G->V[redundant[j]].edgecount;
                            lowest_degree_count = 1;
                        }
                        else if (lowest_degree == G->V[redundant[j]].edgecount)
                        {
                            lowest_degree_count++;
                        }
                    }
                    // randomized choice among the lower degree nodes
                    k = generator.random(0,lowest_degree_count-1);
                    for (j=0;j<redundant_count;j++)
                    {
                        if (lowest_degree == G->V[redundant[j]].edgecount)
                        {
                            if (0 == k)
                            {
                                q = j;
                                break;
                            }
                            else
                            {
                                k--;
                            }
                        }
                    }
                }

                vertex = redundant[q];
                ants[a][vertex] = 0;
                redundant[q] = redundant[redundant_count-1];
                redundant_count--;

                // updating the ways of domination
                dominated_ways[vertex]--;
                for (k=0;k<G->V[vertex].edgecount;k++)
                {
                    dominated_ways[G->V[vertex].sibl[k]]--;
                }

                // some other nodes may have also been rendered non-redundant
                m = 0;
                while (m < redundant_count)
                {
                    j = redundant[m];
                    if (1 >= dominated_ways[j])
                    {
                        is_redundant = false;
                    }
                    else
                    {
                        is_redundant = true;
                        for (k=0;k<G->V[j].edgecount;k++)
                        {
                            if (1 >= dominated_ways[G->V[j].sibl[k]])
                            {
                                is_redundant = false;
                                break;
                            }
                        }
                    }
                    if (! is_redundant)
                    {
                        redundant[m] = redundant[redundant_count-1];
                        redundant_count--;
                    }
                    else
                    {
                        m++;
                    }
                }
            }
        }

        // best ant
        best_ant_cardinality = G->n+1;
        best_ant_index = 0;
        best_ant_count = 0;
        for (a=0;a<max_ants;a++)
        {
            cardinality = 0;
            for (j=0;j<G->n;j++)
            {
                cardinality += ants[a][j];
            }

            if (best_ant_cardinality > cardinality)
            {
                best_ant_cardinality = cardinality;
                best_ant_count = 1;
            }
            else if (best_ant_cardinality == cardinality)
            {
                best_ant_count++;
            }
        }

        q = generator.random(0,best_ant_count-1);
        for (a=0;a<max_ants;a++)
        {
            cardinality = 0;
            for (j=0;j<G->n;j++)
            {
                cardinality += ants[a][j];
            }

            if (best_ant_cardinality == cardinality)
            {
                if (0 == q)
                {
                    best_ant_index = a;
                    break;
                }
                else
                {
                    q--;
                }
            }
        }

        // randomized choice of the best ant
        for (a=0;a<max_ants;a++)
        {
            cardinality = 0;
            for (j=0;j<G->n;j++)
            {
                cardinality += ants[a][j];
            }

            if (best_ant_cardinality == cardinality)
            {
                best_ant_cardinality = cardinality;
            }
        }

        if (F > best_ant_cardinality)
        {
            F = best_ant_cardinality;
            for (j=0;j<G->n;j++)
            {
                result[j] = ants[best_ant_index][j];
            }
            printf("%u: [0,%u]\n", i, F);
            t_stag = 0;
        }

        // pheromone update
        for (j=0;j<G->n;j++)
        {
            // evaporation
            pheromone[j] *= rho;
        }
        // reinforcement of nodes in the best ant
        for (j=0;j<G->n;j++)
        {
            if (1 == ants[best_ant_index][j])
            {
                pheromone[j] += update_param1 / (update_param2 + (double) best_ant_cardinality - (double) F);
            }
            if (pheromone[j] < 0.08)
            {
                pheromone[j] = 0.08;
            }
        }

        t_stag++;
    }

    *output_t = (i * max_ants);
}
