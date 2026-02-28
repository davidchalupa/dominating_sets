#include "algorithm_obdom.h"

#include <QTime>

algorithm_obdom::algorithm_obdom()
{
    problemdom = new problem();
    greedydom = new algorithm_greedydom();
}

algorithm_obdom::~algorithm_obdom()
{
    delete(problemdom);
    delete(greedydom);
}

void algorithm_obdom::greedy_obdom(graph G, refer *permutation, bool *result, bool connected)
{
    refer i,j,v,dominated_count;
    bool ok, ok_connected;

    for (v=0;v<G->n;v++)
    {
        dominated[v] = 0;
        result[v] = 0;
    }
    dominated_count = 0;

    if (! connected)
    {
        for (i=0;i<G->n;i++)
        {
            if (dominated_count >= G->n)
            {
                break;
            }

            v = permutation[i];
            ok = false;
            if (! dominated[v])
            {
                ok = true;
            }
            for (j=0;j<G->V[v].edgecount;j++)
            {
                if (! dominated[G->V[v].sibl[j]])
                {
                    ok = true;
                    break;
                }
            }
            if (ok)
            {
                result[v] = 1;
                if (! dominated[v])
                {
                    dominated[v] = 1;
                    dominated_count++;
                }
                for (j=0;j<G->V[v].edgecount;j++)
                {
                    if (! dominated[G->V[v].sibl[j]])
                    {
                        dominated[G->V[v].sibl[j]] = 1;
                        dominated_count++;
                    }
                }
            }
        }
    }
    else
    {
        while (dominated_count < G->n)
        {
            for (i=0;i<G->n;i++)
            {
                if (dominated_count >= G->n)
                {
                    break;
                }

                v = permutation[i];
                ok = false;
                if (! dominated[v])
                {
                    ok = true;
                }
                for (j=0;j<G->V[v].edgecount;j++)
                {
                    if (! dominated[G->V[v].sibl[j]])
                    {
                        ok = true;
                        break;
                    }
                }
                ok_connected = false;
                // connectedness of the dominating set
                if (dominated_count == 0)
                {
                    ok_connected = true;
                }
                else
                {
                    // at least one neighbor must be in the dominating set for it to remain connected
                    for (j=0;j<G->V[v].edgecount;j++)
                    {
                        if (result[G->V[v].sibl[j]])
                        {
                            ok_connected = true;
                            break;
                        }
                    }
                }
                if (ok && ok_connected)
                {
                    result[v] = 1;
                    if (! dominated[v])
                    {
                        dominated[v] = 1;
                        dominated_count++;
                    }
                    for (j=0;j<G->V[v].edgecount;j++)
                    {
                        if (! dominated[G->V[v].sibl[j]])
                        {
                            dominated[G->V[v].sibl[j]] = 1;
                            dominated_count++;
                        }
                    }
                }
            }
        }
    }
}

void algorithm_obdom::obdom(graph G, bool *result, long long time_limit, refer lower_bound, long long *output_t, bool multi_start, bool init_greedy, bool use_t_stag_max, bool connected)
{
    refer i,q,v,f1_o,f2_o,f1_n,f2_n,f1_b,f2_b,f1_bc,f2_bc,initial_vertices_count,jumped_value;
    unsigned long long t, t_stag;
    unsigned long long t_stag_max = 1000000;
    unsigned long long t_stag_max_multi_start = 10*G->n;
    unsigned long long t_stag_max_multi_start_best = 5000*G->n;
    unsigned long long cycles_max = 20000;
    unsigned long long cycles;
    refer init_selector;
    bool use_current_greedy_init;
    bool use_current_best_init;
    bool extended_waiting_time;

    QTime qTimer;
    qTimer.start();

    refer *best_permutation = new refer[G->n];

    // initial greedy solution
    greedydom->greedydom(G, result);
    possibilities_count = 0;
    initial_vertices_count = 0;
    if (multi_start)
    {
        use_current_greedy_init = (generator.random(0,100) < 0);
    }
    else
    {
        // non-multi-start version (classic RLS_o) should use the greedy permutation at the start
        use_current_greedy_init = true;
    }

    for (v=0;v<G->n;v++)
    {
        if (init_greedy && use_current_greedy_init && result[v])
        {
            current_permutation[initial_vertices_count] = v;
            initial_vertices_count++;
        }
        else
        {
            possibilities[possibilities_count] = v;
            possibilities_count++;
        }
    }

    for (i=initial_vertices_count;i<G->n;i++)
    {
        v = generator.random(0, possibilities_count-1);
        current_permutation[i] = possibilities[v];
        possibilities[v] = possibilities[possibilities_count-1];
        possibilities_count--;
    }

    greedy_obdom(G, current_permutation, result, connected);
    f1_o = problemdom->undominated_vertices_count(G, result);
    f2_o = problemdom->dominating_set_size(G, result);

    for (v=0;v<G->n;v++)
    {
        best_permutation[v] = current_permutation[v];
    }
    f1_b = f1_o;
    f2_b = f2_o;
    f1_bc = f1_o;
    f2_bc = f2_o;

    printf("[%u,%u]\n", f1_b, f2_b);

    t = 0;
    cycles = 0;
    t_stag = 0;
    extended_waiting_time = false;
    while (qTimer.elapsed() < time_limit && cycles < cycles_max && (! use_t_stag_max || t_stag < t_stag_max))
    {
        if (f1_o == 0 && f2_o == lower_bound)
        {
            break;
        }

        // restart if multi-start version is used
        if (multi_start && (extended_waiting_time && t_stag >= t_stag_max_multi_start_best || ! extended_waiting_time && t_stag >= t_stag_max_multi_start))
        {
            //printf("restart\n");
            // initial greedy solution
            greedydom->greedydom(G, result);
            possibilities_count = 0;
            initial_vertices_count = 0;
            init_selector = generator.random(0,99);
            use_current_best_init = (init_selector < 0);
            // 0% chance of restarting from a promising solution
            if (use_current_best_init)
            {
                for (v=0;v<G->n;v++)
                {
                    current_permutation[v] = best_permutation[v];
                }
            }
            else
            {
                // 25% chance of starting greedily
                // 75% chance of starting entirely randomly
                use_current_greedy_init = (init_selector < 0);
                for (v=0;v<G->n;v++)
                {
                    if (init_greedy && use_current_greedy_init && result[v])
                    {
                        current_permutation[initial_vertices_count] = v;
                        initial_vertices_count++;
                    }
                    else
                    {
                        possibilities[possibilities_count] = v;
                        possibilities_count++;
                    }
                }

                for (i=initial_vertices_count;i<G->n;i++)
                {
                    v = generator.random(0, possibilities_count-1);
                    current_permutation[i] = possibilities[v];
                    possibilities[v] = possibilities[possibilities_count-1];
                    possibilities_count--;
                }
            }

            greedy_obdom(G, current_permutation, result, connected);
            f1_o = problemdom->undominated_vertices_count(G, result);
            f2_o = problemdom->dominating_set_size(G, result);
            f1_bc = f1_o;
            f2_bc = f2_o;

            t_stag = 0;
            cycles++;
            extended_waiting_time = false;
        }


        // new permutation with jump operator
        for (i=0;i<G->n;i++)
        {
            new_permutation[i] = current_permutation[i];
        }
        q = generator.random(1, G->n-1);
        jumped_value = new_permutation[q];
        for (i=q;i>0;i--)
        {
            new_permutation[i] = new_permutation[i-1];
        }
        new_permutation[0] = jumped_value;

        greedy_obdom(G, new_permutation, result, connected);
        f1_n = problemdom->undominated_vertices_count(G, result);
        f2_n = problemdom->dominating_set_size(G, result);

        if (f1_n < f1_o || (f1_n == f1_o && f2_n <= f2_o))
        {
            // output only if the best has been improved
            if (f1_n < f1_b || (f1_n == f1_b && f2_n < f2_b))
            {
                printf("[%u,%u]\n", f1_n, f2_n);
                f1_b = f1_n;
                f2_b = f2_n;
                for (i=0;i<G->n;i++)
                {
                    best_permutation[i] = new_permutation[i];
                }
                // if the best solution was hit, then we extend the search this time
                extended_waiting_time = true;
            }
            // stagnation checking against the current run
            if (f1_n < f1_bc || (f1_n == f1_bc && f2_n < f2_bc))
            {
                f1_bc = f1_n;
                f2_bc = f2_n;
                t_stag = 0;
            }
            f1_o = f1_n;
            f2_o = f2_n;
            for (i=0;i<G->n;i++)
            {
                current_permutation[i] = new_permutation[i];
            }
        }

        t++;
        t_stag++;
    }

    greedy_obdom(G, best_permutation, result, connected);
    *output_t = t;

    delete[](best_permutation);
}
