#include "cli.h"
#include "algorithm_greedydom.h"
#include "algorithm_obdom.h"
#include "algorithm_acols.h"

cli::cli()
{
    greedy_max_overall = 0.0;
    greedy_min_overall = 0.0;
    greedy_avg_overall = 0.0;
    undominated_vertices_overall_rls = 0.0;
    dominating_set_size_overall_rls = 0.0;
    undominated_vertices_overall_ma = 0.0;
    dominating_set_size_overall_ma = 0.0;
    undominated_vertices_overall_msob = 0.0;
    dominating_set_size_overall_msob = 0.0;
    undominated_vertices_overall_ob = 0.0;
    dominating_set_size_overall_ob = 0.0;
    undominated_vertices_overall_obma = 0.0;
    dominating_set_size_overall_obma = 0.0;
    undominated_vertices_overall_acols = 0.0;
    dominating_set_size_overall_acols = 0.0;
    undominated_vertices_overall_acoppls = 0.0;
    dominating_set_size_overall_acoppls = 0.0;
    undominated_vertices_overall_acols_s = 0.0;
    dominating_set_size_overall_acols_s = 0.0;
    t_acols = 0;
    t_acoppls = 0;
    t_acols_s = 0;
    t_obdom = 0;
}

void cli::sleep(int milisec)
{
    QTime timer_for_sleep;

    timer_for_sleep.start();

    while (timer_for_sleep.elapsed() < milisec);
}

void cli::component_labeling()
{
    refer v, j, v_central = 0, queue_size, queue_current, current_label;

    refer *queue = new refer[G->n];
    bool *processed = new bool[G->n];

    for (v=0;v<G->n;v++)
    {
        processed[v] = 0;
    }

    current_label = 0;
    for (v_central=0;v_central<G->n;v_central++)
    {
        if (processed[v_central])
        {
            continue;
        }
        // breadth-first search
        queue_current = 0;
        queue_size = 1;
        queue[0] = v_central;
        while (queue_size)
        {
            v = queue[queue_current];
            queue_current++;
            queue_size--;
            processed[v] = 1;
            component_labels[v] = current_label;
            for (j=0;j<G->V[v].edgecount;j++)
            {
                if (! processed[G->V[v].sibl[j]])
                {
                    queue[queue_current+queue_size] = G->V[v].sibl[j];
                    queue_size++;
                    processed[G->V[v].sibl[j]] = 1;
                    component_labels[G->V[v].sibl[j]] = current_label;
                }
            }
        }
        current_label++;
    }

    delete[](queue);
    delete[](processed);

    components = current_label;
}

void cli::exhaustive_search_core(refer current_vertex, refer depth, refer current_size)
{
    refer v;

    if (depth > current_size || found_optimum_exhaustive)
    {
        return;
    }

    if (depth == current_size)
    {
        //if (0 == current_nondominated_vertices)
        if (0 == problemdom_exhaustive->undominated_vertices_count_component(G, current_solution_exhaustive, current_component_label_exhaustive, component_labels))
        {
            printf("Found dominating set of size %ld.\n", depth);
            found_optimum_exhaustive = true;
            current_component_dominating_set_size = depth;
            return;
        }
        scanned++;
        if (0 == scanned % 1000)
        {
            printf("Solutions scanned: %u.\n", scanned);
        }
    }

    for (v=current_vertex;v<G->n;v++)
    {
        if (found_optimum_exhaustive)
        {
            return;
        }
        if (0 == current_solution_exhaustive[v] && current_component_label_exhaustive == component_labels[v])
        {
            current_solution_exhaustive[v] = 1;
            // updating the values
            exhaustive_search_core(v+1, depth+1, current_size);
            current_solution_exhaustive[v] = 0;
            // reverting the values back
        }
    }
}

void cli::exhaustive_search(refer current_size)
{
    refer v;
    refer current_vertex, next_vertex;
    refer depth, current_component_size;
    bool new_vertices_found;

    found_optimum_exhaustive = false;

    problemdom_exhaustive = new problem();

    current_nondominated_vertices = 0;
    current_component_size = 0;
    for (v=0;v<G->n;v++)
    {
        current_solution_exhaustive[v] = 0;
        if (current_component_label_exhaustive == component_labels[v])
        {
            current_component_size++;
        }
    }

    current_vertex_stack[0] = 0;
    current_vertex_stack_size = 1;
    depth_stack[0] = 0;
    depth_stack_size = 1;
    next_vertex_stack[0] = 0;
    next_vertex_stack_size = 1;
    new_vertices_found = false;

    scanned = 0;

    do
    {
        current_vertex = current_vertex_stack[current_vertex_stack_size-1];
        depth = depth_stack[depth_stack_size-1];
        next_vertex = next_vertex_stack[next_vertex_stack_size-1];

        if (depth == current_size)
        {
            //if ((depth < 10 && 0 == problemdom_exhaustive->undominated_vertices_count_component_alternative(G, &current_vertex_stack[1], depth, current_component_size))
            //    || (depth >= 10 && 0 == problemdom_exhaustive->undominated_vertices_count_component(G, current_solution_exhaustive, current_component_label_exhaustive, component_labels)))
            if (0 == problemdom_exhaustive->undominated_vertices_count_component_alternative(G, &current_vertex_stack[1], depth, current_component_size))
            {
                printf("Found dominating set of size %ld.\n", depth);
                current_component_dominating_set_size = depth;
                found_optimum_exhaustive = true;
                break;
            }
            scanned++;
            if (0 == scanned % 10000)
            {
                printf("%u solutions scanned.\n", scanned);
            }
        }

        new_vertices_found = false;
        if (depth < current_size)
        {
            for (v=next_vertex;v<G->n;v++)
            {
                if (0 == current_solution_exhaustive[v] && current_component_label_exhaustive == component_labels[v])
                {
                    current_solution_exhaustive[v] = 1;
                    current_vertex_stack[current_vertex_stack_size] = v;
                    depth_stack[depth_stack_size] = depth+1;
                    next_vertex_stack[next_vertex_stack_size-1] = v+1;
                    next_vertex_stack[next_vertex_stack_size] = v;
                    current_vertex_stack_size++;
                    depth_stack_size++;
                    next_vertex_stack_size++;
                    new_vertices_found = true;
                    break;
                }
            }
        }

        if (! new_vertices_found)
        {
            current_solution_exhaustive[current_vertex] = 0;
            current_vertex_stack_size--;
            depth_stack_size--;
            next_vertex_stack_size--;
        }
    }
    while (0 < next_vertex_stack_size);

    delete(problemdom_exhaustive);
}

int cli::process_from_file()
{
    // instance loaded from file
    choose_instance();
    // shortcut graph
    /*printf("k = ");
    if (EOF == scanf("%u", &k))
    {
        //break;
    }*/
    //generate_shortcut_graph(G, k);
    compute_statistics();

    // full MDS
    //domset_to_mps(false);
    // LP relaxation
    //domset_to_mps(true);

    choose_algorithm();

    printf("===========================================================\n");
    printf("min greedy: %0.2lf\n", greedy_min_overall / (double) MAX_INSTANCE);
    printf("max greedy: %0.2lf\n", greedy_max_overall / (double) MAX_INSTANCE);
    printf("average greedy: %0.2lf\n", greedy_avg_overall / (double) MAX_INSTANCE);
    printf("total order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_ob / (double) MAX_INSTANCE, dominating_set_size_overall_ob / (double) MAX_INSTANCE, (long) (t_obdom / MAX_INSTANCE));
    //printf("total multi-start order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_msob / (double) MAX_INSTANCE, dominating_set_size_overall_msob / (double) MAX_INSTANCE, (long) (t_msobdom / MAX_INSTANCE));
    printf("total ACO-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acols / (double) MAX_INSTANCE, dominating_set_size_overall_acols / (double) MAX_INSTANCE, (long) (t_acols / MAX_INSTANCE));
    printf("total ACO-PP-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acoppls / (double) MAX_INSTANCE, dominating_set_size_overall_acoppls / (double) MAX_INSTANCE, (long) (t_acoppls / MAX_INSTANCE));
    printf("total ACO-LS-S: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acols_s / (double) MAX_INSTANCE, dominating_set_size_overall_acols_s / (double) MAX_INSTANCE, (long) (t_acols_s / MAX_INSTANCE));

    return 0;
}

int cli::process_ba()
{
    int j;

    // instance generated by BA algorithm

    printf("w = "); scanf("%lu",&w);
    printf("n_max = "); scanf("%lu",&n_max);

    for (j=0;j<MAX_INSTANCE;j++)
    {
        generate_instance();
        compute_statistics();
        choose_algorithm();
    }

    printf("===========================================================\n");
    printf("min greedy: %0.2lf\n", greedy_min_overall / (double) MAX_INSTANCE);
    printf("max greedy: %0.2lf\n", greedy_max_overall / (double) MAX_INSTANCE);
    printf("average greedy: %0.2lf\n", greedy_avg_overall / (double) MAX_INSTANCE);
    printf("total order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_ob / (double) MAX_INSTANCE, dominating_set_size_overall_ob / (double) MAX_INSTANCE, (long) (t_obdom / MAX_INSTANCE));
    //printf("total multi-start order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_msob / (double) MAX_INSTANCE, dominating_set_size_overall_msob / (double) MAX_INSTANCE, (long) (t_msobdom / MAX_INSTANCE));
    printf("total ACO-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acols / (double) MAX_INSTANCE, dominating_set_size_overall_acols / (double) MAX_INSTANCE, (long) (t_acols / MAX_INSTANCE));
    printf("total ACO-PP-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acoppls / (double) MAX_INSTANCE, dominating_set_size_overall_acoppls / (double) MAX_INSTANCE, (long) (t_acoppls / MAX_INSTANCE));
    printf("total ACO-LS-S: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acols_s / (double) MAX_INSTANCE, dominating_set_size_overall_acols_s / (double) MAX_INSTANCE, (long) (t_acols_s / MAX_INSTANCE));

    return 0;
}

int cli::process_udg()
{
    int j;

    // instance generated by UDG algorithm

    printf("n_max = "); scanf("%lu",&n_max);
    printf("range = "); scanf("%lu",&range);
    printf("grid = "); scanf("%lu",&grid);

    for (j=0;j<MAX_INSTANCE;j++)
    {
        generate_instance_udg();
        compute_statistics();
        choose_algorithm();
    }

    printf("===========================================================\n");
    printf("min greedy: %0.2lf\n", greedy_min_overall / (double) MAX_INSTANCE);
    printf("max greedy: %0.2lf\n", greedy_max_overall / (double) MAX_INSTANCE);
    printf("average greedy: %0.2lf\n", greedy_avg_overall / (double) MAX_INSTANCE);
    printf("total order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_ob / (double) MAX_INSTANCE, dominating_set_size_overall_ob / (double) MAX_INSTANCE, (long) (t_obdom / MAX_INSTANCE));
    //printf("total multi-start order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_msob / (double) MAX_INSTANCE, dominating_set_size_overall_msob / (double) MAX_INSTANCE, (long) (t_msobdom / MAX_INSTANCE));
    printf("total ACO-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acols / (double) MAX_INSTANCE, dominating_set_size_overall_acols / (double) MAX_INSTANCE, (long) (t_acols / MAX_INSTANCE));
    printf("total ACO-PP-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acoppls / (double) MAX_INSTANCE, dominating_set_size_overall_acoppls / (double) MAX_INSTANCE, (long) (t_acoppls / MAX_INSTANCE));
    printf("total ACO-LS-S: [%0.2lf,%0.2lf] (iter: %ld)\n", undominated_vertices_overall_acols_s / (double) MAX_INSTANCE, dominating_set_size_overall_acols_s / (double) MAX_INSTANCE, (long) (t_acols_s / MAX_INSTANCE));

    return 0;
}

int cli::start_cli()
{
    int j;

    srand((unsigned long long) time(NULL));

    printf("DOMINATING SET PROBLEM SOLVER\nversion 0.1\n");

    return process_from_file();

    //return process_ba();

    //return process_udg();
}

int cli::generate_instance()
{
    if (get_graph() != NULL)
    {
        free_graph();
    }
    generate_graph_BA_model(w,n_max);
    G = get_graph();

    return 0;
}

int cli::generate_instance_udg()
{
    if (get_graph() != NULL)
    {
        free_graph();
    }
    generate_graph_UDG(n_max,range,grid);
    G = get_graph();

    return 0;
}

int cli::choose_instance()
{
    printf("File:");
    scanf("%s",&filename[0]);

    printf("Loading the data...\n");
    if ((source = fopen(filename,"r")) == NULL)
    {
        printf("Error: An error occured during opening of the file.");
        getchar();
        return 1;
    }
    if (get_graph() != NULL)
    {
        free_graph();
    }
    input_graph(source);
    fclose(source);
    G = get_graph();

    return 0;
}

int cli::compute_statistics()
{
    refer i, maxdeg, diameter;
    unsigned long long triangles;

    printf("Computing the basic statistics...\n");

    min_deg = statistics::min_degree(G);
    max_deg = statistics::max_degree(G);
    avg_deg = statistics::average_degree(G);
    stdev_deg = statistics::degree_stdev(G);
    //triangles = statistics::triangles(G);
    //diameter = statistics::diameter(G);
    lower_bound = statistics::domset_lb(G);

    printf("G: %s\n",filename);
    printf("|V| = %u, |E| = %lu, d = %0.3lf\n",G->n,G->m,(double)(2*G->m)/(double)(G->n)/(double)(G->n-1));
    printf("min degree = %ld, max degree = %ld\n",min_deg,max_deg);
    printf("average degree = %0.3lf, degree stdev = %0.3lf\n",avg_deg,stdev_deg);
    printf("componets = %u\n",statistics::components(G));
    //printf("largest component diameter = %u\n", diameter);
    //printf("mean cc = %0.3lf\n",statistics::mean_clustering_coefficient(G));
    //printf("number of triangles = %llu\n",triangles);
    printf("%u <= dominating set size\n",lower_bound);

    /*printf("degree distribution: [");
    maxdeg = statistics::degree_distribution(G, degree_distrib);
    for (i=0;i<=maxdeg;i++)
    {
        printf("%ld", degree_distrib[i]);
        if (i != maxdeg)
        {
            putchar(',');
        }
    }
    printf("]\n");*/

    return 0;
}

int cli::choose_algorithm()
{

    printf("Here will come the run of our dominating set algorithm...\n");

    problem *problemdom = new problem();
    algorithm_greedydom *greedydom = new algorithm_greedydom();
    algorithm_obdom *obdom = new algorithm_obdom();
    algorithm_acols *acols = new algorithm_acols();

    //long max_time = 1 * 1000;
    //long max_time = 60 * 1000;
    //long max_time = 3 * 60 * 1000;
    //long max_time = 5 * 60 * 1000;
    //long max_time = 10 * 60 * 1000;
    long max_time = 60 * 60 * 1000;

    timer.start();

    refer mingreedy = G->n;
    refer maxgreedy = 0;
    double avggreedy = 0.0;
    refer undominated_vertices = 0;
    refer dominating_set_size = G->n*100;
    refer undominated_vertices_msob = 0;
    refer dominating_set_size_msob = G->n*100;
    refer undominated_vertices_ob = 0;
    refer dominating_set_size_ob = G->n*100;
    refer undominated_vertices_acols = 0;
    refer dominating_set_size_acols = G->n*100;
    refer undominated_vertices_acoppls = 0;
    refer dominating_set_size_acoppls = G->n*100;
    refer undominated_vertices_acols_s = 0;
    refer dominating_set_size_acols_s = G->n*100;

    double dominating_set_size_average_msob = 0.0;
    double dominating_set_size_average_ob = 0.0;
    double dominating_set_size_average_acols = 0.0;
    double dominating_set_size_average_acoppls = 0.0;
    double dominating_set_size_average_acols_s = 0.0;

    long long partial_t_obdom = 0;
    long long partial_t_msobdom = 0;
    long long partial_t_acols = 0;
    long long partial_t_acoppls = 0;
    long long partial_t_acols_s = 0;
    refer second_lower_bound;

    refer own_lower_bound;
    printf("own lower bound = ");
    scanf("%u", &own_lower_bound);
    //own_lower_bound = 0;

    // GREEDY
    /*printf("Greedy...\n");
    refer greedy_count = 1000;
    for (refer j=0;j<greedy_count;j++)
    {
        // greedy dominating set
        greedydom->greedydom(G, result);
        undominated_vertices = problemdom->undominated_vertices_count(G, result);
        dominating_set_size = problemdom->dominating_set_size(G, result);
        //printf("[%ld,%ld]\n", undominated_vertices, dominating_set_size);
        if (undominated_vertices == 0 && dominating_set_size < mingreedy)
        {
            for (i=0;i<G->n;i++)
            {
                best_result[i] = result[i];
            }
            mingreedy = dominating_set_size;
        }
        if (undominated_vertices == 0 && dominating_set_size > maxgreedy)
        {
            maxgreedy = dominating_set_size;
        }
        avggreedy += (double) dominating_set_size;

        sleep(10);
    }
    avggreedy /= (double) greedy_count;

    second_lower_bound = ceil((double) maxgreedy / (double) (log(max_deg) + 1));

    printf("min greedy: %u\n", mingreedy);
    printf("max greedy: %u\n", maxgreedy);
    printf("%u <= dominating set size\n",second_lower_bound);*/

    for (int run=0;run<MAX_ALGORITHM;run++)
    {
        // RLS_o
        printf("Order-based RLS:\n");
        for (refer j=0;j<1;j++)
        {
            //obdom->obdom(G, result, max_time, (lower_bound > second_lower_bound) ? lower_bound : second_lower_bound, &partial_t_obdom);
            obdom->obdom(G, result, max_time, own_lower_bound, &partial_t_obdom);
        }
        printf("\n");

        undominated_vertices_ob = problemdom->undominated_vertices_count(G, result);
        if (dominating_set_size_ob >= problemdom->dominating_set_size(G, result))
        {
            dominating_set_size_ob = problemdom->dominating_set_size(G, result);
        }
        t_obdom += partial_t_obdom;
        dominating_set_size_average_ob += (double) problemdom->dominating_set_size(G, result);

        // MSRLS_o
        /*printf("Multi-start order-based RLS:\n");
        for (refer j=0;j<1;j++)
        {
            obdom->obdom(G, result, max_time, own_lower_bound, &partial_t_msobdom, true, true);
            //obdom->obdom(G, result, max_time, own_lower_bound, &partial_t_msobdom, true, false, false, true);
        }
        printf("\n");*/

        undominated_vertices_msob = problemdom->undominated_vertices_count(G, result);
        if (dominating_set_size_msob >= problemdom->dominating_set_size(G, result))
        {
            dominating_set_size_msob = problemdom->dominating_set_size(G, result);
        }
        t_msobdom += partial_t_msobdom;
        dominating_set_size_average_msob += (double) problemdom->dominating_set_size(G, result);

        // ACO-LS original
        printf("ACO-LS:\n");
        acols->acols(G, result, max_time, 10000000, 20, false, false, own_lower_bound, false, 1.0, 10.0, &partial_t_acols);

        undominated_vertices_acols = problemdom->undominated_vertices_count(G, result);
        if (dominating_set_size_acols >= problemdom->dominating_set_size(G, result))
        {
            dominating_set_size_acols = problemdom->dominating_set_size(G, result);
        }
        t_acols += partial_t_acols;
        dominating_set_size_average_acols += (double) problemdom->dominating_set_size(G, result);

        // ACO-PP-LS original
        printf("ACO-PP-LS:\n");
        acols->acols(G, result, max_time, 10000000, 20, false, true, own_lower_bound, false, 2.0, 5.0, &partial_t_acoppls);

        undominated_vertices_acoppls = problemdom->undominated_vertices_count(G, result);
        if (dominating_set_size_acoppls >= problemdom->dominating_set_size(G, result))
        {
            dominating_set_size_acoppls = problemdom->dominating_set_size(G, result);
        }
        t_acoppls += partial_t_acoppls;
        dominating_set_size_average_acoppls += (double) problemdom->dominating_set_size(G, result);

        // ACO-LS sparse
        printf("ACO-LS-S:\n");
        acols->acols(G, result, max_time, 10000000, 20, true, false, own_lower_bound, true, 1.0, 10.0, &partial_t_acols_s);

        undominated_vertices_acols_s = problemdom->undominated_vertices_count(G, result);
        if (dominating_set_size_acols_s >= problemdom->dominating_set_size(G, result))
        {
            dominating_set_size_acols_s = problemdom->dominating_set_size(G, result);
        }
        t_acols_s += partial_t_acols_s;
        dominating_set_size_average_acols_s += (double) problemdom->dominating_set_size(G, result);

        // output - temporary
        /*bool out[MAX_VERTICES],out2[MAX_VERTICES];
        for (i=0;i<G->n;i++) out[i] = false;
        for (i=0;i<G->n;i++) out2[i] = false;
        for (i=0;i<G->n;i++)
        {
            if (result[i] == 1) out[i] = true;
        }
        FILE *f = fopen("./solution.txt", "w");
        for (i=0;i<G->n;i++)
        {
            if (out[i]) { fprintf(f, "%u ", i+1);
            for (refer j=0;j<G->V[i].edgecount;j++)
                if (! out2[G->V[i].sibl[j]] && ! out[G->V[i].sibl[j]])
                { fprintf(f, "%u ", G->V[i].sibl[j]+1);
                    out2[G->V[i].sibl[j]] = true;
                }
            fprintf(f, "\n");}
        }
        fclose(f);*/

        /*printf("%ld <= dominating set size\n",second_lower_bound);
        printf("%ld <= dominating set size\n",statistics::domset_lb(G));*/
    }

    printf("For this instance:\n");
    printf("order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", 0.0, dominating_set_size_average_ob / (double) MAX_ALGORITHM, (long) (t_obdom / MAX_ALGORITHM));
    //printf("multi-start order-based LS: [%0.2lf,%0.2lf] (iter: %ld)\n", 0.0, dominating_set_size_average_msob / (double) MAX_ALGORITHM, (long) (t_msobdom / MAX_ALGORITHM));
    printf("ACO-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", 0.0, dominating_set_size_average_acols / (double) MAX_ALGORITHM, (long) (t_acols / MAX_ALGORITHM));
    printf("ACO-PP-LS: [%0.2lf,%0.2lf] (iter: %ld)\n", 0.0, dominating_set_size_average_acoppls / (double) MAX_ALGORITHM, (long) (t_acoppls / MAX_ALGORITHM));
    printf("ACO-LS-S: [%0.2lf,%0.2lf] (iter: %ld)\n", 0.0, dominating_set_size_average_acols_s / (double) MAX_ALGORITHM, (long) (t_acols_s / MAX_ALGORITHM));

    greedy_min_overall += (double) mingreedy;
    greedy_max_overall += (double) maxgreedy;
    greedy_avg_overall += (double) avggreedy;
    undominated_vertices_overall_msob += undominated_vertices_msob;
    dominating_set_size_overall_msob += dominating_set_size_msob;
    undominated_vertices_overall_ob += undominated_vertices_ob;
    dominating_set_size_overall_ob += dominating_set_size_ob;
    undominated_vertices_overall_acols += undominated_vertices_acols;
    dominating_set_size_overall_acols += dominating_set_size_acols;
    undominated_vertices_overall_acoppls += undominated_vertices_acoppls;
    dominating_set_size_overall_acoppls += dominating_set_size_acoppls;
    undominated_vertices_overall_acols_s += undominated_vertices_acols_s;
    dominating_set_size_overall_acols_s += dominating_set_size_acols_s;


    delete(greedydom);
    delete(obdom);
    delete(acols);
    delete(problemdom);

    free_graph();

    return 0;
}

// converts the dominating set instance to the MPS ILP format
void cli::domset_to_mps(bool relaxation)
{
    FILE *mpsfile;
    refer v, j, i;
    bool *domination_influences = new bool[G->n];
    char filename_mps[PATH_MAX];

    sprintf(filename_mps, "%s_%u.mps", filename, k);

    mpsfile = fopen(filename_mps, "w");
    fprintf(mpsfile, "NAME DOMINATINGSET\n");
    fprintf(mpsfile, "ROWS\n");
    fprintf(mpsfile, " N COST\n");
    for (v=0;v<G->n;v++)
    {
        fprintf(mpsfile, "    G LIM%u\n", v+1);
    }
    fprintf(mpsfile, "COLUMNS\n");
    for (v=0;v<G->n;v++)
    {
        fprintf(mpsfile, "    x%u   COST 1\n", v+1);
        for (j=0;j<G->n;j++)
        {
            domination_influences[j] = false;
        }
        domination_influences[v] = true;
        // other constraints - neighbours of v
        for (j=0;j<G->V[v].edgecount;j++)
        {
            domination_influences[G->V[v].sibl[j]] = true;
        }
        for (j=0;j<G->n;j++)
        {
            if (domination_influences[j])
            {
                fprintf(mpsfile, "    x%u   LIM%u %u\n", v+1, j+1, domination_influences[j]);
            }
        }
    }
    fprintf(mpsfile, "RHS\n");
    for (v=0;v<G->n;v++)
    {
        fprintf(mpsfile, "    RHS1 LIM%u 1\n", v+1);
    }
    if (! relaxation)
    {
        fprintf(mpsfile, "BOUNDS\n");
        for (v=0;v<G->n;v++)
        {
            fprintf(mpsfile, " BV ONE x%u 1\n", v+1, v+1);
        }
    }
    fprintf(mpsfile, "ENDATA\n");
    fclose(mpsfile);

    delete[](domination_influences);
}
