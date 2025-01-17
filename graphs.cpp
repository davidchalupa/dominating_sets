/*
        ===================
        IMPLEMENTACIA GRAFU
        ===================

        G = [V,E]
        |V| = n
        |E| = m
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "graphs.h"
#include "edgetable.h"
#include "random_generator.h"

#define MAX_BA_DEGREE 1000
#define MAX_VERTEX_LABEL 50

graph G;

// auxiliary arrays for edges (to need only one malloc)
refer *E1;
refer *E2;
// auxiliary array for the degrees
refer VEc[MAX_VERTICES];
refer weights[MAX_VERTICES];
random_generator generator;
unsigned long roulette[MAX_VERTICES];
unsigned long roulette_sum;
edgetable *edges;

bool is_labeled;
char vertex_labels[MAX_VERTICES][MAX_VERTEX_LABEL+1];

graph get_graph()
{
    return G;
}

int input_graph(FILE *source)
{
    char ch,chb;
    unsigned long i, j;
    refer v,w;
    long long cnt = 0;
    bool weighted = false;

    // comments
    is_labeled = false;
    while ((ch = fgetc(source)) == 'c')
    {
        chb = fgetc(source);

        if ((ch = fgetc(source)) == 'w')
        {
            weighted = true;
            fscanf(source, "%u", &v);
            fscanf(source, "%u", &w);
            weights[v-1] = w;
            ch = fgetc(source);
        }
        else
        {
            ungetc(ch, source);
            ungetc(chb, source);
            if (EOF != fscanf(source, "%u", &v) && (is_labeled || (! is_labeled && 1 == v)))
            {
                fgets(&vertex_labels[v-1][0], MAX_VERTEX_LABEL, source);
                is_labeled = true;
                if ('\n' != vertex_labels[v-1][strlen(vertex_labels[v-1])-1])
                {
                    while ((ch = fgetc(source)) != '\n');
                }
            }
            else
            {
                while ((ch = fgetc(source)) != '\n');
            }
        }
    }

    // overriding the string "p edge"
    for (i=0;i<1;i++)
    {
        fgetc(source);
    }
    for (i=2;i<6;i++)
    {
        fgetc(source);
    }

    G = (graph) malloc (sizeof(graph_data));
    G->weighted = weighted;

    if ((fscanf(source,"%u",&G->n)) == EOF)
    {
        return 1;
    }
    if ((fscanf(source,"%lu",&G->m)) == EOF)
    {
        return 1;
    }

    // allocation of the auxiliary space
    E1 = new refer[2*G->m+1];
    E2 = new refer[2*G->m+1];

    edges = new edgetable(2*G->m+1);

    // separator line
    fgetc(source);

    // loading the edges to auxiliary arrays
    for (i=0;i<G->m;i++)
    {
        // new line and letter 'e'
        fgetc(source);
        fgetc(source);
        if (fscanf(source,"%u",&v) == EOF)
        {
            break;
        }
        if (fscanf(source,"%u",&w) == EOF)
        {
            break;
        }
        // we accept only one "direction"
        v--; w--;
        if (v != w)
        {
            E1[cnt] = v; E2[cnt] = w;
            VEc[v]++; cnt++;
            E1[cnt] = w; E2[cnt] = v;
            VEc[w]++; cnt++;
            edges->insert(v,w);
            edges->insert(w,v);
        }
        fgetc(source);
    }

    // allocation based on the computed degrees
    for (i=0;i<G->n;i++)
    {
        G->V[i].edgecount = 0;
        G->V[i].sibl = (refer *) malloc (VEc[i]*sizeof(refer));
        G->V[i].weight = weights[i];
    }
    // insertion of data
    unsigned long real_m = 0;
    for (j=0;j<cnt;j++)
    {
        v = E1[j]; w = E2[j];
        G->V[v].edgecount++;
        G->V[v].sibl[G->V[v].edgecount-1] = w;
        real_m++;
    }

    // adjacency lists sorting
    for (i=0;i<G->n;i++)
    {
        QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
    }
    // getting rid of reduntant edges
    bool changed;
    refer decrement;
    for (i=0;i<G->n;i++)
    {
        changed = false;
        decrement = 0;
        if (0 < G->V[i].edgecount)
        {
            for (j=0;j<G->V[i].edgecount-1;j++)
            {
                if (G->V[i].sibl[j] == G->V[i].sibl[j+1])
                {
                    // using a fake vertex to get rid of this
                    // quicksort will then put it to the end
                    G->V[i].sibl[j+1] = G->n+1;
                    decrement++;
                    real_m--;
                    changed = true;
                }
            }
        }
        // sorting once again if needed
        if (changed)
        {
            QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
        }
        // cutting off the tail
        G->V[i].edgecount -= decrement;
        for (j=G->V[i].edgecount;j<G->V[i].edgecount+decrement;j++)
        {
            G->V[i].sibl[j] = 0;
        }
    }

    // the real number of edges cleaned of redundant occurrences
    G->m = real_m / 2;
    G->density = ((double)(G->m))/((double)(G->n*(G->n-1)/2));

    delete[](E1);
    delete[](E2);

    return 0;
}

int output_graph(FILE *target)
{
    refer v,j;

    if (is_labeled)
    {
        for (v=0;v<G->n;v++)
        {
            fprintf(target, "c %u %s", v+1, vertex_labels[v]);
        }
    }
    fprintf(target, "p edge %u %lu\n", G->n, G->m);
    for (v=0;v<G->n;v++)
    {
        for (j=0;j<G->V[v].edgecount;j++)
        {
            if (v <= G->V[v].sibl[j])
            {
                fprintf(target, "e %u %u\n", v+1, G->V[v].sibl[j]+1);
            }
        }
    }

    return 0;
}

// dealokujeme potrebnu pamat
void free_graph()
{
    unsigned long i;

    for (i=0;i<G->n;i++)
    {
        G->V[i].edgecount = 0;
        free(G->V[i].sibl);
    }

    free(G);
    delete(edges);
    G = NULL;
}

bool are_adjacent_binary_search(refer v, refer w)
{
    // the adjacency lists are sorted, so we can use binary search
    if (-1 != BinarySearch(G->V[v].sibl, w, 0, G->V[v].edgecount-1))
    {
        return true;
    }
    if (-1 != BinarySearch(G->V[w].sibl, v, 0, G->V[w].edgecount-1))
    {
        return true;
    }
    return false;
}

bool are_adjacent_edgetable(refer v, refer w)
{
    return edges->isin(v,w);
}

bool are_adjacent(refer v, refer w)
{
    return are_adjacent_edgetable(v,w);
}

void generate_graph_BA_model(unsigned long w, unsigned long n_max)
{
    unsigned long i;
    unsigned long j,q,r;

    srand((unsigned long long) time(NULL));

    // alokujeme priestor pre graf
    G = (graph) malloc (sizeof(graph_data));

    // inicializujeme zakladny graf
    G->n = w;
    G->m = w-1;
    if (0 == w)
    {
        G->m = 0;
    }
    else
    {
        G->m = w-1;
    }
    for (i=0;i<n_max;i++)
    {
        G->V[i].sibl = (refer *) malloc (MAX_BA_DEGREE*sizeof(refer));
    }
    if (w <= 1)
    {
        G->V[0].edgecount = 0;
    }
    else
    {
        for (i=0;i<w;i++)
        {
            if (i == 0 || i == w-1)
            {
                G->V[i].edgecount = 1;
                if (i == 0)
                {
                    G->V[i].sibl[0] = i+1;
                }
                if (i == w - 1)
                {
                    G->V[i].sibl[0] = i-1;
                }
            }
            else
            {
                G->V[i].edgecount = 2;
                G->V[i].sibl[0] = i-1;
                G->V[i].sibl[1] = i+1;
            }
        }
    }

    refer n = G->n;
    for (i=n;i<n_max;i++)
    {
        G->V[i].edgecount = w;
        // roulette wheel preparation
        roulette[0] = 0;
        roulette_sum = 0;
        for (j=1;j<i;j++)
        {
            roulette[j] = roulette[j-1] + G->V[j-1].edgecount;
            roulette_sum += G->V[j-1].edgecount;
        }
        roulette_sum += G->V[i-1].edgecount;
        // roulette wheel selection of new vertices
        for (j=0;j<w;j++)
        {
            do
            {
                r = generator.random(0,roulette_sum-1);
                q = 0;
                while (! (roulette[q] <= r && (q >= i-1 || roulette[q+1] > r)))
                {
                    q++;
                }
            }
            while (are_adjacent(i,q) || are_adjacent(q,i));
            G->V[i].sibl[j] = q;
            if (0 == G->V[q].edgecount % MAX_BA_DEGREE)
            {
                G->V[q].sibl = (refer *) realloc (G->V[q].sibl, (G->V[q].edgecount+1)*MAX_BA_DEGREE*sizeof(refer));
            }
            G->V[q].sibl[G->V[q].edgecount] = i;
            G->V[q].edgecount++;
        }
        G->n++;
        G->m += w;
    }

    edges = new edgetable(2*G->m+1);
    for (i=0;i<G->n;i++)
    {
        for (j=0;j<G->V[i].edgecount;j++)
        {
            edges->insert(i,G->V[i].sibl[j]);
        }
    }

    for (i=0;i<G->n;i++)
    {
        QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
        if (G->V[i].edgecount > MAX_BA_DEGREE)
        {
            // ToDo: overflow handling!
        }
    }
}

// generates a graph with direct shortcuts if distance between two vertices is at most k
void generate_shortcut_graph(graph G, refer k)
{
    // this makes sense only for k >= 2
    if (k < 2)
    {
        return;
    }

    refer v, i, j, v_central = 0, queue_size, queue_current, current_distance;
    refer new_degree, new_vertices_count;

    refer *distances = new refer[G->n];
    refer *queue = new refer[G->n];
    refer *new_vertices = new refer[G->n];
    refer *new_vertex_degrees = new refer[G->n];

    for (v_central=0;v_central<G->n;v_central++)
    {
        // determining the new degree
        new_vertices_count = 0;
        new_degree = G->V[v_central].edgecount;
        for (v=0;v<G->n;v++)
        {
            distances[v] = 0;
        }
        // breadth-first search
        queue_current = 0;
        queue_size = 1;
        queue[0] = v_central;
        while (queue_size)
        {
            v = queue[queue_current];
            current_distance = distances[v];
            queue_current++;
            queue_size--;
            // we cut the path off if we already are at distance k
            if (current_distance >= k)
            {
                continue;
            }
            for (j=0;j<G->V[v].edgecount;j++)
            {
                if (0 == distances[G->V[v].sibl[j]] && G->V[v].sibl[j] != v_central)
                {
                    queue[queue_current+queue_size] = G->V[v].sibl[j];
                    distances[G->V[v].sibl[j]] = current_distance + 1;
                    queue_size++;
                    if (1 < distances[G->V[v].sibl[j]] && distances[G->V[v].sibl[j]] <= k)
                    {
                        new_degree++;
                        new_vertices[new_vertices_count] = G->V[v].sibl[j];
                        new_vertices_count++;
                    }
                }
            }
        }

        // reallocation
        G->V[v_central].sibl = (refer *) realloc (G->V[v_central].sibl, new_degree*sizeof(refer));

        // putting the new vertices in adjacency lists but still not recognizing them
        for (j=0;j<new_vertices_count;j++)
        {
            v = new_vertices[j];
            G->V[v_central].sibl[G->V[v_central].edgecount+j] = v;
        }

        new_vertex_degrees[v_central] = G->V[v_central].edgecount + new_vertices_count;
    }

    for (v_central=0;v_central<G->n;v_central++)
    {
        // recognizing the new vertices by updating the edge count
        G->V[v_central].edgecount = new_vertex_degrees[v_central];

        // sorting the adjacency list
        QuickSort(G->V[v_central].sibl,0,G->V[v_central].edgecount-1);
    }

    // new number of edges
    G->m = 0;
    for (v=0;v<G->n;v++)
    {
        G->m += G->V[v].edgecount;
    }
    G->m /= 2;
    G->density = ((double)(G->m))/((double)(G->n*(G->n-1)/2));

    edges->set_max_edges(2*G->m+1);
    for (i=0;i<G->n;i++)
    {
        for (j=0;j<G->V[i].edgecount;j++)
        {
            edges->insert(i,G->V[i].sibl[j]);
        }
    }

    delete[](distances);
    delete[](queue);
    delete[](new_vertices);
    delete[](new_vertex_degrees);
}

void add_source_to_graph(graph G)
{
    refer v, l;

    G->n++;
    G->V[G->n-1].edgecount = G->n-1;
    G->V[G->n-1].sibl = (refer *) malloc (G->V[G->n-1].edgecount*sizeof(refer));

    for (v=0;v<G->n-1;v++)
    {
        G->V[G->n-1].sibl[v] = v;
        G->V[v].sibl[G->V[v].edgecount] = G->n-1;
        G->V[v].edgecount++;
    }

    // new number of edges
    G->m = 0;
    for (v=0;v<G->n;v++)
    {
        G->m += G->V[v].edgecount;
    }
    G->m /= 2;
    G->density = ((double)(G->m))/((double)(G->n*(G->n-1)/2));
}

void generate_graph_UDG(unsigned long n_max, unsigned long range, unsigned long grid)
{
    unsigned long i,j;
    double range_double = (double) range;
    double distance, distance_x, distance_y;

    double *points_x = new double[n_max+1];
    double *points_y = new double[n_max+1];
    refer *degrees = new refer[n_max+1];

    // alokujeme priestor pre graf
    G = (graph) malloc (sizeof(graph_data));

    // points generation
    for (i=0;i<n_max;i++)
    {
        degrees[i] = 0;
        points_x[i] = generator.random_double()*grid;
        points_y[i] = generator.random_double()*grid;
    }

    for (i=0;i<n_max;i++)
    {
        for (j=i+1;j<n_max;j++)
        {
            // if the two unit disks are within the range, we have an edge between them
            distance_x = points_x[i]-points_x[j];
            distance_y = points_y[i]-points_y[j];
            distance = distance_x*distance_x + distance_y*distance_y;
            if (distance <= range_double*range_double)
            {
                degrees[i]++;
                degrees[j]++;
            }
        }
    }

    // initialization of the base graph
    G->n = n_max;
    G->m = 0;
    for (i=0;i<n_max;i++)
    {
        G->V[i].edgecount = 0;
        G->V[i].sibl = (refer *) malloc (degrees[i]*sizeof(refer));
    }

    for (i=0;i<n_max;i++)
    {
        for (j=i+1;j<n_max;j++)
        {
            // if the two unit disks are within the range, we have an edge between them
            distance_x = points_x[i]-points_x[j];
            distance_y = points_y[i]-points_y[j];
            distance = distance_x*distance_x + distance_y*distance_y;
            if (distance <= range_double*range_double)
            {
                G->V[i].sibl[G->V[i].edgecount] = j;
                G->V[i].edgecount++;
                G->V[j].sibl[G->V[j].edgecount] = i;
                G->V[j].edgecount++;
                G->m++;
            }
        }
    }

    edges = new edgetable(2*G->m+1);
    for (i=0;i<G->n;i++)
    {
        for (j=0;j<G->V[i].edgecount;j++)
        {
            edges->insert(i,G->V[i].sibl[j]);
        }
    }

    for (i=0;i<G->n;i++)
    {
        QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
    }

    delete[](points_x);
    delete[](points_y);
    delete[](degrees);
}

void generate_graph_WS_model(refer n_max, refer k_half, double beta)
{
    refer i, j, l, k_remaining;
    refer old_vertex, new_vertex;
    bool new_vertex_free;
    int q, threshold;

    srand((unsigned long long) time(NULL));

    G = (graph) malloc (sizeof(graph_data));

    // we start only with the upper half of the connections (so that the rewiring is easier)
    G->n = n_max;
    for (i=0;i<n_max;i++)
    {
        G->V[i].sibl = (refer *) malloc (MAX_BA_DEGREE*sizeof(refer));
        G->V[i].edgecount = k_half * 2;
        G->m += k_half;
        // the upper part of the ring lattice
        j = i;
        k_remaining = k_half;
        while (k_remaining > 0)
        {
            j++;
            if (j >= G->n)
            {
                j = 0;
            }
            G->V[i].sibl[k_half-k_remaining] = j;
            k_remaining--;
        }
        j = i;
        k_remaining = k_half;
        while (k_remaining > 0)
        {
            if (0 == j)
            {
                j = G->n-1;
            }
            else
            {
                j--;
            }
            G->V[i].sibl[2*k_half-k_remaining] = j;
            k_remaining--;
        }
    }

    // rewiring: each edge between vertices (i,j), i < j is rewire it to a uniformly random vertex with probability beta
    threshold = (int) (beta * (double) RAND_MAX);
    for (i=0;i<n_max;i++)
    {
        for (j=0;j<G->V[i].edgecount;j++)
        {
            if (i < G->V[i].sibl[j])
            {
                q = generator.random(0,RAND_MAX-1);
                if (q < threshold)
                {
                    // the choice of rewiring
                    do
                    {
                        new_vertex_free = true;
                        new_vertex = generator.random(0,G->n-1);
                        if (new_vertex == i)
                        {
                            new_vertex_free = false;
                        }
                        else
                        {
                            for (l=0;l<G->V[i].edgecount;l++)
                            {
                                // if the vertex is taken (and it is different than the current), we reject it
                                if (new_vertex == G->V[i].sibl[l])
                                {
                                    new_vertex_free = false;
                                    break;
                                }
                            }
                        }
                    }
                    while (! new_vertex_free);
                    old_vertex = G->V[i].sibl[j];
                    G->V[i].sibl[j] = new_vertex;
                    // deleting the old inverse edge
                    for (l=0;l<G->V[old_vertex].edgecount;l++)
                    {
                        if (i == G->V[old_vertex].sibl[l])
                        {
                            G->V[old_vertex].sibl[l] = G->V[old_vertex].sibl[G->V[old_vertex].edgecount-1];
                            G->V[old_vertex].edgecount--;
                            break;
                        }
                    }
                    // inserting the new inverse edge
                    G->V[new_vertex].sibl[G->V[new_vertex].edgecount] = i;
                    G->V[new_vertex].edgecount++;
                }
            }
        }
    }

    edges = new edgetable(2*G->m+1);
    for (i=0;i<G->n;i++)
    {
        for (j=0;j<G->V[i].edgecount;j++)
        {
            edges->insert(i,G->V[i].sibl[j]);
        }
    }


    for (i=0;i<G->n;i++)
    {
        QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
        if (G->V[i].edgecount > MAX_BA_DEGREE)
        {
            // ToDo: overflow handling!
        }
    }
}

