#ifndef GRAPHS_H
#define GRAPHS_H
#include <stdio.h>
#include "common.h"
// maximalny pocet vrcholov grafu
#define MAX_VERTICES 100050
#define MAX_LABELS_CCP 40000

/*
        ----------------
        DEKLARACIA GRAFU
        ----------------
*/

// vrchol grafu
typedef struct VERTEX
{
    refer weight;
    refer edgecount;
    refer *sibl;
} vertex;

// samotny graf
typedef struct GRAPH_DATA
{
    refer n;
    bool weighted;
    unsigned long m;
    double density;
    vertex V[MAX_VERTICES];
} graph_data;
typedef graph_data *graph;

int input_graph(FILE *source);
// dealokujeme potrebnu pamat
void free_graph();
graph get_graph();
long get_problem();
bool are_adjacent(refer v, refer w);
void generate_graph_BA_model(unsigned long w, unsigned long n_max);
void generate_graph_UDG(unsigned long n_max, unsigned long range, unsigned long grid);
void generate_shortcut_graph(graph G, refer k);

void add_source_to_graph(graph G);

#endif // GRAPHS_H
