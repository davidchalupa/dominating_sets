#ifndef CLI_H
#define CLI_H
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <QTime>
#include "graphs.h"
#include "random_generator.h"
#include "time.h"
#include "statistics.h"
#include "problem.h"

#define PATH_MAX 260

class cli
{
private:
    refer i;
    long succ_rate;
    unsigned long max_deg,min_deg;
    unsigned long long current_time,avg_time,avg_iter;
    double avg_deg,stdev_deg;
    refer k;
    char filename[PATH_MAX];
    double avg_brelaz;
    double avg_iggcc;
    double avg_gis;
    double avg_rls_is;
    refer lower_bound;
    refer min_iggcc;
    unsigned long w,n_max,range,grid;
    bool is_optimal;
    QTime timer;    
    FILE *source;
    graph G;
    refer max_label;
    unsigned long long scanned = 0;
    // overall statistics
    double greedy_max_overall;
    double greedy_min_overall;
    double greedy_avg_overall;
    double undominated_vertices_overall_rls;
    double dominating_set_size_overall_rls;
    double undominated_vertices_overall_ma;
    double dominating_set_size_overall_ma;
    double undominated_vertices_overall_msob;
    double dominating_set_size_overall_msob;
    double undominated_vertices_overall_ob;
    double dominating_set_size_overall_ob;
    double undominated_vertices_overall_obma;
    double dominating_set_size_overall_obma;
    double undominated_vertices_overall_acols;
    double dominating_set_size_overall_acols;
    double undominated_vertices_overall_acoppls;
    double dominating_set_size_overall_acoppls;
    double undominated_vertices_overall_acols_s;
    double dominating_set_size_overall_acols_s;
    long long t_obdom;
    long long t_msobdom;
    long long t_acols;
    long long t_acoppls;
    long long t_acols_s;
    int choose_instance();
    int generate_instance();
    int generate_instance_udg();
    int choose_algorithm();
    int compute_statistics();
    void sleep(int milisec);
    bool result[MAX_VERTICES];
    bool best_result[MAX_VERTICES];
    refer histogram[MAX_VERTICES];
    refer sizes[MAX_VERTICES];
    long bad_scenario_counter;
    refer degree_distrib[MAX_VERTICES];
    bool current_solution_exhaustive[MAX_VERTICES];
    refer current_dominating_set_vertices_exhaustive[MAX_VERTICES];
    refer current_vertex_stack[MAX_VERTICES];
    refer current_vertex_stack_size;
    refer depth_stack[MAX_VERTICES];
    refer depth_stack_size;
    refer next_vertex_stack[MAX_VERTICES];
    refer next_vertex_stack_size;
    refer component_labels[MAX_VERTICES];
    problem *problemdom_exhaustive;
    bool found_optimum_exhaustive;
    double approx;
    long discrep;
    refer components;
    refer current_component_label_exhaustive;
    refer current_component_dominating_set_size;
    bool currently_dominated[MAX_VERTICES];
    //bool modified[MAX_VERTICES];
    refer current_nondominated_vertices;
    void domset_to_mps(bool relaxation);
    void exhaustive_search_core(refer current_vertex, refer depth, refer current_size);
    void exhaustive_search(refer current_size);
    void component_labeling();
    int process_from_file();
    int process_ba();
    int process_udg();
public:
    cli();
    int start_cli();
};

#endif // CLI_H
