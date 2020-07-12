#ifndef EDGETABLE_H
#define EDGETABLE_H
#include <stdio.h>
#include <stdlib.h>
#include "graphs.h"

class hashset
{
public:
    bool occ;
    refer v;
    refer w;
    hashset *next;
};

class edgetable
{
private:
    unsigned long long max_edges;
    hashset *H;
public:
    edgetable(unsigned long long max_edges);
    ~edgetable();
    void set_max_edges(unsigned long long max_edges);
    unsigned long long hash(refer v, refer w);
    bool isin(refer v, refer w);
    void insert(refer v, refer w);
    void clear();
};

#endif // EDGETABLE_H
