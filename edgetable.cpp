#include "edgetable.h"

edgetable::edgetable(unsigned long long max_edges)
{
    unsigned long long i;
    this->max_edges = max_edges;
    H = new hashset[max_edges];
    for (i=0;i<max_edges;i++)
    {
        H[i].occ = false;
        H[i].next = NULL;
    }
}

edgetable::~edgetable()
{
    clear();
    delete[](H);
}

void edgetable::set_max_edges(unsigned long long max_edges)
{
    unsigned long long i;
    clear();
    delete(H);
    this->max_edges = max_edges;
    H = new hashset[max_edges];
    for (i=0;i<max_edges;i++)
    {
        H[i].occ = false;
        H[i].next = NULL;
    }
}

// hashovacia funkcia
unsigned long long edgetable::hash(refer v, refer w)
{
    return (((unsigned long long)(v)*MAX_VERTICES+(unsigned long long)w)*2654435761) % (unsigned long long) max_edges;
}

// funkcia: je val elementom mnoziny?
bool edgetable::isin(refer v, refer w)
{
    unsigned long long hsh;
    hsh = hash(v,w);
    // ak este nemame obsadene miesto
    if (! H[hsh].occ)
    {
        return false;
    }
    // ak sa nieco na mieste nachadza
    else
    {
        bool ok = false;
        hashset *actelm;
        actelm = NULL;
        actelm = &H[hsh];
        while (actelm != NULL)
        {
            if (actelm->v == v && actelm->w == w)
            {
                ok = true;
                break;
            }
            actelm = actelm->next;
        }
        return ok;
    }
}

// funkcia: vlozit val do mnoziny
void edgetable::insert(refer v, refer w)
{
    unsigned long long hsh;
    hsh = hash(v,w);
    // ak este nemame obsadene miesto
    if (! H[hsh].occ)
    {
        H[hsh].v = v;
        H[hsh].w = w;
        H[hsh].next = NULL;
        H[hsh].occ = true;
    }
    // ak uz nastala kolizia
    else
    {
        hashset *actelm,*newelm;
        actelm = NULL;
        // vytvorime novy element
        newelm = new hashset();
        newelm->v = v;
        newelm->w = w;
        newelm->next = NULL;
        // vlozime element na koniec zoznamu
        actelm = &H[hsh];
        while (actelm->next != NULL)
        {
            actelm = actelm->next;
        }
        actelm->next = newelm;
    }
}

void edgetable::clear()
{
    unsigned long long i;
    for (i=0;i<max_edges;i++)
    {
        if (H[i].occ && H[i].next != NULL)
        {
            hashset *actelm;
            hashset *delelm;
            actelm = H[i].next;
            while (actelm != NULL)
            {
                delelm = actelm;
                actelm = actelm->next;
                delete(delelm);
            }
            H[i].occ = false;
            H[i].next = NULL;
        }
    }
}
