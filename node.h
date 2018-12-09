#ifndef NODE_H
#define NODE_H
#include <cstdio>
#include <cinttypes>

//That's the data structure for storing a single search node.
//Although one might realize A* pathfinder relying only on g-value,
//it's a good idea to store f- and h-values explicitly for the sake of simplicity
//(and either h- or f-value is definetely needed for realizing different tie-breaking strategies).
//Backpointer is obligatory for any-angle algorithms, e.g. Theta*, and it makes sense to utilize it
//in A*-like algorithms as well for reconstructing path (after the main search phase is finished).

//So, in the end of the day, we have a sort of "universal" search-node structure
//compatable with various types of grid pathfinders (Dijkstra, A*, Jump Point Search, Theta* etc.)
//which means - that's all you need for that project.

struct Node
{
    int     i, j; //grid cell coordinates
    //double * iter_dist;//temporary
    double * F, g; //f-, g- and h-values of the search node
    double H;
    Node    *parent; //backpointer to the predecessor node (e.g. the node which g-value was used to set the g-velue of the current node)

    Node (){}
    Node (int x, int y) {
        this->i = x;
        this->j = y;
    }
    Node (int x, int y, double * dist, double h_func) {
        this->i = x;
        this->j = y;
        this->F = dist;
        this->H = h_func;
    }
    bool operator== (const Node &other) const {
        return i == other.i && j == other.j;
    }
    bool operator!= (const Node &other) const {
        return i != other.i || j != other.j;
    }
};
#endif
