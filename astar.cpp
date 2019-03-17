#include "astar.h"

Astar::Astar(double HW, bool BT, int par)
{
    hweight = HW;
    breakingties = BT;
    search_par = par;
}
