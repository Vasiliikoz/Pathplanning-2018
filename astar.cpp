#include "astar.h"

Astar::Astar(double HW, bool BT, int par)
{
    hweight = HW;
    breakingties = BT;
    search_par = par;
}

double Astar::computeHFromCellToCell(int i1, int j1, int i2, int j2, const EnvironmentOptions &options)
{
    //need to implement
    if (options.metrictype == CN_SP_MT_DIAG) {
        return sqrt(static_cast<double>((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2))) * hweight;
    } else if (options.metrictype == CN_SP_MT_MANH) {
        return (abs(i1 - i2) + abs(j1 - j2)) * hweight;
    } else if (options.metrictype == CN_SP_MT_EUCL) {
        return sqrt(static_cast<double>((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2))) * hweight;
    } else if (options.metrictype == CN_SP_MT_CHEB) {
        if (abs(i1 - i2) > abs(j1 - j2))
            return abs(i1 - i2) * hweight;
        return abs(j1 - j2) * hweight;
    }
    return 0;
}
