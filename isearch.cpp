#include "isearch.h"

ISearch::ISearch()
{
    hweight = 1;
    breakingties = CN_SP_BT_GMAX;
}

ISearch::~ISearch(void) {}


SearchResult ISearch::startSearch(ILogger *Logger, const Map &map, const EnvironmentOptions &options)
{
    //need to implement
    Node start_vertex = map.GetStart();
    Node end_vertex = map.GetEnd();
    int64_t number_of_steps = 0, node_created = 0;
    std::vector<std::vector<double>> matrix_distance(static_cast<size_t>(map.getMapWidth()));
    std::vector<std::vector<bool>> matrix_used(static_cast<size_t>(map.getMapWidth()));
    for (size_t i = 0; i < matrix_distance.size(); i++) {
        matrix_distance[i].resize(static_cast<size_t>(map.getMapHeight()));
        matrix_used[i].resize(static_cast<size_t>(map.getMapHeight()));
    }
    for (size_t i = 0; i < matrix_distance.size(); i++) {
        for (size_t j = 0; j < matrix_distance[0].size(); j++) {
            matrix_used[i][j] = false;
        }
    }
    matrix_distance[start_vertex.i][start_vertex.j] = 0;
    start_vertex.F = &matrix_distance[start_vertex.i][start_vertex.j];
    matrix_used[start_vertex.i][start_vertex.j] = true;
    OPENED opened(start_vertex);
    while (! opened.empty()) {
        number_of_steps++;
        Node minimum = opened.extract_minimum();
        if (minimum == end_vertex)
            break;
        std::list<Node> neighbor = findSuccessors(minimum, map, options);
        for (auto it = neighbor.begin(); it != neighbor.end(); it++) {
            if (!matrix_used[(it->i)][(it->j)]) {
                node_created++;
                matrix_used[(it->i)][(it->j)] = true;
                opened.push_back(Node(it->i, it->j, &matrix_distance[(it->i)][(it->j)],
                                 computeHFromCellToCell(it->i, it->j, end_vertex.i, end_vertex.j, options)));
                matrix_distance[(it->i)][(it->j)] =
                matrix_distance[(minimum.i)][(minimum.j)] + calc_dist(*it, minimum);
            }
            if (matrix_distance[(it->i)][(it->j)] >
                matrix_distance[(minimum.i)][(minimum.j)] + calc_dist(*it, minimum))
                matrix_distance[(it->i)][(it->j)] =
                matrix_distance[(minimum.i)][(minimum.j)] + calc_dist(*it, minimum);
        }
    }
    if (matrix_used[end_vertex.i][end_vertex.j]) {
        sresult.pathfound = true;
    } else {
        sresult.pathfound = 0;
    }
    if (sresult.pathfound) {
        sresult.pathlength = matrix_distance[end_vertex.i][end_vertex.j] ;
    }
    //sresult.nodescreated =  ;
    sresult.numberofsteps = number_of_steps;
    sresult.nodescreated = node_created;
    //sresult.time = ;
    if (sresult.pathfound) {
        Node pos;
        pos.i = end_vertex.i;
        pos.j = end_vertex.j;
        while (pos != start_vertex) {
            lppath.push_front(pos);
            std::list<Node> neighbor = findSuccessors(pos, map, options);
            for (auto it = neighbor.begin(); it != neighbor.end(); it++) {
                if (matrix_used[it->i][it->j] && abs(matrix_distance[it->i][it->j] -
                        matrix_distance[pos.i][pos.j] + calc_dist(*it, pos)) < EPS) {
                    pos = *it;
                }
            }
        }
        lppath.push_front(start_vertex);
        std::pair<int, int> vctr;
        int step = 0;
        Node prev(start_vertex);
        for (auto it = lppath.begin(); it != lppath.end(); it++) {
            step++;
            if (step == 1) {
                vctr.first = prev.i - it->i;
                vctr.second = prev.j - it->j;
                prev = *it;
            } else {
                std::pair<int, int> vctr2;
                vctr2.first = prev.i - it->i;
                vctr2.second = prev.j - it->j;
                if (vctr != vctr2) {
                    hppath.push_back(prev);
                }
                vctr = vctr2;
                prev = *it;
            }
        }
        hppath.push_back(end_vertex);
        sresult.hppath = &hppath; //Here is a constant pointer
        sresult.lppath = &lppath;
    }
    return sresult;
}

std::list<Node> ISearch::findSuccessors(Node curNode, const Map &map, const EnvironmentOptions &options)
{
    std::list<Node> successors;
    bool left = 0, right = 0, up = 0, down = 0;
    if (map.CellOnGrid(curNode.i - 1, curNode.j))
        if (map.CellIsTraversable(curNode.i - 1, curNode.j)) {
            successors.push_back(Node(curNode.i - 1, curNode.j));
            left = 1;
        }
    if (map.CellOnGrid(curNode.i + 1, curNode.j))
        if (map.CellIsTraversable(curNode.i + 1, curNode.j)) {
            successors.push_back(Node(curNode.i + 1, curNode.j));
            right = 1;
        }
    if (map.CellOnGrid(curNode.i, curNode.j + 1))
        if (map.CellIsTraversable(curNode.i, curNode.j + 1)) {
            successors.push_back(Node(curNode.i, curNode.j + 1));
            down = 1;
        }
    if (map.CellOnGrid(curNode.i, curNode.j - 1))
        if (map.CellIsTraversable(curNode.i, curNode.j - 1)) {
            successors.push_back(Node(curNode.i, curNode.j - 1));
            up = 1;
        }
    if (options.allowsqueeze) {
        if (map.CellOnGrid(curNode.i - 1, curNode.j - 1))
            if (map.CellIsTraversable(curNode.i - 1, curNode.j - 1))
                successors.push_back(Node(curNode.i - 1, curNode.j - 1));
        if (map.CellOnGrid(curNode.i + 1, curNode.j + 1))
            if (map.CellIsTraversable(curNode.i + 1, curNode.j + 1))
                successors.push_back(Node(curNode.i + 1, curNode.j + 1));
        if (map.CellOnGrid(curNode.i - 1, curNode.j + 1))
            if (map.CellIsTraversable(curNode.i - 1, curNode.j + 1))
                successors.push_back(Node(curNode.i - 1, curNode.j + 1));
        if (map.CellOnGrid(curNode.i + 1, curNode.j - 1))
            if (map.CellIsTraversable(curNode.i + 1, curNode.j - 1))
                successors.push_back(Node(curNode.i + 1, curNode.j - 1));
    } else {
        if (options.allowdiagonal) {
            if (right && down) {
                if (map.CellOnGrid(curNode.i + 1, curNode.j + 1))
                    if (map.CellIsTraversable(curNode.i + 1, curNode.j + 1))
                        successors.push_back(Node(curNode.i + 1, curNode.j + 1));
            }
            if (right && up) {
                if (map.CellOnGrid(curNode.i + 1, curNode.j - 1))
                    if (map.CellIsTraversable(curNode.i + 1, curNode.j - 1))
                        successors.push_back(Node(curNode.i + 1, curNode.j - 1));
            }
            if (left && down) {
                if (map.CellOnGrid(curNode.i - 1, curNode.j + 1))
                    if (map.CellIsTraversable(curNode.i - 1, curNode.j + 1))
                        successors.push_back(Node(curNode.i - 1, curNode.j + 1));
            }
            if (left && up) {
                if (map.CellOnGrid(curNode.i - 1, curNode.j - 1))
                    if (map.CellIsTraversable(curNode.i - 1, curNode.j - 1))
                        successors.push_back(Node(curNode.i - 1, curNode.j - 1));
            }
        }
    }
    return successors;
}

double ISearch::computeHFromCellToCell(int i1, int j1, int i2, int j2, const EnvironmentOptions &options)
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

double ISearch::calc_dist(Node first_vertex, Node second_vertex)
{
    if (first_vertex.i == second_vertex.i || first_vertex.j == second_vertex.j)
        return 1.0;
    return CN_SQRT_TWO;
}

/*void ISearch::makePrimaryPath(Node curNode)
{
    //need to implement
}*/

/*void ISearch::makeSecondaryPath()
{
    //need to implement
}*/
