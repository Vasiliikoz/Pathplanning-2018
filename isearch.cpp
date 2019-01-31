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
    auto start = std::chrono::system_clock::now();
    std::srand(std::time(nullptr));
    Node start_vertex = map.GetStart();
    Node end_vertex = map.GetEnd();
    OPENED opened;
    CLOSED closed;
    int64_t number_of_steps = 0, node_created = 0;
    opened.insert(start_vertex);
    while (opened.notempty()) {
        Node a = opened.erase_minimum();
        closed.insert(a);
        if (a.i == end_vertex.i && a.j == end_vertex.j) {
            sresult.pathlength = a.F;
            break;
        }
        number_of_steps++;
        auto lst = findSuccessors(a, map, options);
        for (auto i = lst.begin(); i != lst.end(); i++) {
            Node b(i->i, i->j, calc_dist(*i, a) + a.F, computeHFromCellToCell(i->i, i->j, end_vertex.i, end_vertex.j, options), a.i, a.j);
            if (!closed.find(b)){
                opened.new_value(b);
            } else {
                node_created++;
            }
        }
    }
    sresult.pathfound = (closed.find(end_vertex));
    std::cout<< sresult.pathlength << "\n";
    //sresult.nodescreated =  ;
    sresult.numberofsteps = number_of_steps;
    sresult.nodescreated = node_created;
    auto end = std::chrono::system_clock::now();
    sresult.time = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1000000000;
    //sresult.time = ;
    if (sresult.pathfound) {
        Node pos;
        pos.i = end_vertex.i;
        pos.j = end_vertex.j;
        while (pos != start_vertex) {
            lppath.push_front(pos);
            Node nxt = closed.find1(pos);
            pos.i = nxt.previous_i;
            pos.j = nxt.previous_j;
        }
        lppath.push_front(start_vertex);
        auto end = std::chrono::system_clock::now();
        sresult.time = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1000000000;
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
