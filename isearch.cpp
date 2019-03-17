#include "isearch.h"


ISearch::ISearch()
{
    hweight = 1;
    breakingties = CN_SP_BT_GMAX;
}

ISearch::~ISearch(void) {}


SearchResult ISearch::startSearch(ILogger *Logger, const Map &map, const EnvironmentOptions &options)
{
    auto start = std::chrono::system_clock::now();
    std::srand(std::time(nullptr));
    Node start_vertex = map.GetStart();
    start_vertex.brk = breakingties;
    Node end_vertex = map.GetEnd();
    end_vertex.brk = breakingties;
    OPENED opened;
    CLOSED closed((map.getMapWidth() > map.getMapHeight())? map.getMapWidth():map.getMapHeight());
    unsigned int number_of_steps = 0, node_created = 0;
    opened.insert(start_vertex);
    sresult.pathlength = 0;
    if (search_par == 0) {
    while (opened.notempty()) {
        if (Logger->loglevel == CN_LP_LEVEL_FULL_WORD) {
            std::vector<Node> op;
            op.reserve(node_created + 1);
            opened.output(op);
            std::vector<Node> cl;
            cl.reserve(number_of_steps + 1);
            closed.output(cl);
            Logger->writeToLogOpenClose(op, cl, number_of_steps);
        } 
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
            b.brk = breakingties;
            if (!closed.find(b)){
                node_created += opened.new_value(b);
            }
        }
    }
    } else if (search_par == 4) {

        while (opened.notempty()) {
            if (Logger->loglevel == CN_LP_LEVEL_FULL_WORD) {
                std::vector<Node> op;
                op.reserve(node_created + 1);
                opened.output(op);
                std::vector<Node> cl;
                cl.reserve(number_of_steps + 1);
                closed.output(cl);
                Logger->writeToLogOpenClose(op, cl, number_of_steps);
            }
            Node a = opened.erase_minimum();
            closed.insert(a);
            if (a.i == end_vertex.i && a.j == end_vertex.j) {
                sresult.pathlength = a.F;
                //std::cout << sresult.pathlength << "\n";
                break;
            }
            number_of_steps++;
            auto lst = findSuccessors(a, map, options);
            for (auto i = lst.begin(); i != lst.end(); i++) {
                Node b_n(i->i, i->j);
                double d = HF_Theta(i->i, i->j, end_vertex.i, end_vertex.j);
                if (closed.find(b_n))
                    continue;
                Node parent_a(a.previous_i, a.previous_j);
                if (is_parent(parent_a, *i, map)) {
                    parent_a = closed.find1(parent_a);
                    Node b(i->i, i->j, calc_distfortheta(*i, parent_a) + parent_a.F,
                           d, parent_a.i, parent_a.j);
                    b.brk = breakingties;
                    node_created += opened.new_value(b);
                } else {
                    Node b(i->i, i->j, calc_dist(*i, a) + a.F,
                           d, a.i, a.j);
                    b.brk = breakingties;
                    node_created += opened.new_value(b);
                }
            }
        }

    }
    if (Logger->loglevel == CN_LP_LEVEL_MEDIUM_WORD) {
        std::vector<Node> op;
        op.reserve(node_created + 1);
        opened.output(op);
        std::vector<Node> cl;
        cl.reserve(number_of_steps + 1);
        closed.output(cl);
        Logger->writeToLogOpenClose(op, cl, number_of_steps);
    }
    sresult.pathfound = (closed.find(end_vertex));
    sresult.numberofsteps = number_of_steps;
    sresult.nodescreated = node_created;
    auto end = std::chrono::system_clock::now();
    sresult.time = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1000000000;
    if (sresult.pathfound) {
        Node pos;
        pos.i = end_vertex.i;
        pos.j = end_vertex.j;
        //pos.g = sresult.pathlength;
        while (pos != start_vertex) {
            lppath.push_front(pos);
            Node nxt = closed.find1(pos);
            pos.i = nxt.previous_i;
            pos.j = nxt.previous_j;
            //pos.g = nxt.F;
        }
        lppath.push_front(start_vertex);
        auto end = std::chrono::system_clock::now();
        sresult.time = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()) / 1000000000;
        if (search_par == 0) {
            make_hppath(start_vertex, end_vertex);
        } else {
            make_hppath_theta(start_vertex, end_vertex, hppath, lppath);
            swap(lppath, hppath);
        }
        sresult.hppath = &hppath;
        sresult.lppath = &lppath;
    }
    return sresult;
}

std::vector<Node> ISearch::findSuccessors(Node curNode, const Map &map, const EnvironmentOptions &options)
{
    std::vector<Node> successors;
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

double ISearch::HF_Theta(int i1, int j1, int i2, int j2)
{
    return sqrt(static_cast<double>((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2))) * hweight;
}

double ISearch::computeHFromCellToCell(int i1, int j1, int i2, int j2, const EnvironmentOptions &options)
{
    if (options.metrictype == CN_SP_MT_DIAG) {
        return (abs(abs(i1 - i2) - abs(j1 - j2)) + SQRT2 * (std::min(abs(i1 - i2),abs(j1 - j2)))) * hweight;
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

void ISearch::make_hppath(Node start_vertex, Node end_vertex)
{
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
}


void ISearch::make_hppath_theta(Node first, Node second, std::list<Node> &ans, std::list<Node> &cur) {
    Node st = *cur.begin();
    for (auto i = cur.begin(); i != cur.end(); i++) {
        if (((st.i + i->j == st.j + i->i) && (st.i - i->i == 1)) || ((st.i == i->i || st.j == i->j) &&
                abs(st.i - i->i + st.j - i->j) <= 1)) {
            st = *i;
        } else {
            std::list<Node> n;
            ret_lp_parent(st, *i, n);
            for (auto it = n.begin(); it != n.end(); it++)
                ans.push_back(*it);
            st = *i;
        }
        ans.push_back(*i);
    }
}

void ISearch::ret_lp_parent(Node first, Node second, std::list<Node> &ans)
{
    if (first.j == second.j) {
        int min = std::min(first.i, second.i);
        int max = std::max(first.i, second.i);
        for (int i = min; i <= max; i++)
            ans.push_back(Node(i, first.j));
        if (first.i > second.i)
            ans.reverse();
        return;
    }
    if (first.i == second.i) {
        int min = std::min(first.j, second.j);
        int max = std::max(first.j, second.j);
        for (int j = min; j <= max; j++)
            ans.push_back(Node(first.i, j));
        if (first.j > second.j)
            ans.reverse();
        return;
    }
    bool b = false;
    if (first.j > second.j) {
        std::swap(first, second);
        b = true;
    }
    double k = double(second.i - first.i) / double(second.j - first.j);
    Node p_n = first;
    if (first.i < second.i) {
        double b = first.i - k * first.j + k / 2 - 0.5;
        auto func = [k, b](int x) {return k * x + b;};
        while (true) {
            if (p_n.i > func(p_n.j)) {
                p_n.j++;
            } else {
                p_n.i++;
            }
            if (p_n.i == second.i && p_n.j == second.j)
                break;
            ans.push_back(Node(p_n.i, p_n.j));
        }
    } else {
        double b = first.i - k * first.j + k / 2 + 0.5;
        auto func = [k, b](int x) {return k * x + b;};
        while (true) {
            if (p_n.i < func(p_n.j)) {
                p_n.j++;
            } else {
                p_n.i--;
            }
            if (p_n.i == second.i && p_n.j == second.j)
                break;
            ans.push_back(Node(p_n.i, p_n.j));
        }
    }
    if (b) {
        ans.reverse();
    }
    return;
}

double ISearch::calc_dist(Node first_vertex, Node second_vertex)
{
    if (first_vertex.i == second_vertex.i || first_vertex.j == second_vertex.j)
        return 1.0;
    return CN_SQRT_TWO;
}

double ISearch::calc_distfortheta(Node first_vertex, Node second_vertex)
{
    return sqrt((first_vertex.i - second_vertex.i) * (first_vertex.i - second_vertex.i) +
                (first_vertex.j - second_vertex.j) * (first_vertex.j - second_vertex.j));
}


bool ISearch::is_parent(Node first, Node second, const Map &map)
{
    if (first.j == second.j) {
            int min = std::min(first.i, second.i);
            int max = std::max(first.i, second.i);
            for (int i = min; i <= max; i++)
                if (map.CellIsObstacle(i, first.j))
                    return false;
            return true;
        }
        if (first.i == second.i) {
            int min = std::min(first.j, second.j);
            int max = std::max(first.j, second.j);
            for (int j = min; j <= max; j++)
                if (map.CellIsObstacle(first.i, j))
                    return false;
            return true;
        }
    if (first.j > second.j) {
        std::swap(first, second);
    }
    double k = double(second.i - first.i) / double(second.j - first.j);
    if (k == 2.0) {
        int i = std::min(first.i, second.i);
        for (int j = first.j; j <= second.j; j++)
            if (map.CellIsObstacle(i++, j))
                return false;
        return true;
    }
    Node p_n = first;
    if (first.i < second.i) {
        double b = first.i - k * first.j + k / 2 - 0.5;
        auto func = [k, b](int x) {return k * x + b;};
        while (p_n.i != second.i && p_n.j != second.j) {
            if (p_n.i > func(p_n.j)) {
                p_n.j++;
            } else {
                p_n.i++;
            }
            if (map.CellIsObstacle(p_n.i, p_n.j))
                return false;
        }
    } else {
        double b = first.i - k * first.j + k / 2 + 0.5;
        auto func = [k, b](int x) {return k * x + b;};
        while (p_n.i != second.i && p_n.j != second.j) {
            if (p_n.i < func(p_n.j)) {
                p_n.j++;
            } else {
                p_n.i--;
            }
            if (map.CellIsObstacle(p_n.i, p_n.j))
                return false;
        }
    }
    return true;
}
