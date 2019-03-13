#ifndef ISEARCH_H
#define ISEARCH_H
#include "ilogger.h"
#include "searchresult.h"
#include "environmentoptions.h"
#include <algorithm>
#include <list>
#include <vector>
#include <math.h>
#include <limits>
#include <chrono>
#include <functional>
#include <algorithm>
#include <set>
#include <type_traits>
#include <unordered_set>

class ISearch
{
    public:
        ISearch();
        ~ISearch(void);
        SearchResult startSearch(ILogger *Logger, const Map &Map, const EnvironmentOptions &options);

    protected:


        //double computeHFromCellToCell(int start_i, int start_j, int fin_i, int fin_j, const EnvironmentOptions &options) {return 0;}
        std::vector<Node> findSuccessors(Node curNode, const Map &map, const EnvironmentOptions &options);
        double calc_dist(Node first_vertex, Node second_vertex);
        double calc_distfortheta(Node first_vertex, Node second_vertex);
        double computeHFromCellToCell(int i1, int j1, int i2, int j2, const EnvironmentOptions &options);
        void make_hppath(Node start_vertex, Node end_vertex);
        bool is_parent(Node first, Node second, const Map &map);
        void ret_lp_parent(Node first, Node second, std::list<Node> &ans);
        void make_hppath_theta(Node first, Node second, std::list<Node> &ans, std::list<Node> &cur);
        //Node resetParent(Node current, Node parent, const Map &map, const EnvironmentOptions &options) {return current;}//need for Theta*
        int                             search_par;
        SearchResult                    sresult;
        std::list<Node>                 lppath, hppath;
        double                          hweight;//weight of h-value
        bool                            breakingties;//flag that sets the priority of nodes in addOpen function when their F-values is equal
        const double SQRT2 = 1.41421356237;
        const double EPS = 0.000000001;
        struct heap_vertex {
            Node vrtx;
            Node * left, right;
        };
        bool  cmp(Node& a, Node& b){
            if  (a.i == b.i && a.j == b.j)
                return true;
            return false;
        }
        struct cmp_str{
            double F, H;
            bool br;
            cmp_str() {}
            cmp_str(double a, double b, bool brk) {
                F = a;
                H = b;
                br = brk;
            }
            bool operator< (cmp_str other) {
                if (abs(this->F + this->H - (other.F + other.H)) < 0.000000001) {
                    if (this->br)
                        return this->F < other.F;
                    return this->H < other.H;
                } else {
                    return (this->F + this->H < other.F + other.H);
                }
            }
        };

        struct treap {
            Node val;
            size_t prior;
            treap * l, * r;
            cmp_str fun;
            treap() { }
            treap (Node value, size_t p){
                val = value;
                prior = p;
                l = nullptr;
                r = nullptr;
                cmp_str fun(value.F,value.H, value.brk);
            }
        };
        typedef treap * iter;

        void new_fun(iter & it){
            if (it == nullptr)
                return;
            cmp_str mn(it->val.F,it->val.H, it->val.brk);
            if (it->l != nullptr)
                if (it->l->fun < mn)
                    mn = it->l->fun;
            if (it->r != nullptr)
                if (it->r->fun < mn)
                    mn = it->r->fun;
            it->fun = mn;
        }

        void split(iter it, Node val, iter & l, iter & r) {
            if (it == nullptr) {
                l = nullptr;
                r = l;
            } else if (val < it->val) {
                split (it->l, val, l, it->l);
                r = it;
                new_fun(r);
            } else {
                split (it->r, val, it->r, r);
                l = it;
                new_fun(l);
            }
        }

        void insert(iter & t, iter it) {
            if (t == nullptr) {
                t = it;
                //std::cout << t->val.i << " " << t->val.j << "\n";
            } else if (it->prior > t->prior) {
                split (t, it->val, it->l, it->r);
                t = it;
            } else {
                if (it->val < t->val) {
                    insert (t->l, it);
                } else {
                    insert (t->r, it);
                }
            }
            new_fun(t);
        }

        void merge(iter & it, iter l, iter r) {
            if (!l || !r) {
                if (l == nullptr) {
                    it = r;
                } else {
                    it = l;
                }
            } else {
                if (l->prior > r->prior) {
                    merge (l->r, l->r, r);
                    it = l;
                    //delete l;
                } else {
                    merge (r->l, l, r->l);
                    it = r;
                    //delete r;
                }
            }
            new_fun(it);
        }

        void erase(iter & it, Node val) {
            if (it == nullptr) {
                std::cout << "!nl\n";
                return;
            }
            if (cmp(it->val, val)) {
                merge (it, it->l, it->r);
            } else {
                if (it->val < val) {
                    erase (it->r, val);
                } else {
                    erase (it->l, val);
                }
            }
            new_fun(it);
        }
        bool cmp1(cmp_str &d1, cmp_str &d2){
            return (abs(d1.F - d2.F) + abs(d1.H - d2.H)< EPS);
        }
        Node find_min(iter it) {
            cmp_str d = it->fun;
            cmp_str n_w(it->val.F,it->val.H, it->val.brk);
            if (cmp1(it->fun, n_w))
                return it->val;
            if (it->l != nullptr) {
                if (cmp1(it->l->fun,d))
                    return find_min(it->l);
            }
            if (it->r != nullptr)
                if (cmp1(it->r->fun,d))
                    return find_min(it->r);
            std::cout << "\n";
            std::cout << "find_min\n";
            std::cout << "\n";
            return Node();
        }

        bool find(iter it, Node a) {
            if (it == nullptr)
                return false;
            if (cmp(it->val, a))
                return true;
            if (it->val < a)
                return find(it->r, a);
            return find(it->l, a);
        }

        Node find1(iter it, Node a) {
            if (cmp(it->val, a))
                return it->val;
            if (it->val < a)
                return find1(it->r, a);
            return find1(it->l, a);
        }

        bool new_value(iter & it,Node a){
            if (it == nullptr)
                return false;
            if (cmp(it->val, a)) {
                if (it->val.F > a.F) {
                    it->val.F = a.F;
                    new_fun(it);
                }
                return true;
            }
            if (it->val < a) {
                bool b = new_value(it->r, a);
                new_fun(it);
                return b;
            }
            bool b = new_value(it->l, a);
            new_fun(it);
            return b;
        }

        void outp(iter it, std::vector<Node> &v){
            if (it != nullptr) {
                v.push_back(it->val);
                outp(it->l, v);
                outp(it->r, v);
            }
        }

        class OPENED{
        private:
            iter head;
        public:
            OPENED(){
                head = nullptr;
            }
            bool notempty(){
                return head != nullptr;
            }
            void insert(Node a){
                iter it = new treap(a, std::rand());
                ISearch p;
                p.insert(head, it);
                //std::cout << head->val.i << " " << head->val.j << "\n";
            }
            void erase(Node a){
                ISearch p;
                p.erase(head, a);
            }
            Node erase_minimum(){
                ISearch p;
                Node a = p.find_min(head);
                p.erase(head, a);
                return a;
            }
            void minm(){
                ISearch p;
                Node a = p.find_min(head);
                std::cout<< a.i << " " << a.j << "\n";
            }
            int new_value(Node a){
                ISearch p;
                bool b = p.new_value(head, a);
                if (!b) {
                    this->insert(a);
                    return 1;
                }
                return 0;
            }
            bool find(Node a){
                ISearch p;
                return p.find(head, a);
            }
            Node find1(Node a){
                ISearch p;
                return p.find1(head, a);
            }
            void output(std::vector<Node> &v){
                ISearch p;
                p.outp(head, v);
            }
        };

        class CLOSED{
        private:
            std::vector<OPENED> vctr;
        public:
            CLOSED(int a) {
                vctr.resize(a);
            }
            void insert(Node a){
                vctr[a.i].insert(a);
            }
            bool find(Node a){
                return vctr[a.i].find(a);
            }
            Node find1(Node a) {
                return vctr[a.i].find1(a);
            }
            void output(std::vector<Node> &v) {
                for (int i = 0; i < vctr.size(); i++) {
                    vctr[i].output(v);
                }
            }
        };

};
#endif
