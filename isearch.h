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
        //CODE HERE
        //Try to split class functionality to the methods that can be re-used in successors classes,
        //e.g. classes that realize A*, JPS, Theta* algorithms

        //Hint 1. You definetely need class variables for OPEN and CLOSE

        //Hint 2. It's a good idea to define a heuristic calculation function, that will simply return 0
        //for non-heuristic search methods like Dijkstra

        //Hint 3. It's a good idea to define function that given a node (and other stuff needed)
        //will return it's sucessors, e.g. unordered list of nodes

        //Hint 4. It's a good idea to define a resetParent function that will be extensively used for Theta*
        //and for A*/Dijkstra/JPS will exhibit "default" behaviour

        //Hint 5. The last but not the least: working with OPEN and CLOSE is the core
        //so think of the data structures that needed to be used, about the wrap-up classes (if needed)
        //Start with very simple (and ineffective) structures like list or vector and make it work first
        //and only then begin enhancement!



        //double computeHFromCellToCell(int start_i, int start_j, int fin_i, int fin_j, const EnvironmentOptions &options) {return 0;}
        std::list<Node> findSuccessors(Node curNode, const Map &map, const EnvironmentOptions &options);
        double calc_dist(Node first_vertex, Node second_vertex);
        double computeHFromCellToCell(int i1, int j1, int i2, int j2, const EnvironmentOptions &options);
        //void makePrimaryPath(Node curNode);//Makes path using back pointers
        //void makeSecondaryPath();//Makes another type of path(sections or points)
        //Node resetParent(Node current, Node parent, const Map &map, const EnvironmentOptions &options) {return current;}//need for Theta*

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
        struct treap {
            Node val;
            size_t prior;
            treap * l, * r;
            double fun;
            treap() { }
            treap (Node value, size_t p){
                val = value;
                prior = p;
                l = nullptr;
                r = nullptr;
                fun = value.F + value.H;
            }
        };
        typedef treap * iter;

        void new_fun(iter & it){
            if (it == nullptr)
                return;
            double mn = it->val.F + it->val.H;
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
        bool cmp1(double d1, double d2){
            return (abs(d1 - d2) < EPS);
        }
        Node find_min(iter it) {
            double d = it->fun;
            if (cmp1(it->fun,it->val.F + it->val.H))
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

        void outp(iter it){
            if (it == nullptr) {
                std::cout << "null\n";
            } else {
                std::cout << it->val.i << " " << it->val.j << "\n";
                std::cout << "left\n";
                outp(it->l);
                std::cout << "right\n";
                outp(it->r);
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
            void outp(){
                ISearch p;
                p.outp(this->head);
            }
            void mm(){
                if (head != nullptr)
                    std::cout << head->val.i << " " << head->val.j << " " << head->fun << "\n";
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
        };

        class CLOSED{
        private:
            iter head;
        public:
            CLOSED(){
                head = nullptr;
            }
            void insert(Node a){
                iter it = new treap(a, std::rand());
                ISearch p;
                p.insert(head, it);
            }
            void outp(){
                ISearch p;
                p.outp(this->head);
            }
            void mm(){
                if (head != nullptr)
                    std::cout << head->val.i << " " << head->val.j << " " << head->fun << "\n";
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
            void new_value(Node a){
                ISearch p;
                bool b = p.new_value(head, a);
                if (!b)
                    this->insert(a);
            }
            bool find(Node a){
                ISearch p;
                return p.find(head, a);
            }
            Node find1(Node a){
                ISearch p;
                return p.find1(head, a);
            }
        };

};
#endif
