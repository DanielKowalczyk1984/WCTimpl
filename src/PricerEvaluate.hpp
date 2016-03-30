#include "tdzdd/DdEval.hpp"
#include "tdzdd/dd/NodeTable.hpp"
#include <cfloat>
#include <climits>
#include <vector>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>

/**
 * ZDD
 */
template<typename T>
class PricerInfoZDD
{
    public:
        T  *obj;
        int *cost;
        boost::dynamic_bitset<>  *A;

        PricerInfoZDD(): obj(0), cost(0), A(0)
        {
        }

        ~PricerInfoZDD()
        {
            if (obj != 0) {
                delete[] obj;
                obj = 0;
            }

            if (cost != 0) {
                delete[] cost;
                cost = 0;
            }

            if (A != 0) {
                delete[] A;
                A = 0;
            }
        }

        size_t get_max(const int &L)
        {
            T max = obj[0];
            size_t it_max = 0;

            for (size_t i = 0; i < L; ++i) {
                if (max < obj[i]) {
                    max = obj[i];
                    it_max = i;
                }
            }

            return it_max;
        }

        PricerInfoZDD &operator=(const PricerInfoZDD &other)
        {
            obj = other.obj;
            cost = other.cost;
            A = other.A;
            return *this;
        };
};

/**
 * BDD
 */
template<typename T>
class PricerInfoBDD
{
    public:
        T obj;
        std::vector<int> jobs;
        int sum_w;
        int sum_p;
        int cost;

        PricerInfoBDD &operator=(const PricerInfoBDD &other)
        {
            sum_w = other.sum_w;
            sum_p = other.sum_p;
            cost = other.cost;
            obj = other.obj;
            jobs.clear();
            jobs = other.jobs;
            return *this;
        };

        friend std::ostream &operator<<(std::ostream &os, PricerInfoBDD<T> const &o)
        {
            os << "max = " << o.obj << "," << std::endl << "cost = " << o.cost << std::endl;
            return os;
        }
};

template<typename T>
class PricerWeightBDD
{
    public:
        T obj;
        bool take;
        int sum_w;
        int sum_p;
        int cost;

        PricerWeightBDD(): obj(0), take(false), sum_w(0), sum_p(0), cost(0)
        {
        }

        PricerWeightBDD &operator=(const PricerWeightBDD &other)
        {
            sum_w = other.sum_w;
            sum_p = other.sum_p;
            cost = other.cost;
            obj = other.obj;
            take = other.take;
            return *this;
        };

        friend std::ostream &operator<<(std::ostream &os, PricerWeightBDD<T> const &o)
        {
            os << "max = " << o.obj << "," << std::endl << "cost = " << o.cost << std::endl;
            return os;
        }

        void init_terminal_node(int one)
        {
            obj = one ? 0.0 : -1871286761.0;
            sum_p = 0;
            take = false;
        }

        void init_node(int weight)
        {
            obj = 0.0;
            sum_p = weight;
            take = false;
        }


};

template<typename T>
class PricerWeightZDD
{
    public:
        std::unordered_map<int, T> obj;
        std::unordered_map<int, bool> take;
        std::vector<int> weight;

        PricerWeightZDD()
        {
        };

        ~PricerWeightZDD()
        {
        }

        /*void alloc_node(const std::vector<int> *v)
        {
            int max = *std::max_element(v->begin(), v->end());
            obj = new T[max + 1];
            take = new bool[max + 1];
            weight = v;

            for (auto &i : *v) {
                take[i] = false;
            }
        }*/

        /*void alloc_terminal_node(int H_min, int H_max, int H, int one)
        {
            int end = one ? H_max + 1 : H;
            obj = new T[end];
            take = new bool[end];

            for (int i = one ? H_min : 0; i < H_max + 1 ; i++) {
                obj[i] = one ? 0 : -1871286761.0;
            }
        }*/

        void add_weight(int _weight){
            if(obj.find(_weight) == obj.end()) {
                weight.push_back(_weight);
                obj[_weight] = 0.0;
                take[_weight] = false;
            }
        }

        void init_terminal_node(int j, int H_max){
            int end = j ? H_max + 1 : 2*H_max;
            for(size_t i = 0; i < end;i++){
                obj[i] = j ? 0.0 : -1871286761.0;
                take[i] = false;
            }
        }

};

template<typename T>
class PricerFarkasZDD
{
    public:
        T obj;
        bool take;

        PricerFarkasZDD(): obj(0), take(0)
        {
        };

        ~PricerFarkasZDD()
        {
        };

        void init_terminal_node(int one)
        {
            obj = one ? 0.0 : -1871286761.0;
        }

        void init_node()
        {
            take = false;
        }
};


template<typename T>
class Optimal_Solution
{
    public:
        T obj;
        int cost;
        int C_max;
        std::vector<int> jobs;
        Optimal_Solution(PricerInfoZDD<T> *node, const int L)
        {
            int it_max = node->get_max(L);
            boost::dynamic_bitset<> *ptr_max = &node->A[it_max];
            size_t it = ptr_max->find_first();
            obj = node->obj[it_max];
            cost = node->cost[it_max];

            while (it != boost::dynamic_bitset<>::npos) {
                if ((*ptr_max)[it]) {
                    jobs.push_back(it);
                }

                it = ptr_max->find_next(it);
            }
        }

        Optimal_Solution(PricerInfoBDD<T> *node)
        {
            obj = node->obj;
            jobs = node->jobs;
            cost = node->cost;
        }

        Optimal_Solution()
        {
            obj = 0;
            cost = 0;
            C_max = 0;
        }

        Optimal_Solution &operator=(const Optimal_Solution &other)
        {
            obj = other.obj;
            cost = other.cost;
            jobs = other.jobs;
            C_max = other.C_max;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, Optimal_Solution<T> const &o)
        {
            os << "obj = " << o.obj << "," << std::endl << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;
            std::vector<int>::const_iterator it = o.jobs.begin();

            for (; it != o.jobs.end(); ++it) {
                std::cout << *it << " ";
            }

            std::cout << std::endl;
            return os;
        }
};

template<typename E, typename T>
class DurationZDD: public tdzdd::DdEval<E, PricerInfoZDD<T>, Optimal_Solution<T> >
{
        T *pi;
        int *p;
        int *w;
        int nbjobs;
        int L;

    public:
        DurationZDD(T *_pi, int *_p, int *_w, int _nbjobs, int _L)
            : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs), L(_L)
        {
        };

        void evalTerminal(PricerInfoZDD<T> &n, bool one)
        {
            for (unsigned i = 0; i < L + 1; ++i) {
                n.obj[i] = one ? 0 : -15656873.0;
            }
        }

        void evalNode(PricerInfoZDD<T>  &n, int i, tdzdd::DdValues<PricerInfoZDD<T>, 2>    &values) const
        {
            int j = nbjobs - i;
            assert(j >= 0 && j <= nbjobs - 1);
            PricerInfoZDD<T> *n0 = values.get_ptr(0);
            PricerInfoZDD<T> *n1 = values.get_ptr(1);

            for (int k = 0; k < L + 1; ++k) {
                int it = k + p[j];

                if (it <= L) {
                    if (n1->obj[it] < n.obj[k] - w[j]*it + pi[j]) {
                        n1->obj[it] = n.obj[k] - w[j] * it + pi[j];
                        n1->cost[it] = n.cost[k] + w[j] * it;
                        n1->A[it].clear();
                        n1->A[it] = n.A[k];
                        n1->A[it][j] = true;
                    }
                }

                if (n0->obj[k] < n.obj[k]) {
                    n0->obj[k] = n.obj[k];
                    n0->cost[k] = n.cost[k];
                    n0->A[k] = n.A[k];
                }
            }
        }

        void initializenode(PricerInfoZDD<T> &n)
        {
            n.obj = new T [L + 1];
            n.cost = new int [L + 1];
            n.A = new boost::dynamic_bitset<> [L + 1];

            for (unsigned i = 0; i < L + 1; ++i) {
                n.obj[i] = -10983290.0;
                n.A[i].resize(nbjobs);
                n.cost[i] = 0;
            }
        }

        void initializerootnode(PricerInfoZDD<T> &n)
        {
            n.obj = new T [L + 1];
            n.cost = new int [L + 1];
            n.A = new boost::dynamic_bitset<> [L + 1];

            for (unsigned i = 0; i < L + 1; ++i) {
                n.obj[i] = pi[nbjobs];
                n.A[i].resize(nbjobs);
                n.cost[i] = 0;
            }
        }

        Optimal_Solution<T> get_objective(PricerInfoZDD<T> *n)
        {
            Optimal_Solution<T> sol(n, L);
            return sol;
        }
};

template<typename E, typename T>
class DurationBDD: public tdzdd::DdEval<E, PricerInfoBDD<T>, Optimal_Solution<T> >
{
        T *pi;
        int *p;
        int *w;
        int nbjobs;

    public:
        DurationBDD(T *_pi, int *_p, int *_w, int _nbjobs)
            : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs)
        {
        };

        void evalTerminal(PricerInfoBDD<T> &n, bool one)
        {
            n.obj = one ? 0 : -1871286761.0;
            n.cost = 0;
            n.sum_p = 0;
            n.jobs.resize(0);
        }

        void evalNode(PricerInfoBDD<T> &n, int i, tdzdd::DdValues<PricerInfoBDD<T>, 2>    &values) const
        {
            int j = nbjobs - i;
            assert(j >= 0 && j <= nbjobs - 1);
            PricerInfoBDD<T> *n0 = values.get_ptr(0);
            PricerInfoBDD<T> *n1 = values.get_ptr(1);

            if (n0->obj < n.obj) {
                n0->obj = n.obj;
                n0->cost = n.cost;
                n0->sum_p = n.sum_p;
                n0->jobs = n.jobs;
            }

            if (n1->obj < n.obj - (T) w[j] * (n.sum_p + p[j]) +  pi[j]) {
                n1->obj = n.obj - (T) w[j] * (n.sum_p + p[j]) + pi[j];
                n1->cost = n.cost + w[j] * (n.sum_p + p[j]);
                n1->sum_p = n.sum_p + p[j];
                n1->jobs = n.jobs;
                n1->jobs.push_back(j);
            }
        }

        void initializenode(PricerInfoBDD<T> &n)
        {
            n.obj = -1871286761.0;
        }

        void initializerootnode(PricerInfoBDD<T> &n)
        {
            n.obj = pi[nbjobs];
        }

        Optimal_Solution<T> get_objective(PricerInfoBDD<T> *n)
        {
            Optimal_Solution<T> sol(n);
            return sol;
        }
};

template<typename E, typename T>
class WeightBDD: public tdzdd::DdEval<E, PricerWeightBDD<T>, Optimal_Solution<T> >
{
        T *pi;
        int *p;
        int *w;
        int *r;
        int *d;
        int nbjobs;

    public:
        WeightBDD(T *_pi, int *_p, int *_w, int *_r, int *_d, int _nbjobs)
            : pi(_pi), p(_p), w(_w), r(_r), d(_d), nbjobs(_nbjobs)
        {
        };

        void evalTerminal(PricerWeightBDD<T> &n)
        {
            n.obj = pi[nbjobs];
            n.cost = 0.0;
            n.sum_w = 0.0;
            n.sum_p = 0.0;
            n.take = false;
        }

        void evalNode(PricerWeightBDD<T> *n, int i, tdzdd::DdValues<PricerWeightBDD<T>, 2>    &values) const
        {
            int j = nbjobs - i;
            assert(j >= 0 && j <= nbjobs - 1);
            PricerWeightBDD<T> n0 = values.get(0);
            PricerWeightBDD<T> n1 = values.get(1);

            if (n0.obj >= n1.obj - (T) w[j] * (n->sum_p + p[j]) + pi[j]) {
                n->obj = n0.obj;
                n->take = false;
            } else {
                n->obj = n1.obj - (T) w[j] * (n->sum_p + p[j]) + pi[j];
                n->take = true;
            }
        }

        void initializenode(PricerWeightBDD<T> &n)
        {
            n.take = false;
        }

        Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram, tdzdd::DataTable<PricerWeightBDD<T>> *data_table, const tdzdd::NodeId *f)
        {
            Optimal_Solution<T> sol;
            sol.obj = (*data_table)[f->row()][f->col()].obj;
            tdzdd::NodeId cur_node = *f;
            sol.cost = 0;
            sol.C_max = 0;
            int j = nbjobs - cur_node.row();

            while (cur_node.row() != 0 || cur_node.col() != 0) {
                if ((*data_table)[cur_node.row()][cur_node.col()].take &&  r[j] <= sol.C_max && sol.C_max + p[j] <= d[j]) {
                    sol.jobs.push_back(j);
                    cur_node = diagram.privateEntity().child(cur_node, 1);
                    sol.C_max += p[j];
                    sol.cost += w[j] * sol.C_max;
                    j = nbjobs - cur_node.row();
                } else {
                    cur_node = diagram.privateEntity().child(cur_node, 0);
                    j = nbjobs - cur_node.row();
                }
            }
            return sol;
        }
};

template<typename E, typename T>
class WeightZDD: public tdzdd::DdEval<E, PricerWeightZDD<T>, Optimal_Solution<T> >
{
        T *pi;
        int *p;
        int *w;
        int *r;
        int *d;
        int nbjobs;
        int H_min;
        int H_max;

    public:
        WeightZDD(T *_pi, int *_p, int *_w, int *_r, int *_d, int _nbjobs, int Hmin, int Hmax)
            : pi(_pi), p(_p), w(_w), r(_r), d(_d), nbjobs(_nbjobs), H_min(Hmin), H_max(Hmax)
        {
        };

        void evalTerminal(PricerWeightZDD<T> &n)
        {
            for (int i = 0 ; i < 2000; i++) {
                n.obj[i] = pi[nbjobs];
            }
        }

        void evalNode(PricerWeightZDD<T> *n, int i, tdzdd::DdValues<PricerWeightZDD<T>, 2>    &values) const
        {
            int j = nbjobs - i;
            assert(j >= 0 && j <= nbjobs - 1);
            PricerWeightZDD<T> *n0 = values.get_ptr(0);
            PricerWeightZDD<T> *n1 = values.get_ptr(1);

            for (auto &it : (n->weight)) {
                if (n0->obj[it] >= n1->obj[it + p[j]] - (T) w[j] * (it + p[j]) + pi[j]) {
                    n->obj[it] = n0->obj[it];
                    n->take[it] = false;
                } else {
                    n->obj[it] = n1->obj[it + p[j]] - (T) w[j] * (it + p[j]) + pi[j];
                    n->take[it] = true;
                }
            }
        }

        void initializenode(PricerWeightZDD<T> &n)
        {
            for (auto &i : (n.weight)) {
                n.take[i] = false;
                n.obj[i] = 0.0;
            }
        }

        Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram, tdzdd::DataTable<PricerWeightZDD<T>> *data_table, const tdzdd::NodeId *f)
        {
            Optimal_Solution<T> sol;
            sol.C_max = 0;
            sol.obj = (*data_table)[f->row()][f->col()].obj[sol.C_max];
            sol.cost = 0;
            tdzdd::NodeId cur_node = *f;
            int j = nbjobs - cur_node.row();

            while (cur_node.row() != 0) {
                if ((*data_table)[cur_node.row()][cur_node.col()].take[sol.C_max] &&  r[j] <= sol.C_max && sol.C_max + p[j] <= d[j]) {
                    sol.jobs.push_back(j);
                    cur_node = diagram.privateEntity().child(cur_node, 1);
                    sol.C_max += p[j];
                    sol.cost += w[j] * sol.C_max;
                    j = nbjobs - cur_node.row();
                } else {
                    cur_node = diagram.privateEntity().child(cur_node, 0);
                    j = nbjobs - cur_node.row();
                }
            }

            return sol;
        }
};

/**
 * Farkas
 */
template<typename E, typename T>
class Farkas: public tdzdd::DdEval<E, PricerInfoBDD<T> >
{
        T *pi;
        int *p;
        int *w;
        int nbjobs;

    public:
        Farkas(T *_pi, int *_p, int *_w, int _nbjobs)
            : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs)
        {
        };

        void evalTerminal(PricerInfoBDD<T> &n, bool one)
        {
            n.obj = one ? 0 : -1871286761.0;
            n.cost = 0;
            n.sum_p = 0;
            n.jobs.resize(0);
        }

        void evalNode(PricerInfoBDD<T> &n, int i, tdzdd::DdValues<PricerInfoBDD<T>, 2>    &values) const
        {
            int j = nbjobs - i;
            assert(j >= 0 && j <= nbjobs - 1);
            PricerInfoBDD<T> *n0 = values.get_ptr(0);
            PricerInfoBDD<T> *n1 = values.get_ptr(1);

            if (n0->obj <= n.obj) {
                n0->obj = n.obj;
                n0->cost = n.cost;
                n0->sum_p = n.sum_p;
            }

            if (n1->obj < n.obj + pi[j]) {
                n1->obj = n.obj + pi[j];
                n1->cost = n.cost + w[j] * (n.sum_p + p[j]);
                n1->sum_p = n.sum_p + p[j];
            }
        }

        void initializenode(PricerInfoBDD<T> &n)
        {
            n.obj = -1871286761.0;
        }

        void initializerootnode(PricerInfoBDD<T> &n)
        {
            n.obj = pi[nbjobs];
        }

};

template<typename E, typename T>
class FarkasZDD: public tdzdd::DdEval<E, PricerFarkasZDD<T>, Optimal_Solution<T> >
{
        T *pi;
        int *p;
        int *w;
        int *r;
        int *d;
        int nbjobs;
        int H_min;
        int H_max;

    public:
        FarkasZDD(T *_pi, int *_p, int *_w, int *_r, int *_d, int _nbjobs, int Hmin, int Hmax)
            : pi(_pi), p(_p), w(_w), r(_r), d(_d), nbjobs(_nbjobs), H_min(Hmin), H_max(Hmax)
        {
        };

        void evalTerminal(PricerFarkasZDD<T> &n)
        {
            n.obj = pi[nbjobs];
        }

        void evalNode(PricerFarkasZDD<T> *n, int i, tdzdd::DdValues<PricerFarkasZDD<T>, 2>    &values) const
        {
            int j = nbjobs - i;
            assert(j >= 0 && j <= nbjobs - 1);
            PricerFarkasZDD<T> *n0 = values.get_ptr(0);
            PricerFarkasZDD<T> *n1 = values.get_ptr(1);

            if (n0->obj >= n1->obj + pi[j]) {
                n->obj = n0->obj;
                n->take = false;
            } else {
                n->obj = n1->obj + pi[j];
                n->take = true;
            }
        }

        void initializenode(PricerFarkasZDD<T> &n)
        {
            n.take = false;
        }

        Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram, tdzdd::DataTable<PricerFarkasZDD<T>> *data_table, const tdzdd::NodeId *f)
        {
            Optimal_Solution<T> sol;
            sol.obj = (*data_table)[f->row()][f->col()].obj;
            sol.cost = 0;
            sol.C_max = 0;
            tdzdd::NodeId cur_node = *f;
            int j = nbjobs - cur_node.row();

            while (cur_node.row() != 0) {
                if ((*data_table)[cur_node.row()][cur_node.col()].take &&  r[j] <= sol.C_max && sol.C_max + p[j] <= d[j]) {
                    sol.jobs.push_back(j);
                    cur_node = diagram.privateEntity().child(cur_node, 1);
                    sol.C_max += p[j];
                    sol.cost += w[j] * sol.C_max;
                    j = nbjobs - cur_node.row();
                } else {
                    cur_node = diagram.privateEntity().child(cur_node, 0);
                    j = nbjobs - cur_node.row();
                }
            }

            return sol;
        }
};

struct DurationZDDdouble: DurationZDD<DurationZDDdouble, double> {
    DurationZDDdouble(double *_pi, int *_p, int *_w, int _nbjobs, int H_max): DurationZDD<DurationZDDdouble,  double>(_pi, _p, _w, _nbjobs, H_max)
    {
    };
};

struct DurationBDDdouble: DurationBDD<DurationBDDdouble, double> {
    DurationBDDdouble(double *_pi, int *_p, int *_w, int _nbjobs): DurationBDD<DurationBDDdouble, double>(_pi, _p, _w, _nbjobs)
    {
    };
};

struct WeightBDDdouble: WeightBDD<WeightBDDdouble, double> {
    WeightBDDdouble(double *_pi, int *_p, int *_w, int *r, int *d, int _nbjobs): WeightBDD<WeightBDDdouble, double>(_pi, _p, _w, r, d, _nbjobs)
    {
    };
};

struct WeightZDDdouble: WeightZDD<WeightZDDdouble, double> {
    WeightZDDdouble(double *_pi, int *_p, int *_w, int *r, int *d, int _nbjobs, int Hmin, int Hmax): WeightZDD<WeightZDDdouble, double>(_pi, _p, _w, r, d, _nbjobs, Hmin, Hmax)
    {
    };
};

struct FarkasZDDdouble: FarkasZDD<FarkasZDDdouble, double> {
    FarkasZDDdouble(double *_pi, int *_p, int *_w, int *r, int *d, int _nbjobs, int Hmin, int Hmax): FarkasZDD<FarkasZDDdouble, double>(_pi, _p, _w, r, d, _nbjobs, Hmin, Hmax)
    {
    };
};
