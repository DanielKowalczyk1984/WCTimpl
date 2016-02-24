#include "tdzdd/DdEval.hpp"
#include <cfloat>
#include <climits>
#include <vector>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

/**
 * ZDD
 */
template<typename T>
class PricerInfoZDD {
public:
    T  *obj;
    int *cost;
    boost::dynamic_bitset<>  *A;

    PricerInfoZDD():obj(0),cost(0),A(0){

    }

    ~PricerInfoZDD(){
        if(obj != 0) {
            delete[] obj;
            obj = 0;
        }
        if(cost != 0) {
            delete[] cost;
            cost = 0;
        }
        if(A != 0) {
            delete[] A;
            A = 0;
        }
    }

    size_t get_max(const int &L){
        T max = obj[0];
        size_t it_max = 0;
        for(size_t i = 0; i < L; ++i) {
            if(max < obj[i]) {
                max = obj[i];
                it_max = i;
            }
        }
        return it_max;
    }

    PricerInfoZDD& operator=( const PricerInfoZDD& other ) {
        obj = other.obj;
        cost = other.cost;
        A = other.A;
        return *this;
    };
};

template<typename T>
class Optimal_ZDD{
    public:
        T obj;
        std::vector<int> jobs;
        int cost;
        Optimal_ZDD(PricerInfoZDD<T> *node, const int L){
            int it_max = node->get_max(L);
            boost::dynamic_bitset<> *ptr_max = &node->A[it_max];
            size_t it = ptr_max->find_first();
            obj = node->obj[it_max];
            cost = node->cost[it_max];
            while(it != boost::dynamic_bitset<>::npos){
                if((*ptr_max)[it]) {
                    jobs.push_back(it);
                }
                it = ptr_max->find_next(it);
            }
        }

        Optimal_ZDD(){

        }

        Optimal_ZDD& operator=(const Optimal_ZDD& other){
            obj = other.obj;
            cost = other.cost;
            jobs = other.jobs;
        }
};

template<typename E, typename T>
class DurationZDD: public tdzdd::DdEval<E, PricerInfoZDD<T>, Optimal_ZDD<T> > {
    T *pi;
    int *p;
    int *w;
    int nbjobs;
    int L;

public:
    DurationZDD(T *_pi, int *_p, int *_w, int _nbjobs, int _L)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs), L(_L) {
    };

    void evalTerminal( PricerInfoZDD<T>& n, bool one) {
        for(unsigned i = 0; i < L + 1; ++i) {
            n.obj[i] = one ? 0 : -15656873.0;
        }
    }

    void evalNode( PricerInfoZDD<T>  &n, int i, tdzdd::DdValues<PricerInfoZDD<T>, 2>  &  values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfoZDD<T> *n0 = values.get_ptr(0);
        PricerInfoZDD<T> *n1 = values.get_ptr(1);

        for(int k = 0; k < L + 1; ++k) {
            int it = k + p[j];
            if(it <= L) {
                if(n1->obj[it] < n.obj[k] - w[j]*it + pi[j]) {
                    n1->obj[it] = n.obj[k] - w[j]*it + pi[j];
                    n1->cost[it] = n.cost[k] + w[j]*it;
                    n1->A[it].clear();
                    n1->A[it] = n.A[k];
                    n1->A[it][j] = true;
                }
            }
            
            if(n0->obj[k] < n.obj[k]) {
                n0->obj[k] = n.obj[k];
                n0->cost[k] = n.cost[k];
                n0->A[k] = n.A[k];
            }
        }
    }

    void initializenode(PricerInfoZDD<T> &n) {
        n.obj = new T [L + 1];
        n.cost = new int [L + 1];
        n.A = new boost::dynamic_bitset<> [L + 1];
        for(unsigned i = 0; i < L + 1; ++i) {
            n.obj[i] = -10983290.0;
            n.A[i].resize(nbjobs);
            n.cost[i] = 0;
        }

    }

    void initializerootnode(PricerInfoZDD<T> &n) {
        n.obj = new T [L + 1];
        n.cost = new int [L + 1];
        n.A = new boost::dynamic_bitset<> [L + 1];
        for(unsigned i = 0; i < L + 1; ++i) {
            n.obj[i] = pi[nbjobs];
            n.A[i].resize(nbjobs);
            n.cost[i] = 0;
        }
    }

    Optimal_ZDD<T> get_objective(PricerInfoZDD<T> *n){
        Optimal_ZDD<T> sol(n, L);
        return sol;
    }
};

template<typename E, typename T>
class WeightZDD: public tdzdd::DdEval<E, PricerInfoZDD<T>, T> {
    T *pi;
    int *p;
    int *w;
    int nbjobs;
    int L;

public:
    WeightZDD(T *_pi, int *_p, int *_w, int _nbjobs, int _L)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs), L(_L) {
    };

    void evalTerminal( PricerInfoZDD<T>& n, bool one) {
        n.obj = one ? 0 : -1871286761.0;
        n.cost = 0;
        n.sum_w = 0;
        n.jobs.resize(0);
    }

    void evalNode( PricerInfoZDD<T> &n, int i, tdzdd::DdValues<PricerInfoZDD<T>, 2>  &  values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfoZDD<T> *n0 = values.get_ptr(0);
        PricerInfoZDD<T> *n1 = values.get_ptr(1);

        for(unsigned i = 0; i < L + 1; ++i) {
            if(i >= w[j]) {
                if(n0->obj[i - w[j]] > n1->obj[i - w[j]] - p[j]*i + pi[j]) {
                    n.obj[i - w[j]] = n0->obj[j];
                    n.A[i - w[j]][j] = false;
                } else {
                    n.obj[i] = n1->obj[i] - p[j]*i + pi[j];
                    n.A[i][j] = true;
                }
            } else {
                n.obj[i] = n0->obj[i];
                n.A[i][j] = false;
            }
        }
    }

    void initializenode(PricerInfoZDD<T> &n) {
        n.obj = new T [L + 1];
        n.cost = new int [L + 1];
        n.A = new boost::dynamic_bitset<> [L + 1];
        for(unsigned i = 0; i < L + 1; ++i) {
            n.obj[i] = -10983290.0;
            n.A[i].resize(nbjobs);
        }
    }

    /*T get_objective(PricerInfoZDD<T> &n){
        return n.get_max(L);
    }*/
};

/**
 * BDD
 */
template<typename T>
class PricerInfoBDD {
public:
    T obj;
    std::vector<int> jobs;
    int sum_w;
    int sum_p;
    int cost;

    PricerInfoBDD& operator=( const PricerInfoBDD& other ) {
        sum_w = other.sum_w;
        sum_p = other.sum_p;
        cost = other.cost;
        obj = other.obj;
        jobs.clear();
        jobs = other.jobs;
        return *this;
    };

    friend std::ostream& operator<<(std::ostream& os, PricerInfoBDD<T> const& o) {
        os << "max = " << o.obj << "," << std::endl << "cost = " << o.cost << std::endl;
        return os;
    }
};


template<typename E, typename T>
class DurationBDD: public tdzdd::DdEval<E, PricerInfoBDD<T> > {
    T *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    DurationBDD(T *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    };

    void evalTerminal( PricerInfoBDD<T>& n, bool one) {
        n.obj = one ? 0 : -1871286761.0;
        n.cost = 0;
        n.sum_p = 0;
        n.jobs.resize(0);
    }

    void evalNode( PricerInfoBDD<T> &n, int i, tdzdd::DdValues<PricerInfoBDD<T>, 2>  &  values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfoBDD<T> *n0 = values.get_ptr(0);
        PricerInfoBDD<T> *n1 = values.get_ptr(1);

        if (n0->obj < n.obj) {
            n0->obj = n.obj;
            n0->cost = n.cost;
            n0->sum_p = n.sum_p;
            n0->jobs.clear();
            n0->jobs = n.jobs;
        }

        if (n1->obj < n.obj - (T) w[j] * (n.sum_p + p[j]) +  pi[j]) {
            n1->obj = n.obj - (T) w[j] * (n.sum_p + p[j]) + pi[j];
            n1->cost = n.cost + w[j] * (n.sum_p + p[j]);
            n1->sum_p = n.sum_p + p[j];
            n1->jobs.clear();
            n1->jobs = n.jobs;
            n1->jobs.push_back(j);
        }
    }

    void initializenode(PricerInfoBDD<T> &n) {
        n.obj = -1871286761.0;
    }

    void initializerootnode(PricerInfoBDD<T> &n) {
        n.obj = pi[nbjobs];
    }
};

template<typename E, typename T>
class WeightBDD: public tdzdd::DdEval<E, PricerInfoBDD<T> > {
    T *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    WeightBDD(T *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {
    };

    void evalTerminal( PricerInfoBDD<T>& n, bool one) {
        n.obj = one ? 0 : -1871286761.0;
        n.cost = 0;
        n.sum_w = 0;
        n.jobs.resize(0);
    }

    void evalNode( PricerInfoBDD<T> &n, int i, tdzdd::DdValues<PricerInfoBDD<T>, 2>  &  values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfoBDD<T> *n0 = values.get_ptr(0);
        PricerInfoBDD<T> *n1 = values.get_ptr(1);

        if(n0->obj > n1->obj - p[j]*( n1->sum_w + w[j]) + pi[j]) {
            n = *n0;
        } else {
            n.obj = n1->obj - (T)p[j]*( n1->sum_w + w[j]) + pi[j];
            n.sum_p = n1->sum_p + p[j];
            n.sum_w = n1->sum_w + w[j];
            n.cost = n1->cost + (T)p[j]*( n1->sum_w + w[j]);
            n.jobs.clear();
            n.jobs = n1->jobs;
            n.jobs.push_back(j);
        }
    }
};

/**
 * Farkas
 */
template<typename E, typename T>
class Farkas: public tdzdd::DdEval<E, PricerInfoBDD<T> > {
    T *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    Farkas(T *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {
    };

    void evalTerminal( PricerInfoBDD<T>& n, bool one) {
        n.obj = one ? 0 : -1871286761.0;
        n.cost = 0;
        n.sum_p = 0;
        n.jobs.resize(0);
    }

    void evalNode( PricerInfoBDD<T> &n, int i, tdzdd::DdValues<PricerInfoBDD<T>, 2>  &  values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfoBDD<T> *n0 = values.get_ptr(0);
        PricerInfoBDD<T> *n1 = values.get_ptr(1);

        if (n0->obj < n.obj) {
            n0->obj = n.obj;
            n0->cost = n.cost;
            n0->sum_p = n.sum_p;
            n0->jobs.clear();
            n0->jobs = n.jobs;
        }

        if (n1->obj < n.obj - w[j] * (n.sum_p + p[j]) + pi[j]) {
            n1->obj = n.obj + pi[j];
            n1->cost = n.cost + w[j] * (n.sum_p + p[j]);
            n1->sum_p = n.sum_p + p[j];
            n1->jobs.clear();
            n1->jobs = n.jobs;
            n1->jobs.push_back(j);
        }
    }

    void initializenode(PricerInfoBDD<T> &n) {
        n.obj = -1871286761.0;
    }

    void initializerootnode(PricerInfoBDD<T> &n) {
        n.obj = pi[nbjobs];
    }

};

struct DurationZDDdouble: DurationZDD<DurationZDDdouble, double> {
    DurationZDDdouble(double *_pi, int *_p, int *_w, int _nbjobs, int H_max): DurationZDD<DurationZDDdouble,  double>(_pi, _p, _w, _nbjobs, H_max) {
    };
};

struct DurationBDDdouble: DurationBDD<DurationBDDdouble, double> {
    DurationBDDdouble(double *_pi, int *_p, int *_w, int _nbjobs): DurationBDD<DurationBDDdouble, double>(_pi, _p, _w, _nbjobs) {
    };
};

struct DurationBDDint: DurationBDD<DurationBDDint, int> {
    DurationBDDint(int *_pi, int *_p, int *_w, int _nbjobs): DurationBDD<DurationBDDint, int>(_pi, _p, _w, _nbjobs) {
    };
};

struct FarkasBDDdouble: Farkas<FarkasBDDdouble, double> {
    FarkasBDDdouble(double *_pi, int *_p, int *_w, int _nbjobs): Farkas<FarkasBDDdouble, double>(_pi, _p, _w, _nbjobs) {
    };
};

struct FarkasBDDint: Farkas<FarkasBDDint, int> {
    FarkasBDDint(int *_pi, int *_p, int *_w, int _nbjobs): Farkas<FarkasBDDint, int>(_pi, _p, _w, _nbjobs) {
    };
};
