#include "tdzdd/DdEval.hpp"
#include <cfloat>
#include <climits>
#include <vector>
#include <algorithm>


template<typename T>
class PricerInfoZDD {
public:
    T  *obj;
    int *cost;
    bool *A;

    PricerInfoZDD(){
    }

    ~PricerInfoZDD(){
        delete[] obj;
        delete[] cost;
        delete[] A;
    }

    T get_max(int H_max){
        return *std::max_element(obj, obj + (H_max + 1));
    }

    PricerInfoZDD& operator=( const PricerInfoZDD& other ) {
        obj = other.obj;
        cost = other.cost;
        A = other.A;
        return *this;
    };

    friend std::ostream& operator<<(std::ostream& os, PricerInfoZDD<T> const& o) {
        os << "max = " << o.obj << "," << std::endl << "cost = " << o.cost << std::endl;
        return os;
    }
};

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
class DurationZDD: public tdzdd::DdEval<E, PricerInfoZDD<T>, T > {
    T *pi;
    int *p;
    int *w;
    int nbjobs;
    int H_max;

public:
    DurationZDD(T *_pi, int *_p, int *_w, int _nbjobs, int _H_max)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs), H_max(_H_max) {
    };

    void evalTerminal( PricerInfoZDD<T>& n, bool one) {
        for(unsigned i = 0; i < H_max + 1; ++i) {
            n.obj[i] = one ? 0 : -15656873.0;
        }
    }

    void evalNode( PricerInfoZDD<T>  &n, int i, tdzdd::DdValues<PricerInfoZDD<T>, 2>  &  values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfoZDD<T> *n0 = values.get_ptr(0);
        PricerInfoZDD<T> *n1 = values.get_ptr(1);

        for(int k = 0; k < H_max + 1; ++k) {
            int it = k + p[j];
            if(it <= H_max) {
                if(n1->obj[it] < n.obj[k] - w[j]*it + pi[j]) {
                    n1->obj[it] = n.obj[k] - w[j]*it + pi[j];
                    n1->cost[it] = n.cost[k] + w[j]*it;
                    n1->A[it] = true;
                }

            }
            
            if(n0->obj[k] < n.obj[k]) {
                n0->obj[k] = n.obj[k];
                n0->cost[k] = n.cost[k];
                n0->A[k] = false;
            }
        }
    }

    void initializenode(PricerInfoZDD<T> &n) {
        n.obj = new T [H_max + 1];
        n.cost = new int [H_max + 1];
        n.A = new bool [H_max + 1];
        for(unsigned i = 0; i < H_max + 1; ++i) {
            n.obj[i] = -10983290.0;
        }

    }

    void initializerootnode(PricerInfoZDD<T> &n) {
        n.obj = new T [H_max + 1];
        n.cost = new int [H_max + 1];
        n.A = new bool [H_max + 1];
        for(unsigned i = 0; i < H_max + 1; ++i) {
            n.obj[i] = pi[nbjobs];
        }
    }

    T get_objective(PricerInfoZDD<T> &n){
        return n.get_max(H_max);
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
