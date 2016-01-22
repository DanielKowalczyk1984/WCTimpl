#include "tdzdd/DdEval.hpp"
#include <cfloat>
#include <climits>
#include <vector>


template<typename T>
class PricerInfo {
public:
    int sum_w;
    T obj;
    std::vector<int> jobs;
    int cost;

    PricerInfo& operator=( const PricerInfo& other ) {
        sum_w = other.sum_w;
        obj = other.obj;
        jobs = other.jobs;
        cost = other.cost;
        return *this;
    };
};

struct MinNumItems: public tdzdd::DdEval<MinNumItems,int> {
    void evalTerminal(int& n, bool one) const {
        n = one ? 0 : INT_MAX - 1;
    }

    void evalNode(int& n, int, int &v0, int i0, int &v1, int i1) const {
        n = std::min(v0, v1 + 1);
    }
};

struct MaxNumItems: public tdzdd::DdEval<MaxNumItems,int> {
    void evalTerminal(int& n, bool one) const {
        n = one ? 0 : INT_MIN;
    }

    void evalNode(int& n, int, int &v0, int i0, int &v1, int i1) const {
        n = std::max(v0, v1 + 1);
    }
};

template<typename E>
class MaxReducedCostBaseDbl: public tdzdd::DdEval<E, PricerInfo<double> > {
    double *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxReducedCostBaseDbl(double *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    };

    void evalTerminal( PricerInfo<double>& n, bool one) {
        n.obj = one ? pi[nbjobs] : DBL_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);
    }

    void evalNode( PricerInfo<double> &n, int i, tdzdd::DdValues<PricerInfo<double>, 2> const & values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfo<double> n0 = values.get(0);
        PricerInfo<double> n1 = values.get(1);
        if (n0.obj > n1.obj - n1.sum_w * p[j] - (double)w[j] * (double)p[j] + pi[j]) {
            n.obj = n0.obj;
            n.sum_w = n0.sum_w;
            n.jobs = n0.jobs;
            n.cost = n0.cost;
        } else {
            n.obj = n1.obj - n1.sum_w * p[j] - w[j] * p[j] + pi[j];
            n.sum_w = n1.sum_w + w[j];
            n.jobs = n1.jobs;
            n.cost = n1.cost + n1.sum_w * p[j] + w[j] * p[j];
            n.jobs.push_back(j);
        }

    }
};

template<typename E>
class MaxFarkasPricingBaseDbl: public tdzdd::DdEval<E, PricerInfo<double> > {
    double *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxFarkasPricingBaseDbl(double *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {
    };

    void evalTerminal( PricerInfo<double>& n, bool one) {
        n.obj = one ? pi[nbjobs] : DBL_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }



    void evalNode( PricerInfo<double> &n, int i, tdzdd::DdValues<PricerInfo<double>, 2> const & values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfo<double> n0 = values.get(0);
        PricerInfo<double> n1 = values.get(1);
        if (n0.obj > n1.obj + pi[j]) {
            n.obj = n0.obj;
            n.sum_w = n0.sum_w;
            n.jobs = n0.jobs;
            n.cost = n0.cost;
        } else {
            n.obj = n1.obj  + pi[j];
            n.sum_w = n1.sum_w + w[j];
            n.jobs = n1.jobs;
            n.cost = n1.cost + n1.sum_w * p[j] + w[j] * p[j];
            n.jobs.push_back(j);
        }

    }
};

template<typename E>
class MaxReducedCostBaseInt: public tdzdd::DdEval<E, PricerInfo<int> > {
    int *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxReducedCostBaseInt(int *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    };

    void evalTerminal( PricerInfo<int>& n, bool one) {
        n.obj = one ? pi[nbjobs] : INT_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }

    void evalNode( PricerInfo<int> &n, int i, tdzdd::DdValues<PricerInfo<int>, 2> const & values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfo<int> n0 = values.get(0);
        PricerInfo<int> n1 = values.get(1);
        if (n0.obj > n1.obj - n1.sum_w * p[j] - w[j]*p[j] + pi[j]) {
            n.obj = n0.obj;
            n.sum_w = n0.sum_w;
            n.jobs = n0.jobs;
            n.cost = n0.cost;
        } else {
            n.obj = n1.obj - n1.sum_w * p[j] - w[j] * p[j] + pi[j];
            n.sum_w = n1.sum_w + w[j];
            n.jobs = n1.jobs;
            n.cost = n1.cost + n1.sum_w * p[j] + w[j] * p[j];
            n.jobs.push_back(j);
        }
    }
};

template<typename E>
class MaxFarkasPricingBaseInt: public tdzdd::DdEval<E, PricerInfo<int> > {
    int *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxFarkasPricingBaseInt(int *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    };

    void evalTerminal( PricerInfo<int>& n, bool one) {
        n.obj = one ? pi[nbjobs] : INT_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }

    void evalNode( PricerInfo<int> &n, int i, tdzdd::DdValues<PricerInfo<int>, 2> const & values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerInfo<int> n0 = values.get(0);
        PricerInfo<int> n1 = values.get(1);
        if (n0.obj > n1.obj + pi[j]) {
            n.obj = n0.obj;
            n.sum_w = n0.sum_w;
            n.jobs = n0.jobs;
            n.cost = n0.cost;
        } else {
            n.obj = n1.obj + pi[j];
            n.sum_w = n1.sum_w + w[j];
            n.jobs = n1.jobs;
            n.cost = n1.cost + n1.sum_w * p[j] + w[j] * p[j];
            n.jobs.push_back(j);
        }
    }
};



struct MaxReducedCostDbl: MaxReducedCostBaseDbl<MaxReducedCostDbl> {
    MaxReducedCostDbl(double *_pi, int *_p, int *_w, int _nbjobs): MaxReducedCostBaseDbl<MaxReducedCostDbl>(_pi, _p, _w, _nbjobs) {
    };
};

struct MaxReducedCostInt: public MaxReducedCostBaseInt<MaxReducedCostInt> {
    MaxReducedCostInt(int *_pi, int *_p, int *_w, int _nbjobs): MaxReducedCostBaseInt<MaxReducedCostInt>(_pi, _p, _w, _nbjobs) {
    };
};

struct MaxFarkasPricingDbl: public MaxFarkasPricingBaseDbl<MaxFarkasPricingDbl> {
    MaxFarkasPricingDbl(double *_pi, int *_p, int *_w, int _nbjobs): MaxFarkasPricingBaseDbl<MaxFarkasPricingDbl>(_pi, _p, _w, _nbjobs) {
    };
};

struct MaxFarkasPricingInt: public MaxFarkasPricingBaseInt<MaxFarkasPricingInt> {
    MaxFarkasPricingInt(int *_pi, int *_p, int *_w, int _nbjobs): MaxFarkasPricingBaseInt<MaxFarkasPricingInt>(_pi, _p, _w, _nbjobs) {
    };
};

