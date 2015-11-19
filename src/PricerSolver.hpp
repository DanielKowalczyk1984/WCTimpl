#include "PricerConstruct.hpp"
#include "PricerEvaluate.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/eval/Cardinality.hpp"
#include "tdzdd/dd/Node.hpp"
#include <iostream>
#include <vector>

struct MaxNumItems: public tdzdd::DdEval<MaxNumItems,int> {
    void evalTerminal(int& n, bool one) const {
        n = one ? 0 : INT_MIN;
    }

    void evalNode(int& n, int, tdzdd::DdValues<int,2> const& values) const {
        n = std::max(values.get(0), values.get(1) + 1);
    }
};

struct MinNumItems: public tdzdd::DdEval<MinNumItems,int> {
    void evalTerminal(int& n, bool one) const {
        n = one ? 0 : INT_MAX - 1;
    }

    void evalNode(int& n, int, tdzdd::DdValues<int,2> const& values) const {
        n = std::min(values.get(0), values.get(1) + 1);
    }
};

struct IntSubset {
    virtual ~IntSubset() {
    }

    virtual bool contains(int x) const = 0;

    virtual int lowerBound() const {
        return 0;
    }

    virtual int upperBound() const {
        return INT_MAX;
    }
};

class IntRange: public IntSubset {
    int const min;
    int const max;
    int const step;

public:
    IntRange(int min = 0, int max = INT_MAX, int step = 1)
            : min(min), max(max), step(step) {
    }

    bool contains(int x) const {
        if (x < min || max < x) return false;
        return (x - min) % step == 0;
    }

    int lowerBound() const {
        return min;
    }

    int upperBound() const {
        return max;
    }
};

class SizeConstraint: public tdzdd::DdSpec<SizeConstraint,int,2> {
    int const n;
    IntSubset const* const constraint;

public:
    SizeConstraint(int n, IntSubset const& constraint)
            : n(n), constraint(&constraint) {
        assert(n >= 1);
    }

    SizeConstraint(int n, IntSubset const* constraint)
            : n(n), constraint(constraint) {
        assert(n >= 1);
    }

    int getRoot(int& count) const {
        count = 0;
        return (constraint && n < constraint->lowerBound()) ? 0 : n;
    }

    int getChild(int& count, int level, int value) const {
        if (constraint == 0) return (--level >= 1) ? level : -1;

        if (value) {
            if (count >= constraint->upperBound()) return 0;
            ++count;
        }
        else {
            if (count + level <= constraint->lowerBound()) return 0;
        }

        return (--level >= 1) ? level : constraint->contains(count) ? -1 : 0;
    }
};

class PricerSolver
{
public:
    tdzdd::DdStructure<2> dd;
    tdzdd::DdStructure<2> tmp;
    int *p;
    int *w;
    int *r;
    int *d;
    int nbjobs;
    PricerSolver( int *_p, int *_w,  int *_r, int *_d,int njobs, int Hmin, int Hmax):p(_p),w(_w),r(_r),d(_d),nbjobs(njobs){
        PricerSpec ps(p,r,d,nbjobs, Hmin,Hmax);
        dd = tdzdd::DdStructure<2>(ps);
        dd.zddReduce();
        delete [] ps.sum_p;
        delete [] ps.min_p;
    };

    struct PricerInfo<double> solveDbl(double *pi){
        return dd.evaluate(MaxReducedCostDbl(pi,p,w,nbjobs));
    }

    struct PricerInfo<int> solveInt(int* pi){
        return dd.evaluate(MaxReducedCostInt(pi, p, w, nbjobs));
    }

    struct PricerInfo<double> solvefarkasDbl(double *pi){
        return dd.evaluate(MaxFarkasPricingDbl(pi,p,w,nbjobs));
    }

    struct PricerInfo<int> solvefarkasInt(int* pi){
        return dd.evaluate(MaxFarkasPricingInt(pi, p, w, nbjobs));
    }

    void addRestriction(){

    }

    void addDuetimeConstraint(){

    }

    void addReleasetimeConstraint(){

    }

    ~PricerSolver(){
        
    };

    
};


