#include "PricerConstruct.hpp"
#include "PricerEvaluate.hpp"
#include "tdzdd/DdStructure.hpp"
#include <iostream>
#include <vector>

class PricerSolver{
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

   class PricerInfo<double> solveDbl(double *pi){
        return dd.evaluate(MaxReducedCostDbl(pi,p,w,nbjobs));
    }

    class PricerInfo<int> solveInt(int* pi){
        return dd.evaluate(MaxReducedCostInt(pi, p, w, nbjobs));
    }

    class PricerInfo<double> solvefarkasDbl(double *pi){
        return dd.evaluate(MaxFarkasPricingDbl(pi,p,w,nbjobs));
    }

    class PricerInfo<int> solvefarkasInt(int* pi){
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


