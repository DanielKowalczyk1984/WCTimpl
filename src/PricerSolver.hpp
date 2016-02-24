#include "PricerConstruct.hpp"
#include "PricerEvaluate.hpp"
#include "tdzdd/DdStructure.hpp"
#include <iostream>
#include <vector>

class PricerSolver{
public:
    tdzdd::DdStructure<2> dd;
    tdzdd::DdStructure<2> zdd;
    tdzdd::DdStructure<2> tmp;
    int *p;
    int *w;
    int *r;
    int *d;
    const int nbjobs;
    int H_max;

    PricerSolver( int *_p, int *_w,  int *_r, int *_d, int njobs, int Hmin, int Hmax, bool print = false, bool reduce = false): p(_p), w(_w), r(_r), d(_d), nbjobs(njobs), H_max(Hmax) {
        PricerSpec ps(p, r, d, nbjobs, Hmin, Hmax);
        if (print) {
            std::ofstream file;
            file.open("PricerSpec.txt");
            ps.dumpDot(file);
            file.close();
        }
        dd = tdzdd::DdStructure<2>(ps);
        printf("size of dd = %lu\n", dd.size());
        if (reduce) {
            std::cout << "Reducing the size of DD structure" << std::endl;
            zdd = dd;
            zdd.zddReduce();
            if (print) {
                create_dot_zdd("DDStructure.txt");
            }
        }
        delete [] ps.sum_p;
        delete [] ps.min_p;
    };

    void create_dot_zdd(const char* name){
        std::ofstream file;
        file.open(name);    
        dd.dumpDot(file);
        file.close();
    }

    void calculate_reducedcost_iteration(double *pi){
        
        double_t min = .0;
        for (int i = 0; i <= nbjobs; ++i)
        {
            std::cout << pi[i] << " ";
        }
        std::cout << std::endl;
        std::vector<int> v;
        int counter = 0;
        for (tdzdd::DdStructure<2>::const_iterator it = dd.begin(); it != dd.end(); ++it)
        {
            counter++;
            double_t obj = 0.0;
            int time_ = 0;
            for (std::vector<int>::const_iterator i = (*it).begin(); i != (*it).end(); ++i)
            {
                //std::cout << nbjobs - *i << " ";
                time_ += p[nbjobs - *i];
                obj += (double ) w[nbjobs - *i]* (double) time_ - pi[nbjobs - *i];
            }
            if(obj < min) {
                min = obj;
                v.clear();
                for (std::vector<int>::const_iterator i = (*it).begin(); i != (*it).end(); ++i)
                {
                    v.push_back(nbjobs - *i);
                }
            }
        }
        
        for (std::vector<int>::iterator i = v.begin(); i != v.end(); ++i)
        {
            std::cout << *i << " ";
        }
        std::cout << "min = " << min  <<  " count = " << counter << std::endl;
    }

    class PricerInfoBDD<double> solveDbl(double *pi){
        return dd.evaluate(DurationBDDdouble(pi, p, w, nbjobs));
    }

    class PricerInfoBDD<int> solveInt(int* pi){
        return dd.evaluate(DurationBDDint(pi, p, w, nbjobs));
    }

    class PricerInfoBDD<double> solvefarkasDbl(double *pi){
        return dd.evaluate(FarkasBDDdouble(pi, p, w, nbjobs));
    }

    class PricerInfoBDD<int> solvefarkasInt(int* pi){
        return dd.evaluate(FarkasBDDint(pi, p, w, nbjobs));
    }

    class PricerInfoBDD<double> solveDblBDD(double *pi){
        return dd.evaluate(DurationBDDdouble(pi, p, w, nbjobs));
    }

    class PricerInfoBDD<double> solve_duration_bdd_double(double *pi){
        return dd.evaluate_reverse(DurationBDDdouble(pi, p, w, nbjobs));
    }

    class Optimal_ZDD<double> solve_duration_zdd_double(double *pi){
        return zdd.evaluate_forward_DP(DurationZDDdouble(pi, p, w, nbjobs,H_max));
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


