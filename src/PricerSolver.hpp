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
    int H_min;
    int H_max;

    PricerSolver( int *_p, int *_w,  int *_r, int *_d, int njobs, int Hmin, int Hmax, bool print = false, bool reduce = false): p(_p), w(_w), r(_r), d(_d), nbjobs(njobs),H_min(Hmin),H_max(Hmax) {
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

    class Optimal_ZDD<double> dynamic_programming_ahv(double *pi){
        Optimal_ZDD<double> opt_sol;
        opt_sol.cost = 0;
        double **F;
        bool **A;
        int t_min = H_min;
        double min;

        F = new double* [nbjobs + 1];
        A = new bool* [nbjobs + 1];
        for(int i = 0; i < nbjobs + 1; i++) {
            F[i] = new double [H_max + 1];
            A[i] = new bool [H_max + 1];
        }
        /** Initialisation */
        F[0][0] = -pi[nbjobs];
        A[0][0] = false;
        for(int t = 1; t < H_max + 1; t++) {
            F[0][t] = 271827676.0;
            A[0][t] = false;
        }

        for(int i = 1; i < nbjobs + 1; i++) {
            for(int t = 0; t < H_max + 1; t++) {
                F[i][t] = 271827676.0;
                A[i][t] = false;
            }
        }

        /** Recursion */
        for(int i = 1; i < nbjobs + 1; i++) {
            for(int t = 0; t < H_max + 1; t++) {
                if(t >= r[i - 1] + p[i - 1] && t <= d[i - 1]) {
                    if(F[i - 1][t - p[i - 1]] + (double) w[i - 1]*t - pi[i - 1] < F[i - 1][t]) {
                        F[i][t] = F[i - 1][t - p[i - 1]] + (double) w[i - 1]*t - pi[i - 1];
                        A[i][t] = true;
                    } else {
                        F[i][t] = F[i - 1][t];
                        A[i][t] = false;
                    }
                } else {
                    F[i][t] = F[i - 1][t];
                    A[i][t] = false;
                }
            }
        }

        /** Find optimal solution */
        min = F[nbjobs][H_min + 1];
        opt_sol.obj = min;
        for(int i = H_min + 1; i < H_max + 1; ++i) {
            if(F[nbjobs][i] < min) {
                min = F[nbjobs][i];
                t_min = i;
                opt_sol.obj = min;
            }
        }

        /** Construct the solution */
        for(int i = nbjobs; i >= 1; --i) {
            if(A[i][t_min] && r[i - 1] + p[i -1] <= t_min && t_min <= d[i - 1]) {   
                opt_sol.jobs.push_back(i - 1);
                opt_sol.cost += w[i - 1]*t_min;
                t_min -= p[i - 1];
            }
        }

        /** Free the memory */
        for(int i = 0; i < nbjobs + 1; ++i) {
            delete[] A[i];
            delete[] F[i];
        }
        delete[] A;
        delete[] F;
        return opt_sol;
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


