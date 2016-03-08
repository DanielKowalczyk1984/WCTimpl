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
    int *delta;
    const int nbjobs;
    int H_min;
    int H_max;
    tdzdd::DataTable<PricerWeightBDD<double>> zdd_table;
    tdzdd::DataTable<PricerWeightZDD<double>> dd_table;

    PricerSolver( int *_p, int *_w,  int *_r, int *_d, int njobs, int Hmin, int Hmax, bool print = false, bool reduce = false): p(_p), w(_w), r(_r), d(_d), nbjobs(njobs),H_min(Hmin),H_max(Hmax) {
        delta = new int [nbjobs];
        int *sum_p = new int [nbjobs];
        sum_p[0] = p[0];
        int max_d = d[0];
        delta[0] = std::min(sum_p[0],max_d);
        for (int i = 1; i < nbjobs; i++){
            sum_p[i] = sum_p[i - 1] + p[i];
            if (max_d < d[i])
            {
                max_d = d[i];
            }
            delta[i] = std::min(sum_p[i],max_d);
        }

        PricerSpec ps(p, r, d, nbjobs, Hmin, Hmax);
        if (print) {
            std::ofstream file;
            file.open("PricerSpec.txt");
            ps.dumpDot(file);
            file.close();
        }
        dd = tdzdd::DdStructure<2>(ps);
        if (reduce) {
            zdd = dd;
            zdd.zddReduce();
            std::cout << "Reducing the size of DD structure:" <<  std::endl;
            std::cout << "DD = " << dd.size() << " " << "ZDD = " << zdd.size() << std::endl;
            if (print) {
                create_dot_zdd("DDStructure.txt");
            }
        }
        delete [] sum_p;
        delete [] ps.sum_p;
        delete [] ps.min_p;
    };

    void create_dot_zdd(const char* name){
        std::ofstream file;
        file.open(name);    
        zdd.dumpDot(file);
        file.close();
    }

    void init_table(){
        
    }

    class Optimal_Solution<double> dynamic_programming_ahv(double *pi){
        Optimal_Solution<double> opt_sol;
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
            F[0][t] = DBL_MAX/2;
            A[0][t] = false;
        }

        for(int i = 1; i < nbjobs + 1; i++) {
            for(int t = 0; t < H_max + 1; t++) {
                F[i][t] = DBL_MAX/2;
                A[i][t] = false;
            }
        }

        /** Recursion */
        for(int i = 1; i < nbjobs + 1; i++) {
            int j = i - 1;
            for(int t = 0; t < H_max; t++) {
                if(t >= r[j] + p[j] && t <= d[j]) {
                    if(F[j][t - p[j]] + (double) w[j]*t - pi[j] < F[j][t]) {
                        F[i][t] = F[j][t - p[j]] + (double) w[j]*t - pi[j];
                        A[i][t] = true;
                    } else {
                        F[i][t] = F[j][t];
                        A[i][t] = false;
                    }
                } else {
                    F[i][t] = F[j][t];
                    A[i][t] = false;
                }
            }
        }

        /** Find optimal solution */
        min = F[nbjobs][0];
        opt_sol.obj = min;
        for(int i =  H_min; i < H_max + 1; i++) {
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



    class Optimal_Solution<double> solve_duration_bdd_double(double *pi){
        return dd.evaluate_reverse(DurationBDDdouble(pi, p, w, nbjobs));
    }

    class Optimal_Solution<double> solve_duration_zdd_double(double *pi){
        return zdd.evaluate_forward_DP(DurationZDDdouble(pi, p, w, nbjobs,H_max));
    }

    class Optimal_Solution<double> solve_weight_bdd_double(double *pi){
        return dd.evaluate_weight(WeightBDDdouble(pi,p,w,r,d,nbjobs));
    }

    class Optimal_Solution<double> solve_weight_zdd_double(double *pi){
        return zdd.evaluate_weight_ZDD(WeightZDDdouble(pi,p,w,r,d,nbjobs,H_min,H_max));
    }

    void addRestriction(){

    }

    void addDuetimeConstraint(){

    }

    void addReleasetimeConstraint(){

    }

    ~PricerSolver(){
        delete [] delta;
    };

    
};


