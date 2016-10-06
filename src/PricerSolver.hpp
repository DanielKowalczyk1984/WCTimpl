#include "PricerConstruct.hpp"
#include "PricerEvaluate.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/op/Lookahead.hpp"
#include <iostream>
#include <vector>

class PricerSolver {
  public:
    tdzdd::DdStructure<2> *dd;
    tdzdd::DdStructure<2> *zdd;
    int *p;
    int *w;
    int *r;
    int *d;
    int nbjobs;
    int H_min;
    int H_max;
    tdzdd::DataTable<PricerWeightZDD<double> > zdd_table;
    tdzdd::DataTable<PricerWeightBDD<double> > dd_table;
    tdzdd::DataTable<PricerFarkasZDD<double> > farkas_table;
    bool use_zdd;

    PricerSolver(int *_p, int *_w,  int *_r, int *_d, int &njobs, int &Hmin, int &Hmax,
                 bool _use_zdd = true): p(_p), w(_w), r(_r), d(_d), nbjobs(njobs), H_min(Hmin),
        H_max(Hmax), use_zdd(_use_zdd) {
        if (use_zdd) {
            PricerSpec ps(p, r, d, nbjobs, Hmin, Hmax);
            dd = new tdzdd::DdStructure<2>(ps);
            zdd = new tdzdd::DdStructure<2>;
            *zdd = *dd;
            zdd->zddReduce();
            printf("size BDD = %lu, size ZDD= %lu\n", dd->size(), zdd->size());
            init_zdd_table();
            init_bdd_table();
            init_table_farkas();
            delete [] ps.sum_p;
            delete [] ps.min_p;
        }
    };


    PricerSolver(const PricerSolver &other) {
        dd = new tdzdd::DdStructure<2>;
        zdd = new tdzdd::DdStructure<2>;
        *dd = *(other.dd);
        *zdd = *(other.zdd);
        p = other.p;
        w = other.w;
        r = other.r;
        d = other.d;
        nbjobs = other.nbjobs;
        H_min = other.H_min;
        H_max = other.H_max;
        zdd_table.init();
        dd_table.init();
        farkas_table.init();
        use_zdd = true;
    }

    void init_tables() {
        init_bdd_table();
        init_zdd_table();
        init_table_farkas();
    }

    ~PricerSolver() {
        if (use_zdd) {
            delete dd;
            delete zdd;
        }
    }

    void create_dot_zdd(const char *name) {
        std::ofstream file;
        file.open(name);
        zdd->dumpDot(file);
        file.close();
    }

    void init_table_farkas() {
        tdzdd::NodeTableHandler<2> &node_handler_farkas = zdd->getDiagram();
        farkas_table.init(nbjobs + 1);
        size_t const m = node_handler_farkas.privateEntity()[0].size();
        farkas_table[0].resize(m);

        for (size_t i = 0; i < m; ++i) {
            farkas_table[0][i].init_terminal_node(i);
        }

        for (int i = 1; i <= nbjobs; ++i) {
            tdzdd::MyVector<tdzdd::Node<2> > const &node =
                node_handler_farkas.privateEntity()[i];
            size_t const mm = node.size();
            farkas_table[i].resize(mm);

            for (size_t j = 0; j < mm; ++j) {
                farkas_table[i][j].init_node();
            }
        }
    }

    void init_bdd_table() {
        tdzdd::NodeTableHandler<2> &handler = dd->getDiagram();
        tdzdd::NodeId &root = dd->root();
        dd_table.init(root.row() + 1);

        /** init table */
        for (int i = root.row(); i >= 0 ; i--) {
            tdzdd::MyVector<tdzdd::Node<2>> const &node = handler.privateEntity()[i];
            size_t const m = node.size();
            dd_table[i].resize(m);
        }

        /** init root */
        dd_table[root.row()][root.col()].init_node(0);

        for (size_t i = root.row(); i > 0 ; i--) {
            size_t const m = dd_table[i].size();
            int cur_job = nbjobs - i;

            for (size_t j = 0; j < m; j++) {
                int sum_p = dd_table[i][j].sum_p;
                tdzdd::NodeId cur_node = handler.privateEntity().child(i, j, 0);

                if (cur_node.row() != 0) {
                    dd_table[cur_node.row()][cur_node.col()].init_node(sum_p);
                }

                cur_node = handler.privateEntity().child(i, j, 1);

                if (cur_node.row() != 0) {
                    dd_table[cur_node.row()][cur_node.col()].init_node(sum_p + p[cur_job]);
                }
            }
        }

        /** init terminal nodes */
        size_t const mm = handler.privateEntity()[0].size();

        for (size_t j  = 0; j < mm; j++) {
            dd_table[0][j].init_terminal_node(j);
        }
    }

    void init_zdd_table() {
        tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
        tdzdd::NodeId &root = zdd->root();
        zdd_table.init(root.row() + 1);

        /** init table */
        for (int i = root.row(); i >= 0 ; i--) {
            tdzdd::MyVector<tdzdd::Node<2>> const &node = handler.privateEntity()[i];
            size_t const m = node.size();
            zdd_table[i].resize(m);
        }

        /** init root */
        zdd_table[root.row()][root.col()].add_weight(0,0);

        for (int i = root.row(); i > 0 ; i--) {
            size_t const m = zdd_table[i].size();
            int cur_job = nbjobs - i;

            for (size_t j = 0; j < m; j++) {
                tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
                tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

                // for (auto &it : zdd_table[i][j].info_node) {
                //     zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight(it.first);
                //     zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight(it.first + p[cur_job]);
                // }

                for( my_iterator<double> it = zdd_table[i][j].list.begin(); it != zdd_table[i][j].list.end();it++){
                        (*it)->n = zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight((*it)->weight, nbjobs - cur_node_0.row());
                        (*it)->y = zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight((*it)->weight + p[cur_job], nbjobs - cur_node_1.row());
                }
            }
        }

        /** init terminal nodes */
        size_t const mm = handler.privateEntity()[0].size();

        for (size_t j = 0; j < mm ; j++) {
            zdd_table[0][j].init_terminal_node(j);
        }
    }

    void init_bdd_one_conflict(int v1, int v2, int same) {
        int ecount_same = 0;
        int ecount_diff = 0;
        int *elist_same = (int *) NULL;
        int *elist_differ = (int *) NULL;

        if (same) {
            ecount_same = 1;
            elist_same = new int[2];
            elist_same[0] = v1;
            elist_same[1] = v2;
        } else {
            ecount_diff = 1;
            elist_differ = new int[2];
            elist_differ[0] = v1;
            elist_differ[1] = v2;
        }

        ConflictConstraints conflict(nbjobs, elist_same, ecount_same, elist_differ,
                                     ecount_diff);
        dd->zddSubset(conflict);
        *zdd = *dd;
        zdd->zddReduce();
        delete [] elist_same;
        delete [] elist_differ;
    }

    void init_zdd_one_conflict(int v1, int v2, int same) {
        int ecount_same = 0;
        int ecount_diff = 0;
        int *elist_same = (int *) NULL;
        int *elist_differ = (int *) NULL;

        if (same) {
            ecount_same = 1;
            elist_same = new int[2];
            elist_same[0] = v1;
            elist_same[1] = v2;
        } else {
            ecount_diff = 1;
            elist_differ = new int[2];
            elist_differ[0] = v1;
            elist_differ[1] = v2;
        }

        ConflictConstraints conflict(nbjobs, elist_same, ecount_same, elist_differ,
                                     ecount_diff);
        zdd->zddSubset(tdzdd::ZddLookahead<ConflictConstraints>(conflict));
        zdd->zddReduce();
        delete [] elist_same;
        delete [] elist_differ;
    }

    void init_bdd_conflict_solver(int *elist_same, int ecount_same,
                                  int *elist_differ, int ecount_differ) {
        if (ecount_differ + ecount_same > 0) {
            dd = new tdzdd::DdStructure<2>;
            //*dd = root_dd;
            ConflictConstraints conflict(nbjobs, elist_same, ecount_same, elist_differ,
                                         ecount_differ);
            dd->zddSubset(conflict);
            zdd = new tdzdd::DdStructure<2>;
            *zdd = *dd;
            zdd->zddReduce();
        } else {
        }

        if (dd->size() == 0) {
            return;
        }

        init_bdd_table();
        init_table_farkas();
    }

    void init_zdd_conflict_solver(int *elist_same, int ecount_same,
                                  int *elist_differ, int ecount_differ) {
        if (ecount_same + ecount_differ > 0) {
            zdd = new tdzdd::DdStructure<2>;
            ConflictConstraints conflict(nbjobs, elist_same, ecount_same, elist_differ,
                                         ecount_differ);
            zdd->zddSubset(conflict);
            zdd->zddReduce();
        } else {
        }

        init_zdd_table();
        init_table_farkas();
    }

    void free_bdd_solver(int ecount_same, int ecount_differ) {
        if (ecount_same + ecount_differ > 0) {
            delete dd;
            delete zdd;
        }

        dd_table.init();
        farkas_table.init();
    }

    void free_zdd_solver(int ecount_same, int ecount_differ) {
        if (ecount_same + ecount_differ > 0) {
            delete zdd;
        }

        zdd_table.init();
        farkas_table.init();
    }

    class Optimal_Solution<double> dynamic_programming_ahv(double *pi) {
        Optimal_Solution<double> opt_sol;
        opt_sol.cost = 0;
        double **F;
        bool **A;
        int t_min = H_min;
        F = new double* [nbjobs + 1];
        A = new bool* [nbjobs + 1];

        for (int i = 0; i < nbjobs + 1; i++) {
            F[i] = new double [H_max + 1];
            A[i] = new bool [H_max + 1];
        }

        /** Initialisation */
        F[0][0] = -pi[nbjobs];
        A[0][0] = false;

        for (int t = 1; t < H_max + 1; t++) {
            F[0][t] = DBL_MAX / 2;
            A[0][t] = false;
        }

        for (int i = 1; i < nbjobs + 1; i++) {
            for (int t = 0; t < H_max + 1; t++) {
                F[i][t] = DBL_MAX / 2;
                A[i][t] = false;
            }
        }

        /** Recursion */
        for (int i = 1; i < nbjobs + 1; i++) {
            int j = i - 1;

            for (int t = 0; t < H_max; t++) {
                if (t >= r[j] + p[j] && t <= d[j]) {
                    if (F[j][t - p[j]] + (double) w[j]*t - pi[j] < F[j][t]) {
                        F[i][t] = F[j][t - p[j]] + (double) w[j] * t - pi[j];
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
        opt_sol.obj = F[nbjobs][0];

        for (int i =  H_min; i < H_max + 1; i++) {
            if (F[nbjobs][i] < opt_sol.obj) {
                opt_sol.C_max = i;
                opt_sol.obj = F[nbjobs][i];
            }
        }

        t_min = opt_sol.C_max;

        /** Construct the solution */
        for (int i = nbjobs; i >= 1; --i) {
            if (A[i][t_min] && r[i - 1] + p[i - 1] <= t_min && t_min <= d[i - 1]) {
                opt_sol.jobs.push_back(i - 1);
                opt_sol.cost += w[i - 1] * t_min;
                t_min -= p[i - 1];
            }
        }

        /** Free the memory */
        for (int i = 0; i < nbjobs + 1; ++i) {
            delete[] A[i];
            delete[] F[i];
        }

        delete[] A;
        delete[] F;
        return opt_sol;
    }

    class Optimal_Solution<double> dynamic_programming_ahv_farkas(double *pi) {
        Optimal_Solution<double> opt_sol;
        opt_sol.cost = 0;
        double **F;
        bool **A;
        int t_min = H_min;
        F = new double* [nbjobs + 1];
        A = new bool* [nbjobs + 1];

        for (int i = 0; i < nbjobs + 1; i++) {
            F[i] = new double [H_max + 1];
            A[i] = new bool [H_max + 1];
        }

        /** Initialisation */
        F[0][0] = pi[nbjobs];
        A[0][0] = false;

        for (int t = 1; t < H_max + 1; t++) {
            F[0][t] = DBL_MAX / 2;
            A[0][t] = false;
        }

        for (int i = 1; i < nbjobs + 1; i++) {
            for (int t = 0; t < H_max + 1; t++) {
                F[i][t] = DBL_MAX / 2;
                A[i][t] = false;
            }
        }

        /** Recursion */
        for (int i = 1; i < nbjobs + 1; i++) {
            int j = i - 1;

            for (int t = 0; t < H_max; t++) {
                if (t >= r[j] + p[j] && t <= d[j]) {
                    if (F[j][t - p[j]]  + pi[j] < F[j][t]) {
                        F[i][t] = F[j][t - p[j]] + pi[j];
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
        opt_sol.obj = F[nbjobs][0];

        for (int i =  H_min; i < H_max + 1; i++) {
            if (F[nbjobs][i] < opt_sol.obj) {
                opt_sol.C_max = i;
                opt_sol.obj = F[nbjobs][i];
            }
        }

        t_min = opt_sol.C_max;

        /** Construct the solution */
        for (int i = nbjobs; i >= 1; --i) {
            if (A[i][t_min] && r[i - 1] + p[i - 1] <= t_min && t_min <= d[i - 1]) {
                opt_sol.jobs.push_back(i - 1);
                opt_sol.cost += w[i - 1] * t_min;
                t_min -= p[i - 1];
            }
        }

        /** Free the memory */
        for (int i = 0; i < nbjobs + 1; ++i) {
            delete[] A[i];
            delete[] F[i];
        }

        delete[] A;
        delete[] F;
        return opt_sol;
    }



    class Optimal_Solution<double> solve_duration_bdd_double(double *pi) {
        return dd->evaluate_reverse(DurationBDDdouble(pi, p, w, nbjobs));
    }

    class Optimal_Solution<double> solve_duration_zdd_double(double *pi) {
        return zdd->evaluate_forward_DP(DurationZDDdouble(pi, p, w, nbjobs, H_max));
    }

    class Optimal_Solution<double> solve_weight_bdd_double(double *pi) {
        return dd->evaluate_weight(WeightBDDdouble(pi, p, w, r, d, nbjobs), dd_table);
    }

    class Optimal_Solution<double> solve_weight_zdd_double(double *pi) {
        return zdd->evaluate_weight(WeightZDDdouble(pi, p, w, r, d, nbjobs, H_min,
                                    H_max), zdd_table);
    }

    class Optimal_Solution<double> solve_farkas_double(double *pi) {
        return zdd->evaluate_weight(FarkasZDDdouble(pi, p, w, r, d, nbjobs, H_min,
                                    H_max), farkas_table);
    }


    void addDuetimeConstraint() {
    }

    void addReleasetimeConstraint() {
    }

    PricerSolver &operator=(PricerSolver const &other) {
        if (this != &other) {
            dd = new tdzdd::DdStructure<2>;
            zdd = new tdzdd::DdStructure<2>;
            *dd = *(other.dd);
            *zdd = *(other.zdd);
            H_max = other.H_max;
            H_min = other.H_min;
            r = other.r;
            d = other.d;
            w = other.w;
            p = other.p;
            nbjobs = other.nbjobs;
        }

        return *this;
    }

    void set_release_due_time(int *releasetime, int *duetime) {
        r = releasetime;
        d = duetime;
    }

    void iterate_zdd() {
        tdzdd::DdStructure<2>::const_iterator it = zdd->begin();

        for (; it != zdd->end(); ++it) {
            std::vector<int>::const_iterator i = (*it).begin();

            for (; i != (*it).end(); i++) {
                std::cout << nbjobs - *i << " ";
            }

            std::cout << std::endl;
        }
    }


    void iterate_dd() {
        tdzdd::DdStructure<2>::const_iterator it = dd->begin();

        for (; it != dd->end(); ++it) {
            std::vector<int>::const_iterator i = (*it).begin();

            for (; i != (*it).end(); i++) {
                std::cout << nbjobs - *i << " ";
            }

            std::cout << std::endl;
        }
    }


};


