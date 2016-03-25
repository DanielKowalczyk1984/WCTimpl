#include "PricerConstruct.hpp"
#include "PricerEvaluate.hpp"
#include "tdzdd/DdStructure.hpp"
#include <iostream>
#include <vector>

class PricerSolver
{
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
        tdzdd::DataTable<PricerWeightZDD<double>> zdd_table;
        tdzdd::DataTable<PricerWeightBDD<double>> dd_table;
        tdzdd::DataTable<PricerFarkasZDD<double>> farkas_table;

        PricerSolver(int *_p, int *_w,  int *_r, int *_d, int njobs, int Hmin, int Hmax, bool print = false, bool reduce = false): p(_p), w(_w), r(_r), d(_d), nbjobs(njobs), H_min(Hmin), H_max(Hmax)
        {
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
                //if (print) {
                //}
            }

            delete [] ps.sum_p;
            delete [] ps.min_p;
        };

        void create_dot_zdd(const char *name)
        {
            std::ofstream file;
            file.open(name);
            zdd.dumpDot(file);
            file.close();
        }

        void init_table_farkas()
        {
            tdzdd::NodeTableHandler<2> &node_handler_farkas = zdd.getDiagram();
            farkas_table.init(nbjobs + 1);
            size_t const m = node_handler_farkas.privateEntity()[0].size();
            farkas_table[0].resize(m);

            for (size_t i = 0; i < m; ++i) {
                farkas_table[0][i].init_terminal_node(i);
            }

            for (int i = 1; i <= nbjobs; ++i) {
                tdzdd::MyVector<tdzdd::Node<2> > const &node = node_handler_farkas.privateEntity()[i];
                size_t const mm = node.size();
                farkas_table[i].resize(mm);

                for (size_t j = 0; j < mm; ++j) {
                    farkas_table[i][j].init_node();
                }
            }
        }

        void init_bdd_table()
        {
            tdzdd::NodeTableHandler<2> &handler = dd.getDiagram();
            tdzdd::NodeId &root = dd.root();
            dd_table.init(root.row() + 1);

            /** init table */
            for (int i = root.row(); i >= 0 ; i--) {
                tdzdd::MyVector<tdzdd::Node<2>> const &node = handler.privateEntity()[i];
                size_t const m = node.size();
                dd_table[i].resize(m);
            }

            /** init root */
            dd_table[root.row()][root.col()].init_node(0);

            for (size_t i = nbjobs; i > 1 ; i--) {
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

        void init_zdd_table()
        {
            tdzdd::NodeTableHandler<2> &handler = zdd.getDiagram();
            tdzdd::NodeId &root = zdd.root();
            zdd_table.init(root.row() + 1);

            /** init table */
            for (int i = root.row(); i >= 0 ; i--) {
                tdzdd::MyVector<tdzdd::Node<2>> const &node = handler.privateEntity()[i];
                size_t const m = node.size();
                zdd_table[i].resize(m);
            }

            /** init root */
            zdd_table[root.row()][root.col()].add_weight(0);

            for (size_t i = nbjobs; i >= 1 ; i--) {
                size_t const m = zdd_table[i].size();
                int cur_job = nbjobs - i;

                for (size_t j = 0; j < m; j++) {
                    std::vector<int>::iterator it =  zdd_table[i][j].weight.begin();
                    tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
                    tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

                    for (; it != zdd_table[i][j].weight.end(); it++) {
                        zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight(*it);
                        zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight(*it + p[cur_job]);
                    }
                }
            }

            /** init terminal nodes */
            size_t const mm = handler.privateEntity()[0].size();

            for (size_t j = 0; j < mm ; j++) {
                zdd_table[0][j].init_terminal_node(j, H_max);
            }
        }

        class Optimal_Solution<double> dynamic_programming_ahv(double *pi)
        {
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



        class Optimal_Solution<double> solve_duration_bdd_double(double *pi)
        {
                return dd.evaluate_reverse(DurationBDDdouble(pi, p, w, nbjobs));
        }

        class Optimal_Solution<double> solve_duration_zdd_double(double *pi)
        {
                return zdd.evaluate_forward_DP(DurationZDDdouble(pi, p, w, nbjobs, H_max));
        }

        class Optimal_Solution<double> solve_weight_bdd_double(double *pi)
        {
                return dd.evaluate_weight(WeightBDDdouble(pi, p, w, r, d, nbjobs), dd_table);
        }

        class Optimal_Solution<double> solve_weight_zdd_double(double *pi)
        {
                return zdd.evaluate_weight(WeightZDDdouble(pi, p, w, r, d, nbjobs, H_min, H_max), zdd_table);
        }

        class Optimal_Solution<double> solve_farkas_double(double *pi)
        {
                return zdd.evaluate_weight(FarkasZDDdouble(pi, p, w, r, d, nbjobs, H_min, H_max), farkas_table);
        }

        void addConflictConstraints(int *elist_same, int ecount_same, int *elist_differ, int ecount_differ)
        {
            tmp = zdd;
            ConflictConstraints conflict(nbjobs, elist_same, ecount_same, elist_differ, ecount_differ);
            tmp.zddSubset(conflict);
        }

        void addDuetimeConstraint()
        {
        }

        void addReleasetimeConstraint()
        {
        }

        void iterate_zdd()
        {
            tdzdd::DdStructure<2>::const_iterator it = zdd.begin();

            for (; it != zdd.end(); ++it) {
                std::vector<int>::const_iterator i = (*it).begin();

                for (; i != (*it).end(); i++) {
                    std::cout << nbjobs - *i << " ";
                }

                std::cout << std::endl;
            }
        }


        void iterate_dd(){
            tdzdd::DdStructure<2>::const_iterator it = dd.begin();

            for (; it != dd.end(); ++it) {
                std::vector<int>::const_iterator i = (*it).begin();

                for (; i != (*it).end(); i++) {
                    std::cout << nbjobs - *i << " ";
                }

                std::cout << std::endl;
            }
        }

        ~PricerSolver()
        {
        };


};


