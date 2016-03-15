#include "tdzdd/DdSpec.hpp"
#include <vector>


class PricerSpec: public tdzdd::DdSpec<PricerSpec, int, 2>
{
        int nbjobs;
        int *p;
        int *r;
        int *d;
        int Hmin;
        int Hmax;

    public:
        int *sum_p;
        int *min_p;
        PricerSpec(int *_p, int *_r, int *_d, int njobs, int Hmin, int Hmax): p(_p), r(_r), d(_d), Hmin(Hmin), Hmax(Hmax)
        {
            nbjobs = njobs;
            sum_p = new int [nbjobs];
            min_p = new int [nbjobs];
            int end = nbjobs - 1;
            sum_p[end] = p[end];
            min_p[end] = p[end];

            for (int  i = end - 1; i >= 0; i--) {
                sum_p[i] = sum_p[ i + 1] + p[i];

                if (p[i] < min_p[i + 1]) {
                    min_p[i] = p[i];
                } else {
                    min_p[i] = min_p[i + 1];
                }
            }
        }

        ~PricerSpec()
        {
        }

        int getRoot(int &state) const
        {
            state = 0;
            return nbjobs;
        }

        int getChild(int &state, int level, int value) const
        {
            int job = nbjobs - level;
            int _j;
            assert(0 <= job && job <= nbjobs - 1);
            int temp_p = p[job];

            if (level - 1 == 0 && value) {
                return (state + p[job] >= Hmin && state + p[job] <= Hmax) ? -1 : 0;
            } else if (level - 1 == 0) {
                return (state  >= Hmin && state <= Hmax) ? -1 : 0;
            }

            if (value) {
                int sum = state + temp_p;
                _j = min_job(job, sum);

                if (_j < nbjobs) {
                    if ((sum >= Hmin && sum <= Hmax) && (sum + min_p[_j] > Hmax)) {
                        return -1;
                    }
                } else {
                    if ((sum >= Hmin && sum <= Hmax)) {
                        return -1;
                    }

                    return 0;
                }

                state = sum;
            } else {
                _j = min_job(job, state);

                if (_j < nbjobs) {
                    if (state + sum_p[_j] < Hmin) {
                        return 0;
                    }

                    if ((state >= Hmin && state  <= Hmax) && (state + min_p[_j] > Hmax)) {
                        return -1;
                    }
                } else {
                    if (state >= Hmin && state <= Hmax) {
                        return -1;
                    }

                    return 0;
                }
            }

            if (_j == nbjobs && state >= Hmin && state <= Hmax) {
                return -1;
            } else if (_j == nbjobs) {
                return 0;
            }

            assert(_j < nbjobs);
            return nbjobs - _j;
        }

    private:
        int min_job(int j, int state) const
        {
            int i, val = nbjobs;

            for (i = j + 1; i < nbjobs; ++i) {
                if (state >= r[i] && state <= d[i] - p[i]) {
                    val = i;
                    break;
                }
            }

            return val;
        }

};

class Restriction: public tdzdd::DdSpec<Restriction, int, 2>
{
        int nbjobs;
        int n;
        int *jobs;
        int *p;
        int *r;
        int *d;
    public:
        Restriction()
        {
        };

        int getRoot(int   &state) const
        {
            state = 0;
            return nbjobs;
        }
};

class Duedate: public tdzdd::DdSpec<Duedate, int, 2>
{
        int nbjobs;
        int n;
        int *jobs;
        int *p;
        int *r;
        int *d;
    public:
        int *sum_p;
        int *min_p;
        Duedate(int *_p, int *_r, int *_d, int njobs): p(_p), r(_r), d(_d)
        {
            nbjobs = njobs;
            sum_p = new int [nbjobs];
            min_p = new int [nbjobs];
            int end = nbjobs - 1;
            sum_p[end] = p[end];
            min_p[end] = p[end];

            for (int  i = end - 1; i >= 0; i--) {
                sum_p[i] = sum_p[ i + 1] + p[i];

                if (p[i] < min_p[i + 1]) {
                    min_p[i] = p[i];
                } else {
                    min_p[i] = min_p[i + 1];
                }
            }
        }

        int getRoot(int &state) const
        {
            state = 0;
            return nbjobs;
        }

        int getChild(int &state, int level, int value)
        {
            return 0;
        }

    private:
        int min_job(int j, int state) const
        {
            int i, val = nbjobs;

            for (i = j + 1; i < nbjobs; ++i) {
                if (state >= r[i] && state <= d[i] - p[i]) {
                    val = i;
                    break;
                }
            }

            return val;
        }
};

class Releasetime:  public tdzdd::DdSpec<Releasetime, int, 2>
{
        int nbjobs;
        int n;
        int *jobs;
        int *p;
        int *r;
        int *d;
    public:
        Releasetime(int _nbjobs, int _n): nbjobs(_nbjobs), n(_n)
        {
        };

        int getRoot(int &state) const
        {
            state = 0;
            return nbjobs;
        }

        int getChild(int &state, int level, int value)
        {
            return 0;
        }
    private:
        int min_job(int j, int state) const
        {
            int i, val = nbjobs;

            for (i = j + 1; i < nbjobs; ++i) {
                if (state >= r[i] && state <= d[i] - p[i]) {
                    val = i;
                    break;
                }
            }

            return val;
        }
};

