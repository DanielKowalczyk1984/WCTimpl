#include "tdzdd/DdSpec.hpp"
#include <vector>
#include <boost/dynamic_bitset.hpp>

class conflict_state {
  public:
    boost::dynamic_bitset<> add;
    boost::dynamic_bitset<> remove;

    conflict_state() {
    };

    ~conflict_state() {};

};


class PricerSpec: public tdzdd::DdSpec<PricerSpec, int, 2> {
    int nbjobs;
    int *p;
    int *r;
    int *d;
    int Hmin;
    int Hmax;

  public:
    int *sum_p;
    int *min_p;
    PricerSpec(int *_p, int *_r, int *_d, int njobs, int Hmin, int Hmax): p(_p),
        r(_r), d(_d), Hmin(Hmin), Hmax(Hmax) {
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

    ~PricerSpec() {
    }

    int getRoot(int &state) const {
        state = 0;
        return nbjobs;
    }

    int getChild(int &state, int level, int value) const {
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
                if (state + sum_p[_j] < Hmin) {
                    return 0;
                }
  
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

        // if (_j == nbjobs && state >= Hmin && state <= Hmax) {
        //     return -1;
        // } else if (_j == nbjobs) {
        //     return 0;
        // }

        assert(_j < nbjobs);
        return nbjobs - _j;
    }

  private:
    int min_job(int j, int state) const {
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

class ConflictConstraints: public
    tdzdd::DdSpec<ConflictConstraints, conflict_state, 2> {
    int nbjobs;
    std::vector<boost::dynamic_bitset<> > differsets;
    std::vector<boost::dynamic_bitset<> > samesets;

    bool takeable(int job, conflict_state &state) {
        if (state.remove[job]) {
            return false;
        }

        return true;
    }

    bool leaveable(int job, conflict_state &state) {
        if (state.add[job]) {
            return false;
        }

        return true;
    }

  public:
    ConflictConstraints(int _nbjobs, int *elist_same, int ecount_same,
                        int *elist_differ, int ecount_differ): nbjobs(_nbjobs) {
        differsets.resize(_nbjobs);
        samesets.resize(_nbjobs);

        for (int i = 0; i < _nbjobs; i++) {
            differsets[i].resize(_nbjobs);
            samesets[i].resize(_nbjobs);
        }

        for (int i = 0; i < ecount_same; ++i) {
            samesets[elist_same[2 * i]][elist_same[2 * i + 1]] = 1;
        }

        for (int i = 0; i < ecount_differ; ++i) {
            differsets[elist_differ[2 * i]][elist_differ[2 * i + 1]] = 1;
        }
    };

    ~ConflictConstraints() {
    }

    int getRoot(conflict_state &state) const {
        state.add.resize(nbjobs);
        state.remove.resize(nbjobs);
        return nbjobs;
    }

    int getChild(conflict_state &state, int level, int take) {
        int job = nbjobs - level;
        int _j;
        assert(0 <= job && job <= nbjobs - 1);

        if (samesets[job].intersects(differsets[job])) {
            return 0;
        }

        if (level - 1 == 0 && take) {
            return (!state.remove[job]) ? -1 : 0;
        } else if (level - 1 == 0) {
            return (!state.add[job]) ? -1 : 0;
        }

        if (take) {
            if (!takeable(job, state)) {
                return 0;
            }

            state.add |= samesets[job];
            state.remove |= differsets[job];
        } else {
            if (!leaveable(job, state)) {
                return 0;
            }

            state.remove |= samesets[job];
        }

        _j = min_job(job, state);

        if (_j == nbjobs && take) {
            return (!state.remove[job]) ? -1 : 0;
        } else if (_j == nbjobs) {
            return (!state.add[job]) ? -1 : 0;
        }

        assert(_j < nbjobs);
        return nbjobs - _j;
    }

    bool equalTo(conflict_state const &state1, conflict_state const &state2) const {
        if (state2.add != state1.add) {
            return false;
        }

        if (state2.remove != state1.remove) {
            return false;
        }

        return true;
    }

    size_t hashCode(conflict_state const &state) const {
        size_t val = 0;
        size_t it = state.add.find_first();

        while (it != boost::dynamic_bitset<>::npos) {
            val += 1213657 * static_cast<size_t>(it);
            it = state.add.find_next(it);
        }

        it = state.remove.find_first();

        while (it != boost::dynamic_bitset<>::npos) {
            val += 487239 * static_cast<size_t>(it);
            it = state.remove.find_next(it);
        }

        return val;
    }

    int min_job(int j, conflict_state &state) const {
        int i, val = nbjobs;

        for (i = j + 1; i < nbjobs; ++i) {
            if (!state.remove[i]) {
                val = i;
                break;
            }
        }

        return val;
    }
};

