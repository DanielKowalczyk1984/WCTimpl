/*
 * TdZdd: a Top-down/Breadth-first Decision Diagram Manipulation Framework
 * by Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
 * Copyright (c) 2014 ERATO MINATO Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#pragma once


#include "tdzdd/DdEval.hpp"
#include <cfloat>
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
    }
};




namespace tdzdd {
template<typename E>
class MaxReducedCostBaseDbl: public DdEval<E, PricerInfo<double> > {
    double *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxReducedCostBaseDbl(double *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    }

    void evalTerminal( PricerInfo<double>& n, bool one) {
        n.obj = one ? pi[nbjobs] : DBL_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }

    void evalNode( PricerInfo<double> &n, int i, DdValues<PricerInfo<double>, 2> const & values) const {
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
class MaxFarkasPricingBaseDbl: public DdEval<E, PricerInfo<double> > {
    double *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxFarkasPricingBaseDbl(double *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    }

    void evalTerminal( PricerInfo<double>& n, bool one) {
        n.obj = one ? pi[nbjobs] : DBL_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }

    void evalNode( PricerInfo<double> &n, int i, DdValues<PricerInfo<double>, 2> const & values) const {
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
class MaxReducedCostBaseInt: public DdEval<E, PricerInfo<int> > {
    int *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxReducedCostBaseInt(int *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    }

    void evalTerminal( PricerInfo<int>& n, bool one) {
        n.obj = one ? pi[nbjobs] : INT_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }

    void evalNode( PricerInfo<int> &n, int i, DdValues<PricerInfo<int>, 2> const & values) const {
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
class MaxFarkasPricingBaseInt: public DdEval<E, PricerInfo<int> > {
    int *pi;
    int *p;
    int *w;
    int nbjobs;

public:
    MaxFarkasPricingBaseInt(int *_pi, int *_p, int *_w, int _nbjobs)
        : pi(_pi), p(_p), w(_w), nbjobs(_nbjobs) {

    }

    void evalTerminal( PricerInfo<int>& n, bool one) {
        n.obj = one ? pi[nbjobs] : INT_MIN;
        n.sum_w = 0;
        n.cost = 0;
        n.jobs.resize(0);

    }

    void evalNode( PricerInfo<int> &n, int i, DdValues<PricerInfo<int>, 2> const & values) const {
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

}




struct MaxReducedCostDbl: public tdzdd::MaxReducedCostBaseDbl<MaxReducedCostDbl> {
    MaxReducedCostDbl(double *_pi, int *_p, int *_w, int _nbjobs): tdzdd::MaxReducedCostBaseDbl<MaxReducedCostDbl>(_pi, _p, _w, _nbjobs) {

    }
};

struct MaxReducedCostInt: public tdzdd::MaxReducedCostBaseInt<MaxReducedCostInt> {
    MaxReducedCostInt(int *_pi, int *_p, int *_w, int _nbjobs): tdzdd::MaxReducedCostBaseInt<MaxReducedCostInt>(_pi, _p, _w, _nbjobs) {

    }
};

struct MaxFarkasPricingDbl: public tdzdd::MaxFarkasPricingBaseDbl<MaxFarkasPricingDbl> {
    MaxFarkasPricingDbl(double *_pi, int *_p, int *_w, int _nbjobs): tdzdd::MaxFarkasPricingBaseDbl<MaxFarkasPricingDbl>(_pi, _p, _w, _nbjobs) {

    }
};

struct MaxFarkasPricingInt: public tdzdd::MaxFarkasPricingBaseInt<MaxFarkasPricingInt> {
    MaxFarkasPricingInt(int *_pi, int *_p, int *_w, int _nbjobs): tdzdd::MaxFarkasPricingBaseInt<MaxFarkasPricingInt>(_pi, _p, _w, _nbjobs) {

    }
};

