/*
 * Top-Down ZDD Construction Library for Frontier-Based Search
 * by Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
 * Copyright (c) 2012 Japan Science and Technology Agency
 * $Id: ToZBDD.hpp 426 2013-02-26 06:50:04Z iwashita $
 */

#pragma once

#ifndef B_64
#if __SIZEOF_POINTER__ == 8
#define B_64
#endif
#endif
#include "../../SAPPOROBDD/ZBDD.h"

#include <vector>

#include "../DdEval.hpp"

/**
 * Exporter to ZBDD.
 * TdZdd nodes at level @a i are converted to
 * ZBDD nodes at level @a i + @p offset.
 * When the ZBDD variables are not enough, they are
 * created automatically by BDD_NewVar().
 */
template<typename E>
class ToZBDDBase: public tdzdd::DdEval<E, ZBDD>
{
        int const offset;

    public:
        ToZBDDBase(int offset = 0)
            : offset(offset)
        {
        }

        void initialize(int topLevel) const
        {
            while (BDD_VarUsed() < topLevel + offset) {
                BDD_NewVar();
            }
        }

        void evalTerminal(ZBDD &f, bool one) const
        {
            f = ZBDD(one ? 1 : 0);
        }

        void evalNode(ZBDD &f, int level, tdzdd::DdValues<ZBDD, 2>   &values) const
        {
            ZBDD f0 = values.get(0);
            ZBDD f1 = values.get(1);

            if (level + offset <= 0) {
                f = f0;
            } else {
                f = f0 + f1.Change(BDD_VarOfLev(level + offset));
            }
        }
};

struct ToZDBB: public ToZBDDBase<ToZDBB> {
    ToZDBB(int offset = 0): ToZBDDBase<ToZDBB>(offset)
    {
    }
};
