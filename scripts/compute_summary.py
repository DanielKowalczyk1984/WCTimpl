#!/usr/bin/env python2
import csv     # imports the csv module
import sys      # imports the sys module
import re
import os
import scipy.stats
import numpy as np
import argparse


class DataParser(argparse.ArgumentParser):
    """DataParser"""

    def __init__(self):
        super(DataParser, self).__init__()
        self.add_argument('-n', '--jobs', nargs='+', type=int,
                          help='store a list of number of jobs ' +
                          'that you want to calculate',
                          default=[20, 50, 100, 150])
        self.add_argument('-m', '--machines', nargs='+', type=int,
                          help='store a list of  number of jobs that' +
                          ' you want to calculate',
                          default=[3, 5, 8, 10, 12])
        self.add_argument('-c', '--class_inst', type=int,
                          help='store the class that you want to calculate',
                          default=0)
        self.add_argument('-b', '--branch', type=str,
                          help='store the branch method that you want' +
                          ' to calculate',
                          default='CONFLICT')
        self.add_argument('-l', '--list', nargs='+', type=str,
                          help='list of all measurements',
                          default=["tot_real_time", "tot_cputime", "tot_lb",
                                   "tot_lb_root", "tot_lb_lp",
                                   "tot_branch_and_bound", "tot_build_dd",
                                   "tot_pricing", "nb_explored_nodes",
                                   "nb_generated_col"])
        self.add_argument('-H', '--header', nargs='+', type=str,
                          help='list of all measurements',
                          default=["real sec", "cpu", "cpu LB",
                                   "cpu LB root", "cpu LB LP",
                                   "cpu B&B", "cpu build",
                                   "cpu pricing", "#node",
                                   "#col"])
        self.add_argument('-d', '--directory', type=str,
                          help='Give the working directory',
                          default='.')


class Summary(object):
    """docstring for Summary"""

    def __init__(self, f, m, n, cl, br, lst):
        super(Summary, self).__init__()
        self.reader = csv.DictReader(f)
        self.m = m
        self.n = n
        self.br = br
        self.overall = {}
        self.meanlst = {}
        self.geomlst = {}
        self.maxlst = {}
        self.overall["#opt"] = 0
        self.overall["#opt root"] = 0
        self.lst = lst
        for it in self.lst:
            self.overall[it] = []
        if cl == 0:
            self.cl = "all"
            self.p = re.compile('instance')
        else:
            self.cl = str(cl)
            self.p = re.compile('instance' + self.cl)

    def compute(self):
        for row in self.reader:
            if re.match(self.p, row['NameInstance']) and row['status'] == "4":
                self.overall["#opt"] += 1
                if row['solved_at_root'] == "1":
                    self.overall["#opt root"] += 1
                for it in self.lst:
                    self.overall[it].append(float(row[it]))

        if self.overall["#opt"] > 0:
            for it in self.lst:
                self.meanlst[it] = np.mean(self.overall[it])
                self.maxlst[it] = np.amax(self.overall[it])
                # self.geomlst[it] = scipy.stats.gmean(self.overall[it],axis=None)

    def print_to_screen(self):
        output = str()
        if self.cl == "all":
            output += "overall," + str(self.m) + "," + str(self.n)
        else:
            output += self.cl + "," + \
                str(self.m) + "," + str(self.n)

        output += ",%d" % self.overall["#opt"]
        output += ",%d" % self.overall["#opt root"]
        if self.overall["#opt"] > 0:
            for it in self.lst:
                output += ",%f" % self.meanlst[it]
        else:
            for it in self.lst:
                output += ",-"

        print output


def main():
    try:
        parser = DataParser()
    except argparse.ArgumentError:
        sys.exit(2)
    else:
        args = parser.parse_args()

    try:
        directory = args.directory
        os.chdir(directory)
    except OSError:
        print OSError.errno

    header = "class,m,n,#opt,#opt root"
    for it in args.header:
        header += ",%s" % it
    print header

    for m in args.machines:
        for n in args.jobs:
            try:
                fname = 'WCT_%s_%d_%d.csv' % (args.branch, m, n)
                f = open(fname, 'rb+')
            except IOError:
                print "Could not read file:", fname
            else:
                # initialize summary class
                summary = Summary(f, m, n, args.class_inst,
                                  args.branch, args.list)

                # compute everything
                summary.compute()

                # print summary
                summary.print_to_screen()

                f.close()

if __name__ == '__main__':
    main()
