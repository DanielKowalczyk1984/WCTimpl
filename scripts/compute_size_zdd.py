#!/usr/bin/env python3.5
import sys      # imports the sys module
import re
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def diff(x):
    a = x['global_upper_bound']
    b = x['first_lower_bound']

    return a - b


def diff_global(x):
    a = x['global_upper_bound']
    b = x['global_lower_bound']

    return a - b


class DataParser(argparse.ArgumentParser): 

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
                          default=None)
        self.add_argument('-b', '--branch', type=str,
                          help='store the branch method that you want' +
                          ' to calculate',
                          default='CONFLICT')
        self.add_argument('-l', '--list', nargs='+', type=str,
                          help='list of all measurements',
                          default=["size"])
        self.add_argument('-H', '--header', nargs='+', type=str,
                          help='list of all measurements',
                          default=["size"])
        self.add_argument('-d', '--directory', type=str,
                          help='Give the working directory',
                          default='.')
        self.add_argument('-D', '--depth', type=int, nargs='+',
                          help='Give the working directory',
                          default=[0])


class Summary(object):
    """docstring for Summary"""

    def __init__(self, f, m, n, cl, br, lst, depth):
        super(Summary, self).__init__()
        self.df = pd.read_csv(f, sep=";")
        df = self.df[self.df.depth == 0]
        dict_initial = {}
        for idx, row in df.iterrows():
            dict_initial[row.NameInstance] = row['size']

        def initial(a):
            name = a['NameInstance']
            size = a['size']

            return (size - dict_initial[name]) / dict_initial[name]

        self.df['rel'] = self.df.apply(initial, axis=1)
        self.m = m
        self.n = n
        self.br = br
        self.lst = lst
        self.depth = depth
        self.overall = {}
        self.meanlst = {}
        self.geomlst = {}
        self.count = {}
        self.maxlst = {}
        self.overall["opt"] = 0
        self.overall["optroot"] = 0
        for it in self.lst:
            self.overall[it] = []

        self.nb = cl
        if self.nb < 7:
            self.cl = 'instance%d' % self.nb
        else:
            self.cl = 'instance'

    def compute(self):
        df_comp = self.df[self.df['NameInstance'].str.contains(self.cl)]
        for it in self.depth:
            tmp_df = df_comp[df_comp.depth == it]
            self.meanlst[it] = tmp_df['size'].mean()
            self.maxlst[it] = tmp_df['size'].max()
            self.count[it] = tmp_df['NameInstance'].count()

        df_filter = df_comp[df_comp.depth%2 == 0]
        df_filter.boxplot(by='depth', column='rel', whis='range')
        plt.show()

    def print_to_screen(self):
        output = str()
        if self.cl == "instance":
            output += "overall," + str(self.m) + "," + str(self.n)
        else:
            output += str(self.nb) + "," + \
                str(self.m) + "," + str(self.n)

        for it in self.depth:
            if self.count[it] > 0:
                output += ",%f,%d" % (self.meanlst[it], self.maxlst[it])
            else:
                output += ",,"

        print(output)


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
        print(OSError.errno)

    header = "class,m,n"
    for it in args.depth:
        header += ",avg%d,max%d" % (it, it)

    print(header)

    if args.class_inst is None:
        for cl in range(1, 7):
            for n in args.jobs:
                for m in args.machines:
                    try:
                        fname = '%s_%d_%d.csv' % (args.branch, m, n)
                        inst = Summary(fname, m, n, cl,
                                       args.branch, args.list, args.depth)

                    except IOError:
                        print("Could not read file:", fname)
                    else:
                        inst.compute()
                        inst.print_to_screen()
    else:
        for n in args.jobs:
            for m in args.machines:
                try:
                    fname = '%s_%d_%d.csv' % (args.branch, m, n)
                    inst = Summary(fname, m, n, args.class_inst,
                                   args.branch, args.list, args.depth)

                except IOError:
                    print("Could not read file:", fname)
                else:
                    inst.compute()
                    inst.print_to_screen()

if __name__ == '__main__':
    main()
