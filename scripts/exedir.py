#!/usr/bin/env python

import sys
import os
import argparse
import re
from shutil import copyfile


def String_number_bound(n):
    def Result(v):
        v = str(v)
        try:
            return re.match("[0-%d]" % n, v).group(0)
        except:
            raise argparse.ArgumentTypeError(
                "String '%s' does not match required format" % (v,))
    return Result


def String_is_number(s):
    try:
        return re.match("[0-9]+", s).group(0)
    except:
        raise argparse.ArgumentTypeError(
            "String '%s' does not match required format" % (s,))


class OptionParser(argparse.ArgumentParser):
    """DataParser"""

    def __init__(self):
        super(OptionParser, self).__init__()
        test1 = String_number_bound(1)
        test2 = String_number_bound(2)
        test3 = String_number_bound(3)
        self.add_argument('-f', '--nb_solutions',
                          type=String_is_number,
                          help="Number of feasible solutions to construct",
                          default="0")
        self.add_argument('-S', '--stabilization', type=test2,
                          default='1')
        self.add_argument('-z', '--solver', type=test2,
                          help='Pricing solver technique',
                          default='1')
        self.add_argument('-c', '--construct', type=test1,
                          help='Construct heuristically solutions:' +
                          ' 0 = yes(default), 1 = no',
                          default='0')
        self.add_argument('-D', '--diving', type=test1,
                          help='Use diving heuristic: 0 = no(default)' +
                          ', 1 = yes',
                          default='0')
        self.add_argument('-t', '--heur_ahv', type=test1,
                          help='Use ahv test: 0 = no(default), 1 = yes',
                          default='0')
        self.add_argument('-p', '--print', type=test1,
                          help='Print csv-files: 0 = no(default), 1 = yes',
                          default='0')
        self.add_argument('-b', '--branching', type=test3,
                          help='Branching strategy: 0 = conflict(default),' +
                          '1 = ahv,2 = conflict cbfs, 3 = ahv cbfs',
                          default='0')
        self.add_argument('-Z', '--strong_branching',
                          type=test1,
                          help='Use strong branching' +
                          ': 0 = use strong branching(default)' +
                          ', 1 = no strong branching',
                          default='1')
        self.add_argument('-r', '--scatter_search', type=test1,
                          help='Scatter search use' +
                          ': 0 = no scatter search(default)' +
                          ', 1 = scatter search',
                          default='0')
        self.add_argument('-C', '--combine_method', type=test1,
                          help='Combine method scatter search' +
                          ': 0 = Pathrelinking, 1 = PM',
                          default='0')
        self.add_argument('-l', '--time_limit', type=String_is_number,
                          help='Cpu time limit',
                          default='3600')
        self.add_argument('-n', '--number_jobs', nargs='+', type=int,
                          help='List of number of jobs that' +
                          ' has to be processed',
                          default=[20, 50, 100])
        self.add_argument('-m', '--number_machines', nargs='+', type=int,
                          help='List of machines that has to be processed',
                          default=[12, 10, 8, 5, 3])
        self.add_argument('-i', '--instances', type=str,
                          help='directory were all the instances are',
                          default='')
        self.add_argument('-d', '--directory', type=str,
                          help='path of the directory the executable file is',
                          default='.')


def main():
    try:
        parser = OptionParser()
        args = parser.parse_args()
    except argparse.ArgumentError:
        sys.exit(2)

    options_dict = vars(args)
    options_keys = options_dict.keys()
    str_opt = str()
    directory = '/home/daniel/WCTcomputation'
    exe_directory = args.directory
    instances = str()

    for opt in options_keys:
        if opt == "nb_solutions":
            str_opt += "-f %s " % (options_dict[opt])
        elif opt == "branching":
            str_opt += "-b %s " % (options_dict[opt])
        elif opt == "stabilization":
            str_opt += "-S %s " % (options_dict[opt])
        elif opt == "solver":
            str_opt += "-z %s " % (options_dict[opt])
        elif opt == "time_limit":
            str_opt += "-l %s " % (options_dict[opt])
        elif opt == "scatter_search":
            str_opt += "-r %s " % (options_dict[opt])
        elif opt == "construct":
            str_opt += "-c %s " % (options_dict[opt])
        elif opt == "strong_branching":
            str_opt += "-Z %s " % (options_dict[opt])
        elif opt == "heur_ahv":
            str_opt += "-t %s " % (options_dict[opt])
        elif opt == "print":
            str_opt += "-p %s " % (options_dict[opt])
        elif opt == "diving":
            str_opt += "-D %s " % (options_dict[opt])
        elif opt == "combine_method":
            str_opt += "-C %s " % (options_dict[opt])
        elif opt == "instances":
            instances = r'instance%s.*\.txt' % (options_dict[opt])
    try:
        os.chdir(exe_directory)
    except OSError:
        print OSError.errno
        sys.exit(status=2)
    else:
        for n in args.number_jobs:
            for m in args.number_machines:
                print instances
                files = [f for f in os.listdir(
                    directory + '/%d_%d' % (n, m)) if re.match(instances, f)]
                files.sort()

                for f in files:
                    fullpath_instance = directory + '/%d_%d/' % (n, m) + f
                    cmd = './wct ' + str_opt + fullpath_instance + ' ' + str(m)
                    os.system(cmd)

                csv_files = [f for f in os.listdir('.') if re.match(
                    r'.*CONFLICT_%d_%d.csv' % (m, n), f)]
                for f in csv_files:
                    copyfile(f, '/home/daniel/Dropbox/' + f)


if __name__ == "__main__":
    main()
