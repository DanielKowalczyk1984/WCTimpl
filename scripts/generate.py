#!/usr/bin/env python2.7

from random import randint
import argparse
import os


class OptionParser(argparse.ArgumentParser):
    """OptionParser"""

    def __init__(self):
        super(OptionParser, self).__init__()
        self.add_argument('-c', '--class_inst', type=int, nargs='+',
                          help='class instance',
                          default=[1, 2, 3, 4, 5, 6])
        self.add_argument('-n', '--number_jobs', type=int, nargs='+',
                          help='number of jobs', default=[20, 50, 100, 150])
        self.add_argument('-m', '--number_machines', type=int, nargs='+',
                          help='number of machines', default=[3, 5, 8, 10, 12])


class ValidationError(Exception):
    """Vaildation Error"""

    def __init__(self, msg, err):
        super(ValidationError, self).__init__()
        self.msg = msg
        self.err = err


def give_weight_duration(n):
    if n > 6:
        raise ValidationError('Error occured number to large', n)
    elif n == 1:
        return randint(10, 100), randint(1, 10)
    elif n == 2:
        return randint(1, 100), randint(1, 100)
    elif n == 3:
        return randint(10, 20), randint(10, 20)
    elif n == 4:
        return randint(90, 100), randint(90, 100)
    elif n == 5:
        p = randint(90, 100)
        w = randint(p - 5, p + 5)
        return w, p
    elif n == 6:
        p = randint(10, 100)
        w = randint(p - 5, p + 5)
        return w, p


def main():
    try:
        parser = OptionParser()
        args = parser.parse_args()
    except argparse.ArgumentError, e:
        raise e

    for m in args.number_machines:
        for n in args.number_jobs:
            directory = './%d_%d' % (n, m)
            os.mkdir(directory)
            for c in args.class_inst:
                for i in xrange(0, 10):
                    try:
                        f = open(directory + '/instance%d_%d_%d_1%d.txt' %
                                 (c, n, m, i), 'w')
                    except IOError, (ErrorNumber, ErrorMessage):
                        print ErrorMessage
                        break

                    f.write('p %d\n' % n)
                    p = int()
                    w = int()
                    for it in xrange(0, n):
                        try:
                            w, p = give_weight_duration(c)
                        except ValidationError as e:
                            print e.msg
                            break
                        f.write('n %d %d\n' % (w, p))

if __name__ == '__main__':
    main()
