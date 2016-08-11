#!/usr/bin/env python

import getopt
import sys
import glob
import os


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "db:o:O:mr:M:f:w:c:u:s:l:R:S:C:B:W:i:j:z:")
    except getopt.GetoptError as err:
        # print help information and exit:

        print str(err)
        getopt.usage()
        sys.exit(2)
    opt = ''
    for o, a in opts:
        if o == "-d":
            opt += "-d "
        elif o in ("-b"):
            opt += "-b %s " % (a)
        elif o in ("-o"):
            opt += "-o %s " % (a)
        elif o in ("-o"):
            opt += "-o %s " % (a)
        elif o in ("-O"):
            opt += "-O %s " % (a)
        elif o in "-m":
            opt += "-m "
        elif o in ("-r"):
            opt += "-r %s " % (a)
        elif o in ("-M"):
            opt += "-M %s " % (a)
        elif o in ("-f"):
            opt += "-f %s " % (a)
        elif o in ("-R"):
            opt += "-R %s " % (a)
        elif o in ("-w"):
            opt += "-w %s " % (a)
        elif o in ("-z"):
            opt += "-z %s " % (a)
        elif o in ("-u"):
            opt += "-u %s " % (a)
        elif o in ("-s"):
            opt += "-s %s " % (a)
        elif o in ("-l"):
            opt += "-l %s " % (a)
        elif o in ("-S"):
            opt += "-S %s " % (a)
        elif o in ("-C"):
            opt += "-C %s " % (a)
        elif o in ("-B"):
            opt += "-B %s " % (a)
        elif o in ("-i"):
            opt += "-i %s" % (a)
        elif o in ("-j"):
            opt += "-j %s" % (a)
        if o in ("-W"):
            opt += "-W %s " % (a)

    lst = glob.glob("*.txt")
    lst.sort()
    for file in lst:
        a = sys.argv[-1]
        cmd = './wct ' + opt + ' %s %s' % (file, a)
        os.system(cmd)

if __name__ == "__main__":
    main()
