#!/usr/bin/python

import gzip
import sys


def parse(path1, path2):
    contig  = None
    handle1 = None
    handle2 = None
    if path1[:-3] == '.gz':
        handle1 = gzip.open(path1)
    else:
        handle1 = open(path1)
    if path2[:-3] == '.gz':
        handle2 = gzip.open(path2)
    else:
        handle2 = open(path2)

    line1 = handle1.readline().strip()
    line2 = handle2.readline().strip()
    while line1 and line2:
        if line1[0] == '>':
            if line1 != line2:
                raise "Incompatible formats"
            contig = line1[1:].split('|')[-1]
            line1  = handle1.readline().strip()
            line2  = handle2.readline().strip()
            continue
        elif not contig:
            line1  = handle1.readline().strip()
            line2  = handle2.readline().strip()
            continue
        fields1   = line1.split()
        position1 = fields1[0]
        n_meth1   = fields1[1]
        n_unmeth1 = fields1[2]
        fields2   = line2.split()
        position2 = fields2[0]
        n_meth2   = fields2[1]
        n_unmeth2 = fields2[2]
        if position1 != position2:
            raise "Incompatible formats"
        t = (contig, position1, n_meth1, n_unmeth1, n_meth2, n_unmeth2)
        print '\t'.join(t)
        line1 = handle1.readline().strip()
        line2 = handle2.readline().strip()

    handle1.close()
    handle2.close()

def main():
    from optparse import OptionParser
    usage  = 'Usage: %prog [options] FILE1 FILE2'
    parser = OptionParser(usage)
#    parser.add_option('-m', '--min_meth',
#            dest='min_meth',
#            type = 'int',
#            default=0,
#            help='Minimum number of methylated Cs')
#    parser.add_option('-r', '--min_ratio',
#            dest='min_ratio',
#            type = 'float',
#            default=0.0,
#            help='Minimum ratio of methylated Cs over total Cs')

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error('You must pass two file as argument')

    parse(args[0], args[1])


if __name__ == '__main__':
    main()
