#!/usr/bin/python

import gzip
import sys


def parse(path, min_meth, min_ratio, origin):
    handle = None
    contig = None
    if path == '-':
        handle = sys.stdin
    elif path[:-3] == '.gz':
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            contig = line[1:].split('|')[-1]
            continue
        elif not contig:
            continue
        fields   = line.split()
        position = int(fields[0])
        n_meth   = int(fields[1])
        n_unmeth = int(fields[2])
        tot      = n_meth + n_unmeth
        if n_meth < min_meth:
            continue
        if tot > 0:
            ratio = float(n_meth) / tot
            if ratio >= min_ratio:
                t = (contig,
                     'meth_ratio',
                     origin,
                     str(position + 1),
                     str(position + 2),
                     str(int(ratio * 100)),
                     '.',
                     '.',
                     '%s %s:methratio' % (origin, contig))
                print '\t'.join(t)
    if path != '-':
        handle.close()

def main():
    from optparse import OptionParser
    usage  = 'Usage: %prog [options] FILE'
    parser = OptionParser(usage)
    parser.add_option('-m', '--min_meth',
            dest='min_meth',
            type = 'int',
            default=0,
            help='Minimum number of methylated Cs')
    parser.add_option('-r', '--min_ratio',
            dest='min_ratio',
            type = 'float',
            default=0.0,
            help='Minimum ratio of methylated Cs over total Cs')
    parser.add_option('-o', '--origin',
            dest='origin',
            type = 'string',
            default='meth_ratio',
            help='Name for the "origin" field')

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error('You must pass a file as argument')

    parse(args[0], options.min_meth, options.min_ratio, options.origin)


if __name__ == '__main__':
    main()
