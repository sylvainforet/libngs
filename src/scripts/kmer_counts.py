#!/usr/bin/python

import gzip
import logging
import os.path
import string
import sys

rev_table = string.maketrans('ACGTNacgtn', 'TGCANtgcan')

def rev_comp(s):
    return s.translate(rev_table)[::-1]

def getKmerCounts(k, path):
    logging.debug('Parsing %s' % path)
    handle = sys.stdin
    if path != '-':
        if path[:-3] == '.gz':
            handle = gzip.open(path)
        else:
            handle = open(path)
    kmers = {}
    n = 0
    line1 = handle.readline() # @name
    while line1:
        line2 = handle.readline().strip() #seq
        for i in xrange(len(line2) - k + 1):
            word = line2[i:i + k]
            if not 'N' in word:
                if not word in kmers:
                    kmers[word] = 0
                kmers[word] += 1
                word = rev_comp(word)
                if not word in kmers:
                    kmers[word] = 0
                kmers[word] += 1
        n += 1
        if n % 10**6 == 0:
            logging.debug('Parsed %d sequences' % n)
        handle.readline() # +name
        handle.readline() # qual
        line1 = handle.readline() # @name
    if path != '-':
        handle.close()
    return kmers

def printKmerCountsSorted(kmers, path):
    logging.debug('Writing kmer counts to %s' % path)
    handle = sys.stdout
    if path != '-':
        handle = open(path)
    logging.debug('Getting Keys')
    kmerKeys = kmers.keys()
    logging.debug('Sorting Keys')
    kmerKeys.sort()
    logging.debug('Writting Keys')
    for i in kmerKeys:
        handle.write('%s %d\n' % (i, kmers[i]))
    if path != '-':
        handle.close()

def printKmerCounts(kmers, path):
    logging.debug('Writing kmer counts to %s' % path)
    handle = sys.stdout
    if path != '-':
        handle = open(path)
    logging.debug('Writting Keys')
    for i in kmers:
        handle.write('%s %d\n' % (i, kmers[i]))
    if path != '-':
        handle.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s')
    if len(sys.argv) < 3:
        print 'Usage: %s K IN OUT' % os.path.basename(sys.argv[0])
        sys.exit(1)
    kmers = getKmerCounts(int(sys.argv[1]), sys.argv[2])
    printKmerCounts(kmers, sys.argv[3])

if __name__ == '__main__':
    main()
