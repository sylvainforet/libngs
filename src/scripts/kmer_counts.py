#!/usr/bin/python

import cStringIO
import gzip
import logging
import optparse
import os.path
import string
import sys

rev_table = string.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revComp(s):
    return s.translate(rev_table)[::-1]

def countWords(kmers, seq, k, revcomp):
    for i in xrange(len(seq) - k + 1):
        word = seq[i:i + k]
        if not 'N' in word:
            if not word in kmers:
                kmers[word] = 0
            kmers[word] += 1
            if revcomp:
                word = revComp(word)
                if not word in kmers:
                    kmers[word] = 0
                kmers[word] += 1

def getKmerCountsFasta(path, options):
    freqrep = options.freqrep
    revcomp = options.revcomp
    k       = options.k
    logging.debug('Parsing %s' % path)
    handle = sys.stdin
    if path != '-':
        if path[:-3] == '.gz':
            handle = gzip.open(path)
        else:
            handle = open(path)
    kmers = {}
    buf   = None
    n     = 0
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            if buf:
                countWords(kmers, buf.getvalue(), k, revcomp)
                buf.close()
                n += 1
                if n % freqrep == 0:
                    logging.debug('Parsed %d sequences' % n)
            buf = cStringIO.StringIO()
        elif buf:
            buf.write(line)
    if buf:
        countWords(kmers, buf.getvalue(), k, revcomp)
        buf.close()
        n += 1
        if n % freqrep == 0:
            logging.debug('Parsed %d sequences' % n)
    if path != '-':
        handle.close()
    return kmers

def getKmerCountsFastq(path, options):
    freqrep = options.freqrep
    revcomp = options.revcomp
    k       = options.k
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
        countWords(kmers, line2, k, revcomp)
        n += 1
        if n % freqrep == 0:
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
        handle = open(path, "w")
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
        handle = open(path, "w")
    logging.debug('Writting Keys')
    for i in kmers:
        handle.write('%s %d\n' % (i, kmers[i]))
    if path != '-':
        handle.close()

def getOptions():
    '''
    Command line options processing
    '''
    usage = 'Usage: %prog [options] IN OUT'
    parser = optparse.OptionParser(usage)
    # Misc options
    parser.add_option('-l',
                      '--logging',
                      dest='logging',
                      default='WARNING',
                      choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
                      metavar='LEVEL',
                      help='Logging level to display')
    parser.add_option('-r',
                      '--revcomp',
                      dest='revcomp',
                      action='store_true',
                      help='Also scan the reverse complement')
    parser.add_option('-s',
                      '--sort',
                      dest='sort',
                      action='store_true',
                      help='Print kmers in alphabetical order')
    parser.add_option('-q',
                      '--fastq',
                      dest='fastq',
                      action='store_true',
                      help='Input is in fastq format')
    parser.add_option('-e',
                      '--freqrep',
                      dest='freqrep',
                      type='int',
                      default=1000000,
                      help='Verbose-report frequency')
    parser.add_option('-k',
                      '--kmersize',
                      dest='k',
                      type='int',
                      default=0,
                      help='Kmer size')
    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_usage(sys.stderr)
        sys.exit(1)
    if options.k < 1:
        sys.stderr.write("[ERROR] The word size must be strictly positive (found k=%d)\n" % options.k)
        sys.exit(1)
    logging.basicConfig(level=logging.__dict__[options.logging],
                        format='%(asctime)s %(levelname)s %(message)s')
    return options, args

def main():
    options, args = getOptions()
    logging.debug('Starting ...')
    if options.fastq:
        kmers = getKmerCountsFastq(args[0], options)
    else:
        kmers = getKmerCountsFasta(args[0], options)
    if options.sort:
        printKmerCountsSorted(kmers, args[1])
    else:
        printKmerCounts(kmers, args[1])
    logging.debug('All Done')

if __name__ == '__main__':
    main()
