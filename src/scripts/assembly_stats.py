#!/usr/bin/python

import cStringIO

statsNames = ['name',
              'N75', 'N50', 'N25',
              'q75', 'q50', 'q25',
              'm0', 'm100', 'm200', 'm400', 'm800', 'm1000', 'm2000', 'm4000',
              's100', 's500', 's1000', 's2000']

def loadAssembly(path):
    handle = open(path)
    seqs = {}
    name = None
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            name = line
            seqs[name] = cStringIO.StringIO()
        elif name:
            seqs[name].write(line)
    handle.close()
    for name in seqs:
        tmp = seqs[name].getvalue()
        seqs[name].close()
        seqs[name] = tmp
    sizes = [len(x) for x in seqs.values()]
    sizes.sort()
    return sizes

def countLargerThan(sizes, minimum):
    largerSeqs = [x for x in sizes if x >= minimum]
    return len(largerSeqs)

def sumLargerThan(sizes, minimum):
    largerSeqs = [x for x in sizes if x >= minimum]
    return sum(largerSeqs)

def printHeader(handle):
    header         = '\t'.join(statsNames)
    handle.write(header + '\n')

def getN25_50_75(sizes):
    total = sum(sizes)
    tmp   = 0
    n25   = 0
    n50   = 0
    n75   = 0
    for i in sizes:
        tmp += i
        if tmp >= 0.25 * total and not n75:
            n75 = i
        if tmp >= 0.50 * total and not n50:
            n50 = i
        if tmp >= 0.75 * total and not n25:
            n25 = i
    return n25, n50, n75

def printAssemblyStats(sizes, name, handle):
    nContigs       = len(sizes)
    stats          = {}
    stats['name' ] = name
    stats['q75'  ] = sizes[nContigs / 4]
    stats['q50'  ] = sizes[nContigs / 2]
    stats['q25'  ] = sizes[(3 * nContigs) / 4]
    a, b, c        = getN25_50_75(sizes)
    stats['N25'  ] = a
    stats['N50'  ] = b
    stats['N75'  ] = c
    stats['m0'   ] = nContigs
    stats['m100' ] = countLargerThan(sizes, 100)
    stats['m200' ] = countLargerThan(sizes, 200)
    stats['m400' ] = countLargerThan(sizes, 400)
    stats['m800' ] = countLargerThan(sizes, 800)
    stats['m1000'] = countLargerThan(sizes, 1000)
    stats['m2000'] = countLargerThan(sizes, 2000)
    stats['m4000'] = countLargerThan(sizes, 4000)
    stats['s100' ] = sumLargerThan(sizes, 100)
    stats['s500' ] = sumLargerThan(sizes, 500)
    stats['s1000'] = sumLargerThan(sizes, 1000)
    stats['s2000'] = sumLargerThan(sizes, 2000)
    handle.write(name + '\t')
    for i in statsNames[1:-1]:
        handle.write('%d\t' % stats[i])
    handle.write('%d\n' % stats[statsNames[-1]])

def main():
    import sys
    import os.path
    if len(sys.argv) < 2:
        print 'Usage: %s FILE' % os.path.basename(sys.argv[0])
        sys.exit(1)
    printHeader(sys.stdout)
    for i in sys.argv[1:]:
        sizes = loadAssembly(i)
        printAssemblyStats(sizes, i, sys.stdout)

if __name__ == '__main__':
    main()

###
