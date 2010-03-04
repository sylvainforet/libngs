#!/usr/bin/python


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
            seqs[name] = ''
        elif name:
            seqs[name] += line
    handle.close()
    sizes = [len(x) for x in seqs.values()]
    sizes.sort()
    return sizes

def countLargerThan(sizes, minimum):
    largerSeqs = [x for x in sizes if x >= minimum]
    return len(largerSeqs)

def sumLargerThan(sizes, minimum):
    largerSeqs = [x for x in sizes if x >= minimum]
    return sum(largerSeqs)

def printAssemblyStats(sizes, name, handle):
    nContigs       = len(sizes)
    stats          = {}
    statsNames     = ['name',
                      'n75', 'n50', 'n25',
                      'm0', 'm100', 'm200', 'm400', 'm800', 'm1000', 'm2000', 'm4000',
                      's100', 's500', 's1000', 's2000']
    stats['name' ] = name
    stats['n75'  ] = sizes[nContigs / 4]
    stats['n50'  ] = sizes[nContigs / 2]
    stats['n25'  ] = sizes[(3 * nContigs) / 4]
    stats['m0'   ] = nContigs
    stats['m100' ] = countLargerThan(sizes, 100)
    stats['m200' ] = countLargerThan(sizes, 200)
    stats['m400' ] = countLargerThan(sizes, 400)
    stats['m800' ] = countLargerThan(sizes, 800)
    stats['m1000'] = countLargerThan(sizes, 1000)
    stats['m2000'] = countLargerThan(sizes, 2000)
    stats['m4000'] = countLargerThan(sizes, 4000)
    stats['s100']  = sumLargerThan(sizes, 100)
    stats['s500']  = sumLargerThan(sizes, 500)
    stats['s1000'] = sumLargerThan(sizes, 1000)
    stats['s2000'] = sumLargerThan(sizes, 2000)
    header         = '\t'.join(statsNames)
    handle.write(header + '\n')
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
    sizes = loadAssembly(sys.argv[1])
    printAssemblyStats(sizes, sys.argv[1], sys.stdout)

if __name__ == '__main__':
    main()

###
