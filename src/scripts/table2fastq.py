#!/usr/bin/python

import gzip
import sys
import os.path

def parse(path):
    if path[-3:] == '.gz':
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        if fields[10] == '0':
            continue
        name = '%s%s::%s:%s:%s:%s#%s/%s' % tuple(fields[:8])
        print '@%s' % name
        print fields[8].replace('.', 'N')
        print '+%s' % name
        print fields[9].replace('.', 'N')
    handle.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: %s FILE' % os.path.basename(sys.argv[0])
        sys.exit(1)
    parse(sys.argv[1])
