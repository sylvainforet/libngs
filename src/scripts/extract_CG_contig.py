#!/usr/bin/python

import sys
import os.path


def parseCG(path, name):
    handle = open(path)
    inContig = False
    for line in handle:
        if line[0] == '>':
            if inContig:
                #inContig = False
                break
            elif line.find(name) > -1:
                inContig = True
                print line.strip()
        else:
            if inContig:
                print line.strip()
    handle.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage: %s FILE NAME" % os.path.basename(sys.argv[0])
        sys.exit(1)
    parseCG(sys.argv[1], sys.argv[2])
