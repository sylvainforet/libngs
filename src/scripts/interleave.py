#!/usr/bin/python

import sys
import os.path

# This is a very primitive script without any ckecking or security

def interleave(path1, path2, pathOut):
    f1 = open(path1)
    f2 = open(path2)
    fo = open(pathOut, "w")
    while True:
        l1 = f1.readline()
        if not l1:
            break
        fo.write(l1)
        fo.write(f1.readline())
        fo.write(f1.readline())
        fo.write(f1.readline())
        fo.write(f2.readline())
        fo.write(f2.readline())
        fo.write(f2.readline())
        fo.write(f2.readline())
    fo.close()

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "Usage: %s IN1 IN2 OUT" % os.path.basename(sys.argv[0])
        sys.exit(1)
    interleave(sys.argv[1], sys.argv[2], sys.argv[3])
