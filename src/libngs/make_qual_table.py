#!/usr/bin/python

for i in xrange(128):
    if i >= 64 and i < 104:
        print '%-2d, /* %-3d -- %c */' % (i - 64 , i, i)
    else:
        print '0 , /* %-3d */' % i
