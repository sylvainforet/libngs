#!/usr/bin/python

letters = ['A', 'C', 'G', 'T']

def make_table():
    table = [0] * 256
    for i in xrange(256):
        bin_vals    = [0] * 4
        rev_vals    = [0] * 4
        bin_vals[0] = i & 3
        bin_vals[1] = (i >> 2) & 3
        bin_vals[2] = (i >> 4) & 3
        bin_vals[3] = (i >> 6) & 3
        rev_vals[0] = (~bin_vals[3]) & 3
        rev_vals[1] = (~bin_vals[2]) & 3
        rev_vals[2] = (~bin_vals[1]) & 3
        rev_vals[3] = (~bin_vals[0]) & 3
        #print '>>> --- ', bin_vals, rev_vals, [letters[x] for x in bin_vals], [letters[x] for x in rev_vals]
        revcomp     = rev_vals[0] | (rev_vals[1] << 2) | (rev_vals[2] << 4) | (rev_vals[3] << 6)
        table[i]    = revcomp
    return table

def print_table(table):
    print 'const unsigned char revcomp_table[256] ='
    print '{'
    for i in xrange(len(table)):
        if i == len(table) - 1:
            print '  [%3d] = %3d' % (i, table[i])
        else:
            print '  [%3d] = %3d,' % (i, table[i])
    print '};'

def main():
    table = make_table()
    print_table(table)

if __name__ == '__main__':
    main()
