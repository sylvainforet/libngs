#!/bin/sh

set -x
aclocal -I build
libtoolize --force --copy
autoheader
automake --add-missing --copy
autoconf

