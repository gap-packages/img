[![CI](https://github.com/gap-packages/img/actions/workflows/CI.yml/badge.svg)](https://github.com/gap-packages/img/actions/workflows/CI.yml)
[![Code Coverage](https://codecov.io/github/gap-packages/img/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/img)

# The IMG package

This is the README file for the GAP package "IMG".

This package implements Iterated Monodromy groups.

This release (0.2.3) contains basic functionality for computing with IMGs,
though some rough edges will still need to be smoothed. In particular,
the TODO file describes plans for future development, including
computations of presentations for self-similar groups, interval
arithmetic for complex dynamics, etc.

The package is distributed in source form, and does not require
anything else than a running GAP 4.9 or later. For updates, check
<https://github.com/gap-packages/img/>

The package requires some compiled modules. For this, the commands
'./configure' and 'make' in the package's main directory should create
the required binary files.

To visualize complex maps and marked spheres, a web-based interface
may be installed. Again in the package's main directly, issue

    git clone https://github.com/laurentbartholdi/rsserver.git

This service will require a functional 'node' (javascript runtime).

To use the package, start GAP and type

    LoadPackage("IMG");

The "IMG" package banner should appear on the screen.

For details on how to use the IMG package, please consult the 
documentation. Though this is usually not necessary, it may be
recompiled, after the command

    LoadPackage("IMG");
    
by invoking

    DOC@IMG();

at the GAP prompt. The documentation will then be available in the
`doc` subdirectory (view the file `manual.pdf` via a PDF viewer).

An external module (DLL) provides floating-point functions used for
complex dynamics calculations. To compile it, first make sure that
your system has C and JAVA compilers. IMG needs some external
libraries, but they are included in the subdirectory `extern`.
If your system has them, their location can be specified using the
arguments `--with-levmar` and `--with-libdogleg` to the `./configure`
script. The option `--with-gaproot` tells IMG where to look for the
GAP installation; by default, IMG searches for it in `../..`.
The library `cholmod` is also required (often part of `suitesparse`),
but is not provided.

In all cases, then run

    ./configure && make
    
in IMG's root directory.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or any
later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program, in the file COPYING.  If not, see
<https://www.gnu.org/licenses/>.

  Laurent Bartholdi, Göttingen, 19 November 2012
