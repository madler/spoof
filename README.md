Synopsis
--------

_spoof_ will help you modify a message to make the CRC of that message equal
any desired value. _spoof_ does not need the message itself, but just the
length of the message, the exclusive-or of the message's current CRC and the
desired new CRC, and a set of bit locations in the message to potentially
modify. _spoof_ will then deliver a subset of those locations whose bits should
be inverted. Then the modified message will have the desired CRC.

Motivation
----------

The purpose of _spoof_ is to illustrate the extremely non-cryptographic nature
of a Cyclic Redundancy Check (CRC) as a signature. Since a CRC is a linear
operation on the message, it is easy to invert the operation to construct
messages with any desired CRC. _spoof_ runs in O(log(_n_)) time, where _n_ is
the length of the message.

Installation
------------

Simply compile and link spoof.c and fline.c with a standard C compiler and
library. "standard" here is defined as C99. spoof.c is a command line program
that takes input from stdin and produces output on stdout.

Tests
-----

This input, specifying a 4-bit CRC with polynomial x<sup>4</sup>+x+1 reflected,
an exclusive-or of the current and desired CRC of 1111, a message length of 89
(decimal) bytes, and four candidate bit locations to change in byte offset and
bit number:

    4 1 c
    f 89
    37 0
    41 0
    45 0
    49 0

will produce this output:

    invert these bits in the sequence:
    offset bit
        41 0

License
-------

This code is under the zlib license, permitting free commercial use.
