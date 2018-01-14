Synopsis
--------

_spoof_ will assist you in modifying a message to make the CRC of that message
equal any desired value. _spoof_ does not need the message itself, but just the
length of the message, the exclusive-or of the message's current CRC and the
desired new CRC, and a set of bit locations in the message to potentially
modify. _spoof_ will then deliver a subset of those locations whose bits should
be inverted. Then the modified message will have the desired CRC.

_ruse_ will modify a small number of bits in a file such that the specified CRC
of the file is unchanged. The least number of bits possible will be changed,
starting at a randomly chosen location in the file. To do this, _ruse_ uses a
provided dictionary of CRC codewords in the file codewords.txt, pulled from the
work of Philip Koopman, which he makes available online. If the specified CRC
is not covered in codewords.txt, then the CRC polynomial is used as the
codeword.

Motivation
----------

The purpose of _spoof_ is to illustrate the extremely non-cryptographic nature
of a Cyclic Redundancy Check (CRC) as a signature. Since a CRC is a linear
operation on the message, it is easy to invert the operation to construct
messages with any desired CRC. _spoof_ runs in O(log(_n_)) time, where _n_ is
the length of the message. Similarly, _ruse_ applies codewords found through
brute-force searches for bit patterns that leave CRCs with given polynomials
unchanged.

Installation
------------

Simply compile and link spoof.c and fline.c with a standard C99 compiler, and
compile ruse.cc with a standard C++11 compiler. spoof.c is a command-line
program that takes input from stdin and produces output on stdout. ruse is a
command-line program that modifies the named file. The instructions for each
are near the start of the respective source files.

Test
----

This spoof input, specifying a 4-bit CRC with polynomial x<sup>4</sup>+x+1
reflected, an exclusive-or of the current and desired CRC of 1111, a message
length of 89 (decimal) bytes, and four candidate bit locations to change in
byte offset and bit number:

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

This ruse command:

    ./ruse 32 4c11db7 1 file

will modify some number of bits in file, leaving the standard CRC-32 unchanged.
The script getcodes can be used to pull the latest codewords from Philip Koopman's
website:

    ./getcodes > codewords.txt

License
-------

This code is under the zlib license, permitting free commercial use.
