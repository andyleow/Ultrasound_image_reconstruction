Matlab serialization
====================

Matlab object serialization functions built with undocumented mex functions
`mxSerialize` and `mxDeserialize`. The function can convert any ordinary matlab
variable into a uint8 array. These functions are unsupported and may change at
any time without notice in the future Matlab release. As of Matlab R2014a,
there was a big change in these hidden APIs, and this package supports only
`C++` MEX interface.

Build
-----

Use MEX command in Matlab.

    >> cd /path/to/matlab-serialization
    >> mex serialize.cc
    >> mex deserialize.cc

Alternatively, use the attached Makefile in UNIX.

Usage
-----

Add path to the `matlab-serialization` and compile before use.

    >> addpath('/path/to/matlab-serialization');

Use `serialize` to encode arbitrary matlab variable. The function returns
encoded variable as a `uint8` array.

    >> x = randn(1,4)

    x =

        0.7147   -0.2050   -0.1241    1.4897

    >> y = serialize(x)

    y =

      Columns 1 through 21

        0    1   73   77    0    0    0    0   14    0    0    0   80 ...

Use `deserialize` to retrieve the encoded variable.

    >> z = deserialize(y)

    z =

        0.7147   -0.2050   -0.1241    1.4897

License
-------

The code may be redistributed under the BSD clause 3 license.
