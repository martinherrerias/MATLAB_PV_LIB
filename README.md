# MATLAB_PV_LIB

An ongoing effort to shape my own copy of [sandialabs/MATLAB_PV_LIB](https://github.com/sandialabs/MATLAB_PV_LIB), which had diverged over the years (before I had heard of version control), as a fork of the original. 

Functions that are different are currently just named differently, and placed in the `mod` folder. E.g. `mod/pvlmod_x` for the original `pvl_x`. Most of it are trivial changes (input format and parsing), that should be merged and marked as changes from to the originals.

Cases where changes are more dramatic (e.g. modifications to the Perez model) should perhaps be moved outside of the library, or renamed.
