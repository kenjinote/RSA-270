# GNFS Educational Framework (C++)

This repository contains an educational C++ framework modeling the General Number Field Sieve (GNFS), the state-of-the-art algorithm for integer factorization.

## Overview
The GNFS algorithm is responsible for breaking large RSA challenges. This framework breaks down the complex mathematical operations into the five distinct phases used by modern factoring tools (such as CADO-NFS):

1. **Polynomial Selection**: Generating optimal algebraic and rational polynomials using base-m expansion.
2. **Special-q Lattice Sieving**: Efficient caching and relation finding through factor base grids.
3. **Filtering**: Singleton removal and clique merging for matrix reduction.
4. **Block Wiedemann Algorithm (Linear Algebra over GF(2))**: Finding null spaces (perfect squares) in sparse matrices.
5. **Algebraic Square Root**: Computing roots in the algebraic number field and mapping back via homomorphism to find the GCD.

## Dependencies
- C++14 or higher
- [Boost C++ Libraries](https://www.boost.org/) (Specifically `boost::multiprecision` for large integers)

## Disclaimer
This is an educational mock/framework designed to demonstrate the data structures and pipeline architecture of GNFS. It uses dummy relations for the sieving phase to allow the pipeline to execute instantly for demonstration purposes.

## Usage
Compile with any standard C++ compiler linking Boost:
```bash
cl.exe /EHsc /O2 /utf-8 -I<path_to_boost> gnfs.cpp
```
Execute the resulting binary to see the pipeline simulation.
