# Euler1D

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fergu.github.io/Euler1D.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fergu.github.io/Euler1D.jl/dev/)
[![Build Status](https://github.com/fergu/Euler1D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fergu/Euler1D.jl/actions/workflows/CI.yml?query=branch%3Amain)

Caution: The pre-1.0 version of this package is undergoing testing. While no major changes to functionality are planned, the API, particularly the names of output fields, may not be completely stable.

This is a package that solves the one-dimensional compressible Euler equations on a Lagrangian (moving) mesh and assuming an ideal gas equation of state. 

The primary contribution of this package is intended to be as a demonstration of how a hydrodynamics code can be implemented entirely in Julia, including problem setup, the simulation code itself, and results postprocessing. In this way, there is no need to maintain a separate parser to interpret "input decks" for configuring the intial problem state. The language used to set up the problem is the same as the one used to perform the simulation. As long as the language supports it, it can be used with no modification to this package. Ideally, this can all be done while maintaining relative parity in execution time with similar codes written in compiled languages (C, C++, Fortran, etc).
