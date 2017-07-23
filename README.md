# S2

[![Build Status](https://travis-ci.org/danhhz/s2.svg?branch=master)](https://travis-ci.org/danhhz/s2)

This is a library for manipulating geometric shapes. Unlike many geometry
libraries, S2 is primarily designed to work with _spherical geometry_, i.e.,
shapes drawn on a sphere rather than on a planar 2D map. (In fact, the name S2
is derived from the mathematical notation for the unit sphere.) This makes it
especially suitable for working with geographic data.

# Status of the Rust Library

This library is principally a port of [the C++ S2 library], adapting to Rust
idioms where it makes sense. We detail the progress of this port below relative
to that C++ library.

## [ℝ¹](https://godoc.org/github.com/golang/geo/r1) - One-dimensional Cartesian coordinates

Not started.

## [ℝ²](https://godoc.org/github.com/golang/geo/r2) - Two-dimensional Cartesian coordinates

Not started.

## [ℝ³](https://godoc.org/github.com/golang/geo/r3) - Three-dimensional Cartesian coordinates

Partial implemention of R3Vector.

## [S¹](https://godoc.org/github.com/golang/geo/s1) - Circular Geometry

Not started.

## [S²](https://godoc.org/github.com/golang/geo/s2) - Spherical Geometry

Partial implementations of:

- S2CellId
- S2LatLng
- S2Point

[the C++ S2 library]: https://code.google.com/archive/p/s2-geometry-library/
