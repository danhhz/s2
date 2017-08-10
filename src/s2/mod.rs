// Copyright 2005 Google Inc. All Rights Reserved.
// Copyright 2017 Daniel Harrison. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! Module s2 implements types and functions for working with geometry in S²
//! (spherical geometry).
//!
//! Its related modules, parallel to this one, are s1 (operates on S¹), r1
//! (operates on ℝ¹) and r3 (operates on ℝ³).
//!
//! This module provides types and functions for the S2 cell hierarchy and
//! coordinate systems. The S2 cell hierarchy is a hierarchical decomposition of
//! the surface of a unit sphere (S²) into 'cells'; it is highly efficient,
//! scales from continental size to under 1cm² and preserves spatial locality
//! (nearby cells have close IDs).
//!
//! A presentation that gives an overview of S2 is
//! https://docs.google.com/presentation/d/1Hl4KapfAENAOf4gv-pSngKwvS_jwNVHRPZTTDzXXn6Q/view.

#![deny(missing_docs)]

extern crate rand;

pub mod s2cellid;
pub mod s2latlng;

pub use self::s2cellid::*;
pub use self::s2latlng::*;
use super::r3;
use super::s1;

/// An S2Point represents a point on the unit sphere as a 3D vector. Usually
/// points are normalized to be unit length, but some methods do not require
/// this. Among other things, there are overloaded operators that make it
/// convenient to write arithmetic expressions (e.g. (1-x)*p1 + x*p2).
pub type S2Point = r3::Vector;

/// The S2 class is simply a namespace for constants and static utility
/// functions related to spherical geometry, such as area calculations and edge
/// intersection tests. The name "S2" is derived from the mathematical symbol
/// for the two-dimensional unit sphere (note that the "2" refers to the
/// dimension of the surface, not the space it is embedded in).
///
/// This class also defines a framework for decomposing the unit sphere into a
/// hierarchy of "cells". Each cell is a quadrilateral bounded by four
/// geodesics. The top level of the hierarchy is obtained by projecting the six
/// faces of a cube onto the unit sphere, and lower levels are obtained by
/// subdividing each cell into four children recursively.
///
/// This class specifies the details of how the cube faces are projected onto
/// the unit sphere. This includes getting the face ordering and orientation
/// correct so that sequentially increasing cell ids follow a continuous
/// space-filling curve over the entire sphere, and defining the transformation
/// from cell-space to cube-space in order to make the cells more uniform in
/// size.
///
/// This file also contains documentation of the various coordinate systems and
/// conventions used.
pub struct S2;

impl S2 {
    /// Return true if two points are within the given distance of each other
    /// (this is mainly useful for testing).
    pub fn approx_equals(a: &S2Point, b: &S2Point) -> bool {
        return a.angle(b) <= 1e-15;
    }

    /// A fast routine for converting floating-point numbers to integers. This
    /// mirrors the FastIntRound method in the c++ libaray, which works around
    /// the unavailability of llrint at the time the library was written.
    pub fn fast_int_round(x: f64) -> i64 {
        unsafe {
            return llrint(x);
        }
    }

    /// Convert an s or t value  to the corresponding u or v value. This is a
    /// non-linear transformation from [-1,1] to [-1,1] that attempts to make
    /// the cell sizes more uniform. This implemention matches
    /// S2_QUADRATIC_PROJECTION in the c++ library.
    pub fn st_to_uv(s: f64) -> f64 {
        if s >= 0.5 {
            return (1.0 / 3.0) * (4.0 * s * s - 1.0);
        } else {
            return (1.0 / 3.0) * (1.0 - 4.0 * (1.0 - s) * (1.0 - s));
        }
    }

    /// The inverse of the uv_to_st transformation. Note that it is not always
    /// true that uv_to_st(st_to_uv(x)) == x due to numerical errors.
    pub fn uv_to_st(u: f64) -> f64 {
        // TODO(dan): This uses the default S2_QUADRATIC_PROJECTION. Support the
        // others?
        if u >= 0.0 {
            return 0.5 * (1.0 + 3.0 * u).sqrt();
        } else {
            return 1.0 - 0.5 * (1.0 - 3.0 * u).sqrt();
        }
    }

    fn face(r: &r3::Vector) -> usize {
        let mut f = r.largest_component();
        match f {
            r3::X_AXIS if r.x < 0.0 => f += 3,
            r3::Y_AXIS if r.y < 0.0 => f += 3,
            r3::Z_AXIS if r.z < 0.0 => f += 3,
            _ => {}
        }
        return f;
    }

    /// Given a valid face for the given point p (meaning that dot product of p
    /// with the face normal is positive), return the corresponding u and v
    /// values (which may lie outside the range [-1,1]).
    fn valid_face_xyz_to_uv(face: usize, r: r3::Vector) -> (f64, f64) {
        // debug_assert!(r.dot(&S2::face_uv_to_xyz(face, 0.0, 0.0)) == 0.0);
        return match face {
            0 => (r.y / r.x, r.z / r.x),
            1 => (-r.x / r.y, r.z / r.y),
            2 => (-r.x / r.z, -r.y / r.z),
            3 => (r.z / r.x, r.y / r.x),
            4 => (r.z / r.y, -r.x / r.y),
            5 => (-r.y / r.z, -r.x / r.z),
            _ => unreachable!(),
        };
    }

    /// Convert (face, u, v) coordinates to a direction vector (not necessarily
    /// unit length).
    pub fn face_uv_to_xyz(face: usize, u: f64, v: f64) -> r3::Vector {
        return match face {
            0 => S2Point { x: 1.0, y: u, z: v },
            1 => S2Point {
                x: -u,
                y: 1.0,
                z: v,
            },
            2 => S2Point {
                x: -u,
                y: -v,
                z: 1.0,
            },
            3 => S2Point {
                x: -1.0,
                y: -v,
                z: -u,
            },
            4 => S2Point {
                x: v,
                y: -1.0,
                z: -u,
            },
            5 => S2Point {
                x: v,
                y: u,
                z: -1.0,
            },
            _ => unreachable!(),
        };
    }

    /// Convert a direction vector (not necessarily unit length) to (face, u, v)
    /// coordinates.
    pub fn xyz_to_face_uv(r: r3::Vector) -> (usize, f64, f64) {
        let f = S2::face(&r);
        let (u, v) = S2::valid_face_xyz_to_uv(f, r);
        return (f, u, v);
    }
}

extern "C" {
    fn llrint(a: f64) -> i64;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// This class defines various static functions that are useful for writing unit
    /// tests.
    pub struct S2Testing;

    impl S2Testing {
        /// Return a random unit-length vector.
        pub fn random_point() -> S2Point {
            return S2Point {
                x: 2.0 * rand::random::<f64>() - 1.0,
                y: 2.0 * rand::random::<f64>() - 1.0,
                z: 2.0 * rand::random::<f64>() - 1.0,
            }.normalize();
        }

        // Return a random cell id at the given level. The distribution is
        // uniform over the space of cell ids, but only approximately uniform
        // over the surface of the sphere.
        pub fn get_random_cell_id_of_level(level: usize) -> S2CellId {
            let face = rand::random::<usize>() % S2CellId::NUM_FACES;
            let pos = rand::random::<u64>() & ((1 << (2 * S2CellId::MAX_LEVEL)) - 1);
            return S2CellId::from_face_pos_level(face, pos, level);
        }
    }
}
