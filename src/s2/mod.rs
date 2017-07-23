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

extern crate rand;

pub mod s2cellid;
pub mod s2latlng;

pub use self::s2cellid::*;
pub use self::s2latlng::*;
use super::r3;
use super::s1;


pub type S2Point = r3::Vector;

pub struct S2;

impl S2 {
    pub const MAX_CELL_LEVEL: usize = 30;

    pub fn approx_equals(a: &S2Point, b: &S2Point) -> bool {
        return a.angle(b) <= 1e-15;
    }

    pub fn fast_int_round(x: f64) -> i64 {
        // TODO(dan): This is not the same logic.
        return x.round() as i64;
    }

    pub fn st_to_uv(s: f64) -> f64 {
        if s >= 0.5 {
            return (1.0 / 3.0) * (4.0 * s * s - 1.0);
        } else {
            return (1.0 / 3.0) * (1.0 - 4.0 * (1.0 - s) * (1.0 - s));
        }
    }

    pub fn uv_to_st(u: f64) -> f64 {
        // TODO(dan): This uses the default S2_QUADRATIC_PROJECTION. Support the
        // others?
        if u >= 0.0 {
            return 0.5 * (1.0 + 3.0 * u).sqrt();
        } else {
            return 1.0 - 0.5 * (1.0 - 3.0 * u).sqrt();
        }
    }

    pub fn face(r: &r3::Vector) -> usize {
        let mut f = r.largest_component();
        match f {
            r3::X_AXIS if r.x < 0.0 => f += 3,
            r3::Y_AXIS if r.y < 0.0 => f += 3,
            r3::Z_AXIS if r.z < 0.0 => f += 3,
            _ => {}
        }
        return f;
    }

    pub fn valid_face_xyz_to_uv(face: usize, r: r3::Vector) -> (f64, f64) {
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

    pub fn xyz_to_face_uv(r: r3::Vector) -> (usize, f64, f64) {
        let f = S2::face(&r);
        let (u, v) = S2::valid_face_xyz_to_uv(f, r);
        return (f, u, v);
    }
}

pub struct S2Testing;

impl S2Testing {
    pub fn random_point() -> S2Point {
        return S2Point {
            x: 2.0 * rand::random::<f64>() - 1.0,
            y: 2.0 * rand::random::<f64>() - 1.0,
            z: 2.0 * rand::random::<f64>() - 1.0,
        }.normalize();
    }
    pub fn get_random_cell_id_of_level(level: usize) -> S2CellId {
        let face = rand::random::<usize>() % S2CellId::NUM_FACES;
        let pos = rand::random::<u64>() & ((1 << (2 * S2CellId::MAX_LEVEL)) - 1);
        return S2CellId::from_face_pos_level(face, pos, level);
    }
}
