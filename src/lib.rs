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

#![feature(test)]

#[macro_use]
extern crate lazy_static;
extern crate test;

pub mod s2;

pub mod r3 {
    pub const X_AXIS: usize = 0;
    pub const Y_AXIS: usize = 1;
    pub const Z_AXIS: usize = 2;

    #[derive(Debug)]
    pub struct Vector {
        pub x: f64,
        pub y: f64,
        pub z: f64,
    }

    impl Vector {
        pub fn abs(&self) -> Vector {
            return Vector {
                x: self.x.abs(),
                y: self.y.abs(),
                z: self.z.abs(),
            };
        }
        pub fn largest_component(&self) -> usize {
            let t = self.abs();
            if t.x > t.y {
                if t.x > t.z {
                    return X_AXIS;
                }
                return Z_AXIS;
            }
            if t.y > t.z {
                return Y_AXIS;
            }
            return Z_AXIS;
        }
        pub fn norm(&self) -> f64 {
            return self.dot(self).sqrt();
        }
        pub fn normalize(&self) -> Vector {
            if self.x == 0.0 && self.y == 0.0 && self.z == 0.0 {
                return Vector {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                };
            }
            return self.mul(1.0 / self.norm());
        }
        pub fn mul(&self, m: f64) -> Vector {
            return Vector {
                x: self.x * m,
                y: self.y * m,
                z: self.z * m,
            };
        }
        pub fn dot(&self, other: &Vector) -> f64 {
            return self.x * other.x + self.y * other.y + self.z * other.z;
        }

        pub fn angle(&self, other: &Vector) -> f64 {
            return self.cross_prod(other).norm().atan2(self.dot(self));
        }
        pub fn cross_prod(&self, other: &Vector) -> Vector {
            return Vector {
                x: self.y * other.z - self.z * other.y,
                y: self.z * other.x - self.x * other.z,
                z: self.x * other.y - self.y * other.x,
            };
        }
    }
}

pub mod s1 {
    use std::f64::consts::PI;

    #[derive(Debug, PartialEq)]
    pub struct S1Angle(pub f64);

    impl S1Angle {
        pub fn from_radians(radians: f64) -> S1Angle {
            return S1Angle(radians);
        }
        pub fn from_degrees(degrees: f64) -> S1Angle {
            return S1Angle(degrees * (PI / 180.0));
        }

        pub fn radians(&self) -> f64 {
            return self.0;
        }
        pub fn degrees(&self) -> f64 {
            return self.0 * (180.0 / PI);
        }

        pub fn eq(&self, other: &S1Angle) -> bool {
            return self.0 == other.0;
        }
    }
}
