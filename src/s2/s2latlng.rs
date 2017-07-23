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

use super::{r3, s1};

use std::f64::consts::PI;

#[derive(Debug, PartialEq)]
pub struct S2LatLng {
    pub coords: [f64; 2],
}

impl S2LatLng {
    pub fn from_point(p: &super::S2Point) -> S2LatLng {
        let ll = S2LatLng {
            coords: [
                S2LatLng::latitude(&p).radians(),
                S2LatLng::longitude(&p).radians(),
            ],
        };
        debug_assert!(ll.is_valid());
        return ll;
    }
    pub fn from_radians(lat: f64, lng: f64) -> S2LatLng {
        return S2LatLng { coords: [lat, lng] };
    }
    pub fn from_degrees(lat: f64, lng: f64) -> S2LatLng {
        return S2LatLng {
            coords: [
                s1::S1Angle::from_degrees(lat).0,
                s1::S1Angle::from_degrees(lng).0,
            ],
        };
    }

    pub fn latitude(p: &super::S2Point) -> s1::S1Angle {
        // We use atan2 rather than asin because the input vector is not
        // necessarily unit length, and atan2 is much more accurate than asin
        // near the poles.
        return s1::S1Angle::from_radians(p.z.atan2((p.x * p.x + p.y * p.y).sqrt()));
    }
    pub fn longitude(p: &super::S2Point) -> s1::S1Angle {
        // Note that atan2(0, 0) is defined to be zero.
        return s1::S1Angle::from_radians(p.y.atan2(p.x));
    }

    pub fn lat(&self) -> s1::S1Angle {
        return s1::S1Angle::from_radians(self.coords[0]);
    }
    pub fn lng(&self) -> s1::S1Angle {
        return s1::S1Angle::from_radians(self.coords[1]);
    }

    pub fn is_valid(&self) -> bool {
        return self.lat().radians().abs() <= PI / 2.0 && self.lng().radians().abs() <= PI;
    }

    pub fn normalized(&self) -> S2LatLng {
        let lat = (-PI / 2.0).max((PI / 2.0).min(self.lat().radians()));
        let lng: f64;
        unsafe {
            // drem(x, 2 * M_PI) reduces its argument to the range [-M_PI, M_PI]
            // inclusive, which is what we want here.
            lng = drem(self.lng().radians(), PI * 2.0);
        }
        return S2LatLng { coords: [lat, lng] };
    }

    pub fn to_point(&self) -> super::S2Point {
        let phi = self.lat().radians();
        let theta = self.lng().radians();
        let cosphi = phi.cos();
        return r3::Vector {
            x: theta.cos() * cosphi,
            y: theta.sin() * cosphi,
            z: phi.sin(),
        };
    }

    pub fn get_distance(&self, other: S2LatLng) -> s1::S1Angle {
        // This implements the Haversine formula, which is numerically stable
        // for small distances but only gets about 8 digits of precision for
        // very large distances (e.g. antipodal points).  Note that 8 digits is
        // still accurate to within about 10cm for a sphere the size of the
        // Earth.
        //
        // This could be fixed with another sin() and cos() below, but at that
        // point you might as well just convert both arguments to S2Points and
        // compute the distance that way (which gives about 15 digits of
        // accuracy for all distances).

        debug_assert!(self.is_valid());
        debug_assert!(other.is_valid());
        let lat1 = self.lat().radians();
        let lat2 = other.lat().radians();
        let lng1 = self.lng().radians();
        let lng2 = other.lng().radians();
        let dlat = (0.5 * (lat2 - lat1)).sin();
        let dlng = (0.5 * (lng2 - lng1)).sin();
        let x = dlat * dlat + dlng * dlng * lat1.cos() * lat2.cos();
        return s1::S1Angle::from_radians(2.0 * x.sqrt().atan2(0f64.max(1.0 - x).sqrt()));
    }
}

extern "C" {
    fn drem(a: f64, b: f64) -> f64;
}

#[cfg(test)]
mod tests {
    extern crate test;
    use super::S2LatLng;
    use super::s1::S1Angle;
    use super::super::{S2, S2Testing};
    use std::f64::consts::PI;

    #[test]
    fn basic() {
        let ll_rad = S2LatLng::from_radians(PI / 4.0, PI / 2.0);
        assert_eq!(PI / 4.0, ll_rad.lat().radians());
        assert_eq!(PI / 2.0, ll_rad.lng().radians());
        assert_eq!(true, ll_rad.is_valid());
        let ll_deg = S2LatLng::from_degrees(45.0, 90.0);
        assert_eq!(ll_rad, ll_deg);
        assert_eq!(true, ll_deg.is_valid());
        assert_eq!(false, S2LatLng::from_degrees(-91.0, 0.0).is_valid());
        assert_eq!(false, S2LatLng::from_degrees(0.0, 181.0).is_valid());

        let mut bad = S2LatLng::from_degrees(120.0, 200.0);
        assert_eq!(false, bad.is_valid());
        let mut better = bad.normalized();
        assert_eq!(true, better.is_valid());
        assert_eq!(S1Angle::from_degrees(90.0), better.lat());
        assert_eq!(
            S1Angle::from_degrees(-160.0).radians(),
            better.lng().radians()
        );

        bad = S2LatLng::from_degrees(-100.0, -360.0);
        assert_eq!(false, bad.is_valid());
        better = bad.normalized();
        assert_eq!(true, better.is_valid());
        assert_eq!(S1Angle::from_degrees(-90.0), better.lat());
        assert_eq!(0.0, better.lng().radians());
    }

    #[test]
    fn conversion() {
        // Test special cases: poles, "date line"
        assert_eq!(
            90.0,
            S2LatLng::from_point(&S2LatLng::from_degrees(90.0, 65.0).to_point())
                .lat()
                .degrees()
        );
        assert_eq!(
            -PI / 2.0,
            S2LatLng::from_point(&S2LatLng::from_radians(-PI / 2.0, 1.0).to_point())
                .lat()
                .radians()
        );
        assert_eq!(
            180.0,
            S2LatLng::from_point(&S2LatLng::from_degrees(12.2, 180.0).to_point())
                .lng()
                .degrees()
                .abs()
        );
        assert_eq!(
            PI,
            S2LatLng::from_point(&S2LatLng::from_radians(0.1, -PI).to_point())
                .lng()
                .radians()
                .abs()
        );

        // Test a bunch of random points.
        for _ in 0..100000 {
            let p = S2Testing::random_point();
            assert!(S2::approx_equals(&p, &S2LatLng::from_point(&p).to_point()));
        }
    }

    fn assert_f64_near(a: f64, b: f64, epsilon: f64) {
        assert!((a - b).abs() < epsilon);
    }

    #[test]
    fn distance() {
        assert_eq!(
            0.0,
            S2LatLng::from_degrees(90.0, 0.0)
                .get_distance(S2LatLng::from_degrees(90.0, 0.0))
                .radians()
        );
        assert_f64_near(
            77.0,
            S2LatLng::from_degrees(-37.0, 25.0)
                .get_distance(S2LatLng::from_degrees(-66.0, -155.0))
                .degrees(),
            1e-13,
        );
        assert_f64_near(
            115.0,
            S2LatLng::from_degrees(0.0, 165.0)
                .get_distance(S2LatLng::from_degrees(0.0, -80.0))
                .degrees(),
            1e-13,
        );
        assert_f64_near(
            180.0,
            S2LatLng::from_degrees(47.0, -127.0)
                .get_distance(S2LatLng::from_degrees(-47.0, 53.0))
                .degrees(),
            2e-6,
        );
    }

    #[bench]
    fn to_point(b: &mut test::Bencher) {
        let ll = S2LatLng::from_degrees(0x150bc888 as f64 * 1e-7, 0x5099d63f as f64 * 1e-7);
        b.iter(|| { return ll.to_point(); })
    }
}
