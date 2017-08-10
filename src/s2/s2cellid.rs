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

//! An S2CellId is a 64-bit unsigned integer that uniquely identifies a cell in
//! the S2 cell decomposition.

extern crate rand;

use super::super::s2;

use std::cmp;

/// An S2CellId is a 64-bit unsigned integer that uniquely identifies a cell in
/// the S2 cell decomposition. It has the following format:
///
/// ```text
/// id = [face][face_pos]
///
/// face: a 3-bit number (range 0..5) encoding the cube face.
///
/// face_pos: a 61-bit number encoding the position of the center of this cell
/// along the Hilbert curve over this face.
/// ```
///
/// Sequentially increasing cell ids follow a continuous space-filling curve
/// over the entire sphere. They have the following properties:
///
/// - The id of a cell at level k consists of a 3-bit face number followed
///    by k bit pairs that recursively select one of the four children of
///    each cell. The next bit is always 1, and all other bits are 0.
///    Therefore, the level of a cell is determined by the position of its
///    lowest-numbered bit that is turned on (for a cell at level k, this
///    position is 2 * (kMaxLevel - k).)
///
/// - The id of a parent cell is at the midpoint of the range of ids spanned
///    by its children (or by its descendants at any level).
///
/// Leaf cells are often used to represent points on the unit sphere, and this
/// class provides methods for converting directly between these two
/// representations. For cells that represent 2D regions rather than discrete
/// points, it is better to use the S2Cell class.
///
/// This class is intended to be copied by value as desired.
#[derive(Debug, PartialEq)]
pub struct S2CellId {
    id: u64,
}

impl S2CellId {
    // Although only 60 bits are needed to represent the index of a leaf cell,
    // we need an extra bit in order to represent the position of the center of
    // the leaf cell along the Hilbert curve.

    /// The number of faces on the unit cube.
    pub const NUM_FACES: usize = 6;
    /// The number of levels needed to specify a leaf cell.
    pub const MAX_LEVEL: usize = 30;
    const FACE_BITS: usize = 3;
    const POS_BITS: usize = 2 * S2CellId::MAX_LEVEL + 1;
    const MAX_SIZE: i64 = 1 << S2CellId::MAX_LEVEL;

    /// Return an S2CellId with the given id.
    pub fn new(id: u64) -> S2CellId {
        return S2CellId { id: id };
    }

    /// Return a cell given its face (range 0..5), 61-bit Hilbert curve position
    /// within that face, and level (range 0..kMaxLevel). The given position
    /// will be modified to correspond to the Hilbert curve position at the
    /// center of the returned cell.
    pub fn from_face_pos_level(face: usize, pos: u64, level: usize) -> S2CellId {
        let cell = S2CellId { id: ((face as u64) << S2CellId::POS_BITS) + (pos | 1) };
        return cell.parent(level);
    }

    /// Return the i- or j-index of the leaf cell containing the given s- or
    /// t-value. Values are clamped appropriately.
    fn st_to_ij(s: f64) -> i64 {
        return cmp::max(0,
                        cmp::min(S2CellId::MAX_SIZE - 1,
                                 s2::S2::fast_int_round((S2CellId::MAX_SIZE as f64) * s - 0.5)));
    }

    /// Return a leaf cell given its cube face (range 0..5) and i- and
    /// j-coordinates (see s2.h).
    pub fn from_face_ij(face: usize, i: i64, j: i64) -> S2CellId {
        let mut n = (face as u64) << (S2CellId::POS_BITS - 1);
        let mut bits = face & SWAP_MASK;
        for k in (0..8).rev() {
            let mask = (1 << LOOKUP_BITS) - 1;
            bits += (((i >> k * LOOKUP_BITS) & mask) as usize) << (LOOKUP_BITS + 2);
            bits += (((j >> k * LOOKUP_BITS) & mask) as usize) << 2;
            bits = LOOKUP_POS[bits] as usize;
            n |= ((bits >> 2) as u64) << (k * 2 * LOOKUP_BITS);
            bits &= SWAP_MASK | INVERT_MASK;
        }
        return S2CellId { id: n * 2 + 1 };
    }

    /// Return the leaf cell containing the given point (a direction vector, not
    /// necessarily unit length).
    pub fn from_point(p: s2::S2Point) -> S2CellId {
        let (face, u, v) = s2::S2::xyz_to_face_uv(p);
        let i = S2CellId::st_to_ij(s2::S2::uv_to_st(u));
        let j = S2CellId::st_to_ij(s2::S2::uv_to_st(v));
        return S2CellId::from_face_ij(face, i, j);
    }

    /// Return the leaf cell containing the given normalized S2LatLng.
    pub fn from_lat_lng(ll: &s2::S2LatLng) -> S2CellId {
        return S2CellId::from_point(ll.to_point());
    }

    /// The 64-bit unique identifier for this cell.
    pub fn id(&self) -> u64 {
        return self.id;
    }

    /// Return true if id() represents a valid cell.
    pub fn is_valid(&self) -> bool {
        return self.face() < S2CellId::NUM_FACES && (self.lsb() & 0x1555555555555555) > 0;
    }

    /// Which cube face this cell belongs to, in the range 0..5.
    pub fn is_face(&self) -> bool {
        return (self.id & (S2CellId::lsb_for_level(0) - 1)) == 0;
    }

    /// Which cube face this cell belongs to, in the range 0..5.
    pub fn face(&self) -> usize {
        return ((self.id as u64) >> S2CellId::POS_BITS) as usize;
    }

    /// The position of the cell center along the Hilbert curve over this face,
    /// in the range 0..(2**kPosBits-1).
    pub fn pos(&self) -> u64 {
        return self.id as u64 & (!0u64 >> S2CellId::FACE_BITS);
    }

    /// Return the subdivision level of the cell (range 0..MAX_LEVEL).
    pub fn level(&self) -> usize {
        if self.is_leaf() {
            return S2CellId::MAX_LEVEL;
        }

        let mut x = self.id as i64; // TODO(dan): The c++ impl uses u32
        let mut level: isize = -1;
        if x != 0 {
            level += 16
        } else {
            x = self.id as i64 >> 32;
        }

        x &= -x;
        if x & 0x00005555 > 0 {
            level += 8;
        }
        if x & 0x00550055 > 0 {
            level += 4;
        }
        if x & 0x05050505 > 0 {
            level += 2;
        }
        if x & 0x11111111 > 0 {
            level += 1;
        }
        debug_assert!(level >= 0);
        debug_assert!(level <= S2CellId::MAX_LEVEL as isize);
        return level as usize;
    }

    /// Return true if this is a leaf cell (more efficient than checking whether
    /// level() == MAX_LEVEL).
    pub fn is_leaf(&self) -> bool {
        return self.id & 1 > 0;
    }

    /// Return the (face, i, j) coordinates for the leaf cell corresponding to
    /// this cell id. Since cells are represented by the Hilbert curve position
    /// at the center of the cell, the returned (i,j) for non-leaf cells will be
    /// a leaf cell adjacent to the cell center. If "orientation" is non-empty,
    /// also return the Hilbert curve orientation for the current cell.
    pub fn to_face_ij_orientation(&self) -> (usize, usize, usize, usize) {
        let f = self.face();
        let mut i = 0usize;
        let mut j = 0usize;
        let mut orientation = f & SWAP_MASK;
        let mut nbits = S2CellId::MAX_LEVEL - 7 * LOOKUP_BITS; // first iteration

        for k in (0..8).rev() {
            orientation += (((self.id >> (k * 2 * LOOKUP_BITS + 1)) &
                             ((1 << (2 * nbits)) - 1)) << 2) as usize;
            orientation = LOOKUP_IJ[orientation] as usize;
            i += (orientation >> (LOOKUP_BITS + 2)) << (k * LOOKUP_BITS);
            j += ((orientation >> 2) & ((1 << LOOKUP_BITS) - 1)) << (k * LOOKUP_BITS);
            orientation &= SWAP_MASK | INVERT_MASK;
            nbits = LOOKUP_BITS // following iterations;
        }

        if self.lsb() & 0x1111111111111110 != 0 {
            orientation &= SWAP_MASK;
        }
        return (f, i, j, orientation);
    }

    /// Return the (face, si, ti) coordinates of the center of the cell. Note
    /// that although (si,ti) coordinates span the range [0,2**31] in general,
    /// the cell center coordinates are always in the range [1,2**31-1] and
    /// therefore can be represented using a signed 32-bit integer.
    pub fn get_center_si_ti(&self) -> (usize, usize, usize) {
        let (face, i, j, _) = self.to_face_ij_orientation();
        let delta = if self.is_leaf() {
            1
        } else {
            unimplemented!();
        };
        return (face, 2 * i + delta, 2 * j + delta);
    }

    /// Return the direction vector corresponding to the center of the given
    /// cell. The vector returned by to_point_raw is not necessarily unit length.
    pub fn to_point_raw(&self) -> s2::S2Point {
        let (face, si, ti) = self.get_center_si_ti();
        return s2::S2::face_uv_to_xyz(face,
                                      s2::S2::st_to_uv((0.5 / S2CellId::MAX_SIZE as f64) *
                                                       si as f64),
                                      s2::S2::st_to_uv((0.5 / S2CellId::MAX_SIZE as f64) *
                                                       ti as f64));
    }

    /// Return the S2LatLng corresponding to the center of the given cell.
    pub fn to_lat_lng(&self) -> s2::S2LatLng {
        return s2::S2LatLng::from_point(&self.to_point_raw());
    }

    /// Return the lowest-numbered bit that is on for this cell id, which is
    /// equal to (uint64(1) << (2 * (kMaxLevel - level))). So for example,
    /// a.lsb() <= b.lsb() if and only if a.level() >= b.level(), but the first
    /// test is more efficient.
    pub fn lsb(&self) -> u64 {
        return ((self.id as i64) & -(self.id as i64)) as u64;
    }

    /// Return the lowest-numbered bit that is on for cells at the given level.
    fn lsb_for_level(level: usize) -> u64 {
        return 1 << 2 * (S2CellId::MAX_LEVEL - level);
    }

    /// Return the cell at the given level (which must be less than or equal to
    /// the current level).
    pub fn parent(&self, level: usize) -> S2CellId {
        debug_assert!(self.is_valid());
        debug_assert!(level <= self.level());
        let lsb = S2CellId::lsb_for_level(level);
        return S2CellId { id: (self.id & -(lsb as i64) as u64) | lsb };
    }

    /// Return the cell at the previous level (which must be less than or equal
    /// to the current level).
    pub fn immediate_parent(&self) -> S2CellId {
        debug_assert!(self.is_valid());
        debug_assert!(!self.is_face());
        let new_lsb = self.lsb() << 2;
        return S2CellId { id: (self.id & -(new_lsb as i64) as u64) | new_lsb };
    }
}

const LOOKUP_BITS: usize = 4;
const SWAP_MASK: usize = 0x01;
const INVERT_MASK: usize = 0x02;
const LOOKUP_SIZE: usize = 1 << (2 * LOOKUP_BITS + 2);
const POS_TO_IJ: [[usize; 4]; 4] = [
    [0, 1, 3, 2], // canonical order:    (0,0), (0,1), (1,1), (1,0)
    [0, 2, 3, 1], // axes swapped:       (0,0), (1,0), (1,1), (0,1)
    [3, 2, 0, 1], // bits inverted:      (1,1), (1,0), (0,0), (0,1)
    [3, 1, 0, 2], // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
];
const POS_TO_ORIENTATION: [usize; 4] = [SWAP_MASK, 0, 0, INVERT_MASK | SWAP_MASK];

lazy_static! {
    // TODO(dan): These each end up being 4KB. Consider codegen-ing them so
    // clients don't have to deal with the lazy init.
    static ref LOOKUP_POS: [u16; LOOKUP_SIZE] = {
       let mut lookup_pos = [0u16; LOOKUP_SIZE];
       let mut lookup_ij = [0u16; LOOKUP_SIZE];
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, 0, 0, 0);
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, SWAP_MASK, 0, SWAP_MASK);
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, INVERT_MASK, 0, INVERT_MASK);
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, SWAP_MASK|INVERT_MASK, 0, SWAP_MASK|INVERT_MASK);
       return lookup_pos;
    };
    // TODO(dan): Stop doing the lazy init twice.
    static ref LOOKUP_IJ: [u16; LOOKUP_SIZE] = {
       let mut lookup_pos = [0u16; LOOKUP_SIZE];
       let mut lookup_ij = [0u16; LOOKUP_SIZE];
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, 0, 0, 0);
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, SWAP_MASK, 0, SWAP_MASK);
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, INVERT_MASK, 0, INVERT_MASK);
       init_lookup_cell(&mut lookup_pos, &mut lookup_ij, 0, 0, 0, SWAP_MASK|INVERT_MASK, 0, SWAP_MASK|INVERT_MASK);
       return lookup_ij;
    };
}

fn init_lookup_cell(lookup_pos: &mut [u16; LOOKUP_SIZE],
                    lookup_ij: &mut [u16; LOOKUP_SIZE],
                    level: usize,
                    i: usize,
                    j: usize,
                    orig_orientation: usize,
                    pos: usize,
                    orientation: usize) {
    if level == LOOKUP_BITS {
        let ij = (i << LOOKUP_BITS) + j;
        lookup_pos[(ij << 2) + orig_orientation] = ((pos << 2) + orientation) as u16;
        lookup_ij[(pos << 2) + orig_orientation] = ((ij << 2) + orientation) as u16;
        return;
    }

    let (nlevel, ni, nj, npos) = (level + 1, i << 1, j << 1, pos << 2);
    let r = POS_TO_IJ[orientation];
    init_lookup_cell(lookup_pos,
                     lookup_ij,
                     nlevel,
                     ni + (r[0] >> 1),
                     nj + (r[0] & 1),
                     orig_orientation,
                     npos,
                     orientation ^ POS_TO_ORIENTATION[0]);
    init_lookup_cell(lookup_pos,
                     lookup_ij,
                     nlevel,
                     ni + (r[1] >> 1),
                     nj + (r[1] & 1),
                     orig_orientation,
                     npos + 1,
                     orientation ^ POS_TO_ORIENTATION[1]);
    init_lookup_cell(lookup_pos,
                     lookup_ij,
                     nlevel,
                     ni + (r[2] >> 1),
                     nj + (r[2] & 1),
                     orig_orientation,
                     npos + 2,
                     orientation ^ POS_TO_ORIENTATION[2]);
    init_lookup_cell(lookup_pos,
                     lookup_ij,
                     nlevel,
                     ni + (r[3] >> 1),
                     nj + (r[3] & 1),
                     orig_orientation,
                     npos + 3,
                     orientation ^ POS_TO_ORIENTATION[3]);
}

#[cfg(test)]
mod tests {
    use super::S2CellId;
    use super::super::S2LatLng;
    use super::super::tests::S2Testing;

    fn get_cellid(lat_degrees: f64, lng_degrees: f64) -> S2CellId {
        let id = S2CellId::from_lat_lng(&S2LatLng::from_degrees(lat_degrees, lng_degrees));
        return id;
    }

    #[test]
    fn default_constructor() {
        let id = S2CellId::new(0);
        assert_eq!(0, id.id());
        assert_eq!(false, id.is_valid());
    }

    #[test]
    fn face_definitions() {
        assert_eq!(0, get_cellid(0.0, 0.0).face());
        assert_eq!(1, get_cellid(0.0, 90.0).face());
        assert_eq!(2, get_cellid(90.0, 0.0).face());
        assert_eq!(3, get_cellid(0.0, 180.0).face());
        assert_eq!(4, get_cellid(0.0, -90.0).face());
        assert_eq!(5, get_cellid(-90.0, 0.0).face());
    }

    #[test]
    fn parent_child_relationships() {
        let id = S2CellId::from_face_pos_level(3, 0x12345678, S2CellId::MAX_LEVEL - 4);
        assert_eq!(true, id.is_valid());
        assert_eq!(3, id.face());
        assert_eq!(0x12345700, id.pos());
        assert_eq!(S2CellId::MAX_LEVEL - 4, id.level());
        assert_eq!(false, id.is_leaf());

        // assert_eq!(0x12345610, id.child_begin(id.level() + 2).pos());
        // assert_eq!(0x12345640, id.child_begin().pos());
        assert_eq!(0x12345400, id.immediate_parent().pos());
        assert_eq!(0x12345000, id.parent(id.level() - 2).pos());

        // TODO(dan): Port the rest of this test.
    }

    #[test]
    fn inverses() {
        // Check the conversion of random leaf cells to S2LatLngs and back.
        for _ in 0..10000 {
            let id = S2Testing::get_random_cell_id_of_level(S2CellId::MAX_LEVEL);
            assert_eq!(true, id.is_leaf());
            assert_eq!(S2CellId::MAX_LEVEL, id.level());
            let center = id.to_lat_lng();

            assert_eq!(id.id(), S2CellId::from_lat_lng(&center).id());
        }
    }

    #[test]
    fn cellid_latlng() {
        let tests = [(0x47a1cbd595522b39, 49.703498679, 11.770681595),
                     (0x46525318b63be0f9, 55.685376759, 12.588490937),
                     (0x52b30b71698e729d, 45.486546517, -93.449700022),
                     (0x46ed8886cfadda85, 58.299984854, 23.049300056),
                     (0x3663f18a24cbe857, 34.364439040, 108.330699969),
                     (0x10a06c0a948cf5d, -30.694551352, -30.048758753),
                     (0x2b2bfd076787c5df, -25.285264027, 133.823116966),
                     (0xb09dff882a7809e1, -75.000000031, 0.000000133),
                     (0x94daa3d000000001, -24.694439215, -47.537363213),
                     (0x87a1000000000001, 38.899730392, -99.901813021),
                     (0x4fc76d5000000001, 81.647200334, -55.631712940),
                     (0x3b00955555555555, 10.050986518, 78.293170610),
                     (0x1dcc469991555555, -34.055420593, 18.551140038),
                     (0xb112966aaaaaaaab, -69.219262171, 49.670072392)];
        for test in tests.iter() {
            let id = S2CellId { id: test.0 };
            let ll = S2LatLng::from_degrees(test.1, test.2);
            assert_eq!(id, S2CellId::from_lat_lng(&ll));
            assert!(id.to_lat_lng().get_distance(ll).0 < 1e-9);
        }
    }
}
