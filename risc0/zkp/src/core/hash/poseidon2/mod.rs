// Copyright 2023 RISC Zero, Inc.
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

//! An implementation of Poseidon2 targeting the Baby Bear.

pub(crate) mod consts;
mod rng;

use alloc::{boxed::Box, rc::Rc, vec::Vec};

use risc0_core::field::{
    baby_bear::{BabyBear, BabyBearElem, BabyBearExtElem, Elem},
    ExtElem,
};

use self::consts::{
    M_INT_DIAG, M_EXT, ROUNDS_HALF_FULL, ROUNDS_PARTIAL, ROUND_CONSTANTS,
};
pub use self::{consts::CELLS, rng::Poseidon2Rng};
use super::{HashFn, HashSuite, Rng, RngFactory};
use crate::core::digest::{Digest, DIGEST_WORDS};

/// The 'rate' of the sponge, i.e. how much we can safely add/remove per mixing.
pub const CELLS_RATE: usize = 16;

/// The size of the hash output in cells (~ 248 bits)
pub const CELLS_OUT: usize = 8;

/// A hash implemention for Poseidon
struct Poseidon2HashFn;

impl HashFn<BabyBear> for Poseidon2HashFn {
    fn hash_pair(&self, a: &Digest, b: &Digest) -> Box<Digest> {
        let both: Vec<BabyBearElem> = a
            .as_words()
            .iter()
            .chain(b.as_words())
            .map(|w| BabyBearElem::new_raw(*w))
            .collect();
        assert!(both.len() == 16);
        to_digest(unpadded_hash(both.iter()))
    }

    fn hash_elem_slice(&self, slice: &[BabyBearElem]) -> Box<Digest> {
        to_digest(unpadded_hash(slice.iter()))
    }

    fn hash_ext_elem_slice(&self, slice: &[BabyBearExtElem]) -> Box<Digest> {
        to_digest(unpadded_hash(
            slice.iter().flat_map(|ee| ee.subelems().iter()),
        ))
    }
}

struct Poseidon2RngFactory;

impl RngFactory<BabyBear> for Poseidon2RngFactory {
    fn new_rng(&self) -> Box<dyn Rng<BabyBear>> {
        Box::new(Poseidon2Rng::new())
    }
}

/// A hash suite using Poseidon2 for both MT hashes and RNG
pub struct Poseidon2HashSuite;

impl Poseidon2HashSuite {
    /// Construct a new Poseidon2HashSuite
    pub fn new_suite() -> HashSuite<BabyBear> {
        HashSuite {
            name: "poseidon2".into(),
            hashfn: Rc::new(Poseidon2HashFn {}),
            rng: Rc::new(Poseidon2RngFactory {}),
        }
    }
}

fn to_digest(elems: [BabyBearElem; CELLS_OUT]) -> Box<Digest> {
    let mut state: [u32; DIGEST_WORDS] = [0; DIGEST_WORDS];
    for i in 0..DIGEST_WORDS {
        state[i] = elems[i].as_u32_montgomery();
    }
    Box::new(Digest::from(state))
}

fn add_round_constants_full(cells: &mut [Elem; CELLS], round: usize) {
    for i in 0..CELLS {
        cells[i] += ROUND_CONSTANTS[round * CELLS + i];
    }
}

fn add_round_constants_partial(cells: &mut [Elem; CELLS], round: usize) {
    cells[0] += ROUND_CONSTANTS[round * CELLS];
}

fn sbox(x: Elem) -> Elem {
    let x2 = x * x;
    let x4 = x2 * x2;
    let x6 = x4 * x2;
    x6 * x
}

fn do_full_sboxes(cells: &mut [Elem; CELLS]) {
    for cell in cells.iter_mut() {
        *cell = sbox(*cell);
    }
}

fn do_partial_sboxes(cells: &mut [Elem; CELLS]) {
    cells[0] = sbox(cells[0]);
}

fn multiply_by_m_int(cells: &mut [Elem; CELLS]) {
    // Exploit the fact that off-diagonal entries of M_INT are all 1.
    let mut sum = Elem::new(0);
    for i in 0..CELLS {
        sum += cells[i];
    }
    for i in 0..CELLS {
        cells[i] = sum + M_INT_DIAG[i] * cells[i];
    }
}

fn multiply_by_m_ext(cells: &mut [Elem; CELLS]) {
    let old_cells = *cells;
    for i in 0..CELLS {
        let mut tot = Elem::new(0);
        for j in 0..CELLS {
            tot += M_EXT[i * CELLS + j] * old_cells[j];
        }
        cells[i] = tot;
    }
}

fn full_round(cells: &mut [Elem; CELLS], round: usize) {
    add_round_constants_full(cells, round);
    do_full_sboxes(cells);
    multiply_by_m_ext(cells);
}

fn partial_round(cells: &mut [Elem; CELLS], round: usize) {
    add_round_constants_partial(cells, round);
    do_partial_sboxes(cells);
    multiply_by_m_int(cells);
}

/// The raw sponge mixing function
pub fn poseidon2_mix(cells: &mut [Elem; CELLS]) {
    let mut round = 0;

    // First linear layer.
    multiply_by_m_ext(cells);

    // Do initial full rounds
    for _i in 0..ROUNDS_HALF_FULL {
        full_round(cells, round);
        round += 1;
    }
    // Do partial rounds
    for _i in 0..ROUNDS_PARTIAL {
        partial_round(cells, round);
        round += 1;
    }
    // Do remaining full rounds
    for _i in 0..ROUNDS_HALF_FULL {
        full_round(cells, round);
        round += 1;
    }
}

/// Perform a unpadded hash of a vector of elements.  Because this is unpadded
/// collision resistance is only true for vectors of the same size.  If the size
/// is variable, this is subject to length extension attacks.
pub fn unpadded_hash<'a, I>(iter: I) -> [Elem; CELLS_OUT]
where
    I: Iterator<Item = &'a Elem>,
{
    let mut state = [Elem::new(0); CELLS];
    let mut count = 0;
    let mut unmixed = 0;
    for val in iter {
        state[unmixed] += *val;
        count += 1;
        unmixed += 1;
        if unmixed == CELLS_RATE {
            poseidon2_mix(&mut state);
            unmixed = 0;
        }
    }
    if unmixed != 0 || count == 0 {
        poseidon2_mix(&mut state);
    }
    state.as_slice()[0..CELLS_OUT].try_into().unwrap()
}

#[cfg(test)]
mod tests {
    use test_log::test;

    use super::*;

    fn do_partial_sboxes(cells: &mut [Elem; CELLS]) {
        cells[0] = sbox(cells[0]);
    }

    fn partial_round_naive(cells: &mut [Elem; CELLS], round: usize) {
        add_round_constants_partial(cells, round);
        do_partial_sboxes(cells);
        multiply_by_m_int_naive(cells);
    }

    fn multiply_by_m_int_naive(cells: &mut [Elem; CELLS]) {
        let old_cells = *cells;
        for i in 0..CELLS {
            let mut tot = Elem::new(0);
            for j in 0..CELLS {
                if i == j {
                    tot += (M_INT_DIAG[i] + Elem::new(1)) * old_cells[j];
                } else {
                    tot += old_cells[j];
                }
            }
            cells[i] = tot;
        }
    }

    // Naive version of poseidon
    fn poseidon2_mix_naive(cells: &mut [Elem; CELLS]) {
        let mut round = 0;
        multiply_by_m_ext(cells);
        for _i in 0..ROUNDS_HALF_FULL {
            full_round(cells, round);
            round += 1;
        }
        for _i in 0..ROUNDS_PARTIAL {
            partial_round_naive(cells, round);
            round += 1;
        }
        for _i in 0..ROUNDS_HALF_FULL {
            full_round(cells, round);
            round += 1;
        }
    }

    #[test]
    fn compare_naive() {
        // Make a fixed input
        let mut test_in_1 = [Elem::new(1); CELLS];
        // Copy it
        let mut test_in_2 = test_in_1;
        // Try two versions
        poseidon2_mix_naive(&mut test_in_1);
        poseidon2_mix(&mut test_in_2);
        log::debug!("test_in_1: {:?}", test_in_1);
        log::debug!("test_in_2: {:?}", test_in_2);
        // Verify they are the same
        assert_eq!(test_in_1, test_in_2);
    }

    macro_rules! baby_bear_array {
        [$($x:literal),* $(,)?] => {
            [$(Elem::new($x)),* ]
        }
    }

    #[test]
    fn poseidon2_test_vectors() {
        let mut buf: &mut [Elem; CELLS] = &mut baby_bear_array![
            0x00000000, 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007,
            0x00000008, 0x00000009, 0x0000000A, 0x0000000B, 0x0000000C, 0x0000000D, 0x0000000E, 0x0000000F,
            0x00000010, 0x00000011, 0x00000012, 0x00000013, 0x00000014, 0x00000015, 0x00000016, 0x00000017,
        ];
        log::debug!("input: {:?}", buf);
        poseidon2_mix(&mut buf);
        let goal: [u32; CELLS] = [
            0x2ED3E23D, 0x12921FB0, 0x0E659E79, 0x61D81DC9, 0x32BAE33B, 0x62486AE3, 0x1E681B60, 0x24B91325,
            0x2A2EF5B9, 0x50E8593E, 0x5BC818EC, 0x10691997, 0x35A14520, 0x2BA6A3C5, 0x279D47EC, 0x55014E81,
            0x5953A67F, 0x2F403111, 0x6B8828FF, 0x1801301F, 0x2749207A, 0x3DC9CF21, 0x3C985BA2, 0x57A99864,
        ];
        for i in 0..CELLS {
            assert_eq!(buf[i].as_u32(), goal[i]);
        }

        log::debug!("output: {:?}", buf);
    }
}
