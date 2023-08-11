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

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use risc0_core::field::{baby_bear::BabyBearElem, Elem};
use risc0_zkp::core::hash::poseidon::{poseidon_mix, CELLS};

fn init_cells(n: usize) -> Vec<[BabyBearElem; CELLS]> {
    let mut vec = Vec::new();
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        vec.push([BabyBearElem::random(&mut rng); CELLS])
    }
    vec
}

fn benchmark_poseidon_mix(c: &mut Criterion) {
    let mut group = c.benchmark_group("poseidon_mix");
    for n in [1024, 4096, 8192].iter() {
        let mut cells = init_cells(*n);
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            b.iter(|| {
                for i in 0..cells.len() {
                    poseidon_mix(&mut cells[i])
                }
            })
        });
    }

    group.finish();
}

criterion_group!(benches, benchmark_poseidon_mix);
criterion_main!(benches);

