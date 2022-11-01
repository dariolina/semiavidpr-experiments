use semiavidpr::SemiAvidPr;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

fn bench_disperse_compute_column_commitments(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let mut group = c.benchmark_group("disperse_compute_column_commitments".to_string());
    group.sample_size(10);

    let n = 1024;
    #[allow(non_snake_case)]
    for L in vec![256, 512, 1024, 2048] {
        let k = n / 3;
        let scheme = SemiAvidPr::setup(&mut rng, n, k, L);

        group.throughput(Throughput::Bytes(scheme.get_filesize_in_bytes()));
        group.bench_with_input(format!("n={} L={}", n, L), &(n, L), |b, (_n, _L)| {
            let file_uncoded = scheme.generate_random_file(&mut rng);
            b.iter(|| {
                black_box(scheme.disperse_compute_column_commitments(&file_uncoded));
            })
        });
    }

    group.finish();
}

fn bench_disperse_encode_rows(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let mut group = c.benchmark_group("disperse_encode_rows".to_string());
    group.sample_size(10);

    for n in vec![256, 512, 1024, 2048] {
        #[allow(non_snake_case)]
        for L in vec![256, 512, 1024, 2048] {
            let k = n / 3;
            let scheme = SemiAvidPr::setup(&mut rng, n, k, L);

            group.throughput(Throughput::Bytes(scheme.get_filesize_in_bytes()));
            group.bench_with_input(format!("n={} L={}", n, L), &(n, L), |b, (_n, _L)| {
                let file_uncoded = scheme.generate_random_file(&mut rng);
                b.iter(|| {
                    black_box(scheme.disperse_encode_rows_lagrange(&file_uncoded));
                })
            });
        }
    }

    group.finish();
}

fn bench_bls12_381(c: &mut Criterion) {
    bench_disperse_compute_column_commitments(c);
    bench_disperse_encode_rows(c);
}

criterion_group!(benches, bench_bls12_381);
criterion_main!(benches);
