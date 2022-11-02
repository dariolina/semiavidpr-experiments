use semiavidpr::SemiAvidPr;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

fn bench_disperse_compute_column_commitments(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let mut group = c.benchmark_group("disperse_compute_column_commitments".to_string());
    group.sample_size(10);

    let coded_chunks = 1024;
    for chunk_length in vec![256, 512, 1024, 2048] {
        let uncoded_chunks = coded_chunks / 3;
        let scheme = SemiAvidPr::setup(&mut rng, coded_chunks, uncoded_chunks, chunk_length);

        group.throughput(Throughput::Bytes(scheme.get_filesize_in_bytes()));
        group.bench_with_input(
            format!("n={} L={}", coded_chunks, chunk_length),
            &(coded_chunks, chunk_length),
            |b, (_coded_chunks, _chunk_length)| {
                let file_uncoded = scheme.generate_random_file(&mut rng);
                b.iter(|| {
                    black_box(scheme.disperse_compute_column_commitments(&file_uncoded));
                })
            },
        );
    }

    group.finish();
}

fn bench_disperse_encode_rows(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let mut group = c.benchmark_group("disperse_encode_rows".to_string());
    group.sample_size(10);

    for coded_chunks in vec![256, 512, 1024, 2048] {
        for chunk_length in vec![256, 512, 1024, 2048] {
            let uncoded_chunks = coded_chunks / 3;
            let scheme = SemiAvidPr::setup(&mut rng, coded_chunks, uncoded_chunks, chunk_length);

            group.throughput(Throughput::Bytes(scheme.get_filesize_in_bytes()));
            group.bench_with_input(
                format!("n={} L={}", coded_chunks, chunk_length),
                &(coded_chunks, chunk_length),
                |b, (_coded_chunks, _chunk_length)| {
                    let file_uncoded = scheme.generate_random_file(&mut rng);
                    b.iter(|| {
                        black_box(scheme.disperse_encode_rows_lagrange(&file_uncoded));
                    })
                },
            );
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
