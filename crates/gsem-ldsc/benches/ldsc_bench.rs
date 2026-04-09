use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_weights(c: &mut Criterion) {
    let n = 10_000;
    let chi2: Vec<f64> = (0..n).map(|i| 1.0 + 0.5 * (i as f64 / n as f64)).collect();
    let ld: Vec<f64> = (0..n).map(|i| 5.0 + 20.0 * (i as f64 / n as f64)).collect();
    let w_ld = ld.clone();
    let n_vec = vec![50000.0; n];

    c.bench_function("h2_weights_10k", |b| {
        b.iter(|| {
            gsem_ldsc::weights::compute_h2_weights(
                black_box(&chi2),
                black_box(&ld),
                black_box(&w_ld),
                black_box(&n_vec),
                1_000_000.0,
            )
        })
    });
}

fn bench_regression(c: &mut Criterion) {
    let n = 10_000;
    let ld: Vec<f64> = (0..n).map(|i| 5.0 + 20.0 * (i as f64 / n as f64)).collect();
    let y: Vec<f64> = (0..n)
        .map(|i| 1.0 + 0.01 * ld[i] + 0.1 * ((i * 7) as f64 % 3.0))
        .collect();
    let w = vec![1.0 / n as f64; n];

    c.bench_function("block_regression_10k_200blocks", |b| {
        b.iter(|| {
            gsem_ldsc::regression::block_regression(
                black_box(&ld),
                black_box(&y),
                black_box(&w),
                200,
            )
        })
    });
}

fn bench_h2_estimation(c: &mut Criterion) {
    let n = 50_000;
    let z: Vec<f64> = (0..n)
        .map(|i| {
            let x = (i as f64 * 2.71828).sin();
            x * 1.3
        })
        .collect();
    let n_vec = vec![50000.0; n];
    let ld: Vec<f64> = (0..n).map(|i| 5.0 + 30.0 * (i as f64 / n as f64)).collect();
    let w_ld = ld.clone();

    c.bench_function("estimate_h2_50k_snps", |b| {
        b.iter(|| {
            gsem_ldsc::heritability::estimate_h2(
                black_box(&z),
                black_box(&n_vec),
                black_box(&ld),
                black_box(&w_ld),
                1_000_000.0,
                200,
            )
        })
    });
}

criterion_group!(
    benches,
    bench_weights,
    bench_regression,
    bench_h2_estimation
);
criterion_main!(benches);
