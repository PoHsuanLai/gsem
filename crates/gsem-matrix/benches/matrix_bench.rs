use criterion::{Criterion, criterion_group, criterion_main};
use faer::Mat;
use std::hint::black_box;

fn bench_near_pd(c: &mut Criterion) {
    // 10x10 non-PD matrix
    let mut mat = Mat::from_fn(10, 10, |i, j| {
        if i == j {
            1.0
        } else {
            0.9 - 0.1 * (i as f64 - j as f64).abs()
        }
    });
    // Make it non-PD by forcing a negative eigenvalue structure
    mat[(0, 9)] = 0.99;
    mat[(9, 0)] = 0.99;

    c.bench_function("near_pd_10x10", |b| {
        b.iter(|| gsem_matrix::near_pd::nearest_pd(black_box(&mat), true, 100, 1e-8))
    });

    // 50x50
    let mat50 = Mat::from_fn(50, 50, |i, j| {
        if i == j {
            1.0
        } else {
            0.5 * (-(((i as f64 - j as f64) / 10.0).powi(2))).exp()
        }
    });
    c.bench_function("near_pd_50x50", |b| {
        b.iter(|| gsem_matrix::near_pd::nearest_pd(black_box(&mat50), false, 100, 1e-8))
    });
}

fn bench_vech(c: &mut Criterion) {
    let mat = Mat::from_fn(20, 20, |i, j| (i + j) as f64 * 0.1);
    c.bench_function("vech_20x20", |b| {
        b.iter(|| gsem_matrix::vech::vech(black_box(&mat)).unwrap())
    });

    let v = gsem_matrix::vech::vech(&mat).unwrap();
    c.bench_function("vech_reverse_20x20", |b| {
        b.iter(|| gsem_matrix::vech::vech_reverse(black_box(&v), 20).unwrap())
    });
}

fn bench_smooth(c: &mut Criterion) {
    let mat = Mat::from_fn(20, 20, |i, j| {
        if i == j {
            1.0
        } else {
            0.3 / (1.0 + (i as f64 - j as f64).abs())
        }
    });
    c.bench_function("is_pd_20x20", |b| {
        b.iter(|| gsem_matrix::smooth::is_pd(black_box(&mat)))
    });

    c.bench_function("cov_to_cor_20x20", |b| {
        let cov = Mat::from_fn(20, 20, |i, j| if i == j { (i + 1) as f64 } else { 0.5 });
        b.iter(|| gsem_matrix::smooth::cov_to_cor(black_box(&cov)))
    });
}

criterion_group!(benches, bench_near_pd, bench_vech, bench_smooth);
criterion_main!(benches);
