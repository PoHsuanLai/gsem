use criterion::{Criterion, black_box, criterion_group, criterion_main};
use faer::Mat;

fn bench_parse_model(c: &mut Criterion) {
    let model = "F1 =~ NA*V1 + V2 + V3 + V4 + V5\n\
                  F2 =~ NA*V3 + V4 + V5 + V6 + V7\n\
                  F1 ~~ 1*F1\nF2 ~~ 1*F2\nF1 ~~ F2\n\
                  V1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3\nV4 ~~ V4\n\
                  V5 ~~ V5\nV6 ~~ V6\nV7 ~~ V7";

    c.bench_function("parse_model_2factor_7var", |b| {
        b.iter(|| gsem_sem::syntax::parse_model(black_box(model), false))
    });
}

fn bench_implied_cov(c: &mut Criterion) {
    let model_str = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3";
    let pt = gsem_sem::syntax::parse_model(model_str, false).unwrap();
    let obs: Vec<String> = (1..=3).map(|i| format!("V{i}")).collect();
    let model = gsem_sem::model::Model::from_partable(&pt, &obs);

    c.bench_function("implied_cov_1factor_3var", |b| {
        b.iter(|| black_box(&model).implied_cov())
    });
}

fn bench_fit_dwls(c: &mut Criterion) {
    let s = faer::mat![[0.60, 0.42, 0.35], [0.42, 0.50, 0.30], [0.35, 0.30, 0.40],];
    let v_diag = vec![0.001; 6];
    let model_str = "F1 =~ NA*V1 + V2 + V3\nF1 ~~ 1*F1\nV1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3";
    let pt = gsem_sem::syntax::parse_model(model_str, false).unwrap();
    let obs: Vec<String> = (1..=3).map(|i| format!("V{i}")).collect();

    c.bench_function("fit_dwls_1factor_3var", |b| {
        b.iter(|| {
            let mut model = gsem_sem::model::Model::from_partable(&pt, &obs);
            gsem_sem::estimator::fit_dwls(
                black_box(&mut model),
                black_box(&s),
                black_box(&v_diag),
                1000,
                None,
            )
        })
    });

    // Larger model: 2 factors, 6 vars
    let s6 = Mat::from_fn(6, 6, |i, j| {
        if i == j {
            0.5
        } else {
            0.2 / (1.0 + (i as f64 - j as f64).abs())
        }
    });
    let v_diag6 = vec![0.001; 21]; // 6*7/2
    let model6 = "F1 =~ NA*V1 + V2 + V3\nF2 =~ NA*V4 + V5 + V6\n\
                   F1 ~~ 1*F1\nF2 ~~ 1*F2\nF1 ~~ F2\n\
                   V1 ~~ V1\nV2 ~~ V2\nV3 ~~ V3\nV4 ~~ V4\nV5 ~~ V5\nV6 ~~ V6";
    let pt6 = gsem_sem::syntax::parse_model(model6, false).unwrap();
    let obs6: Vec<String> = (1..=6).map(|i| format!("V{i}")).collect();

    c.bench_function("fit_dwls_2factor_6var", |b| {
        b.iter(|| {
            let mut model = gsem_sem::model::Model::from_partable(&pt6, &obs6);
            gsem_sem::estimator::fit_dwls(
                black_box(&mut model),
                black_box(&s6),
                black_box(&v_diag6),
                1000,
                None,
            )
        })
    });
}

criterion_group!(
    benches,
    bench_parse_model,
    bench_implied_cov,
    bench_fit_dwls
);
criterion_main!(benches);
