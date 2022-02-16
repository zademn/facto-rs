use criterion::{black_box, criterion_group, criterion_main, Criterion};
use facto_rs::{Factorizer, QuadraticSieveBuilder};
use rug::Integer;

fn qs_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("qs");
    let ps = [
        Integer::from_str_radix("54941", 10).unwrap(),
        Integer::from_str_radix("2410055309", 10).unwrap(),
        Integer::from_str_radix("9270271495854497069", 10).unwrap(),
        Integer::from_str_radix("629449496776611089395111", 10).unwrap(),
    ];
    let qs = [
        Integer::from_str_radix("58271", 10).unwrap(),
        Integer::from_str_radix("3406603901", 10).unwrap(),
        Integer::from_str_radix("15658485116535926689", 10).unwrap(),
        Integer::from_str_radix("1157146192290239025720557", 10).unwrap(),
    ];
    let ns = [
        Integer::from_str_radix("3201467011", 10).unwrap(),
        Integer::from_str_radix("8210103817265160409", 10).unwrap(),
        Integer::from_str_radix("145158408244084883965506502843949374541", 10).unwrap(),
        Integer::from_str_radix("728365088434062605443262763506401413148047996827", 10).unwrap(),
    ];

    let bits = [32, 64, 128, 160];
    for (((p, q), n), bit) in ps.iter().zip(qs).zip(ns).zip(bits) {
        let qsieve_builder = QuadraticSieveBuilder::new(n).verbose(false);
        let qsieve = qsieve_builder.build();
        group
            .sample_size(10)
            .bench_function(format!("qs_{bit}"), |b| {
                b.iter(|| {
                    let res = qsieve.factor();
                    assert_eq!(res, Some((p.clone(), q.clone())));
                })
            });
    }
}

fn qs_threaded_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("qs");
    let ps = [
        Integer::from_str_radix("54941", 10).unwrap(),
        Integer::from_str_radix("2410055309", 10).unwrap(),
        Integer::from_str_radix("9270271495854497069", 10).unwrap(),
        Integer::from_str_radix("629449496776611089395111", 10).unwrap(),
    ];
    let qs = [
        Integer::from_str_radix("58271", 10).unwrap(),
        Integer::from_str_radix("3406603901", 10).unwrap(),
        Integer::from_str_radix("15658485116535926689", 10).unwrap(),
        Integer::from_str_radix("1157146192290239025720557", 10).unwrap(),
    ];
    let ns = [
        Integer::from_str_radix("3201467011", 10).unwrap(),
        Integer::from_str_radix("8210103817265160409", 10).unwrap(),
        Integer::from_str_radix("145158408244084883965506502843949374541", 10).unwrap(),
        Integer::from_str_radix("728365088434062605443262763506401413148047996827", 10).unwrap(),
    ];
    let bits = [32, 64, 128, 160];

    let num_cores = 4;
    for (((p, q), n), bit) in ps.iter().zip(qs).zip(ns).zip(bits) {
        let qsieve_builder = QuadraticSieveBuilder::new(n)
            .num_cores(num_cores)
            .sieve_size(100000)
            .verbose(false);
        let qsieve = qsieve_builder.build();
        group
            .sample_size(10)
            .bench_function(format!("qs_threaded_{bit}"), |b| {
                b.iter(|| {
                    let res = qsieve.factor();
                    assert_eq!(res, Some((p.clone(), q.clone())));
                })
            });
    }
}
criterion_group!(benches, qs_benchmark, qs_threaded_benchmark);
criterion_main!(benches);
