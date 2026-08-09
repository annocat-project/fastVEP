use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{raw_json_blob_fields, Osa2Metadata, Osa2Record, Osa2Writer};

fn record(position: u32, ref_allele: &[u8], alt_allele: &[u8], row: u32) -> Osa2Record {
    Osa2Record {
        chrom: "1".into(),
        position,
        ref_allele: ref_allele.to_vec(),
        alt_allele: alt_allele.to_vec(),
        values: Vec::new(),
        json_blob: Some(format!(r#"{{"row":{row}}}"#)),
    }
}

fn json(value: AnnotationValue) -> String {
    match value {
        AnnotationValue::Json(json) => json,
        other => panic!("expected allele-specific JSON, got {other:?}"),
    }
}

fn write_fixture(path: &std::path::Path, record_list: bool) {
    let records = vec![
        record(10_000, b"A", b"G", 1),
        record(10_000, b"A", b"G", 2),
        record(10_100, b"ACGTA", b"A", 3),
        record(10_100, b"ACGTA", b"A", 4),
        record(10_200, b"N", b"A", 5),
        record(10_200, b"N", b"A", 6),
        record(10_300, b"ANNNN", b"A", 7),
        record(10_300, b"ANNNN", b"A", 8),
    ];
    let metadata = Osa2Metadata {
        format_version: 2,
        name: "dbNSFP".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: "dbnsfp".into(),
        match_by_allele: true,
        is_array: false,
        record_list,
        is_positional: false,
        chunk_bits: 20,
        description: "duplicate-key regression".into(),
    };

    Osa2Writer::new(metadata, raw_json_blob_fields())
        .write_all(std::fs::File::create(path).unwrap(), &records)
        .unwrap();
}

#[test]
fn record_list_returns_all_duplicate_compact_long_and_raw_records() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("duplicate_keys.osa2");
    write_fixture(&path, true);

    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(reader.verify(Some("1")).unwrap().record_count, 8);
    assert_eq!(
        json(
            reader
                .annotate_position("1", 10_000, "A", "G")
                .unwrap()
                .unwrap()
        ),
        r#"[{"row":1},{"row":2}]"#
    );
    assert_eq!(
        json(
            reader
                .annotate_position("1", 10_100, "ACGTA", "A")
                .unwrap()
                .unwrap()
        ),
        r#"[{"row":3},{"row":4}]"#
    );
    assert_eq!(
        json(
            reader
                .annotate_position("1", 10_200, "N", "A")
                .unwrap()
                .unwrap()
        ),
        r#"[{"row":5},{"row":6}]"#
    );
    assert_eq!(
        json(
            reader
                .annotate_position("1", 10_300, "ANNNN", "A")
                .unwrap()
                .unwrap()
        ),
        r#"[{"row":7},{"row":8}]"#
    );
}

#[test]
fn ordinary_lookup_keeps_first_record_behavior() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("single_record.osa2");
    write_fixture(&path, false);

    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(
        json(
            reader
                .annotate_position("1", 10_000, "A", "G")
                .unwrap()
                .unwrap()
        ),
        r#"{"row":1}"#
    );
}
