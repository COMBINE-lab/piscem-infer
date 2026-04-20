use std::env;
use std::fs::File;
use std::io::BufReader;

use arrow2::array::{ListArray, PrimitiveArray};
use arrow2::io::parquet::read::{self as pq_read, FileReader};

fn main() {
    let path = env::args().nth(1).expect("need path arg");
    let file = File::open(&path).unwrap();
    let reader = BufReader::new(file);
    let metadata = pq_read::read_metadata(&mut BufReader::new(File::open(&path).unwrap())).unwrap();
    let schema = pq_read::infer_schema(&metadata).unwrap();
    let mut file_reader = FileReader::new(reader, metadata.row_groups, schema, None, None, None);
    let chunk = file_reader.next().unwrap().unwrap();
    let columns = chunk.columns();

    let labels_array = columns[0]
        .as_any()
        .downcast_ref::<ListArray<i32>>()
        .unwrap();
    let values = labels_array
        .values()
        .as_any()
        .downcast_ref::<PrimitiveArray<u32>>()
        .unwrap();
    let counts_array = columns[1]
        .as_any()
        .downcast_ref::<PrimitiveArray<u64>>()
        .unwrap();

    let n = labels_array.len();
    println!("Total EQCs: {}", n);

    println!("\nFirst 15 EQ classes:");
    for i in 0..std::cmp::min(15, n) {
        let start = labels_array.offsets()[i] as usize;
        let end = labels_array.offsets()[i + 1] as usize;
        let labels: Vec<u32> = values
            .values_iter()
            .skip(start)
            .take(end - start)
            .copied()
            .collect();
        let count = counts_array.value(i);
        println!("  EQC {}: labels={:?}, count={}", i, labels, count);
    }

    let mut len_counts = std::collections::HashMap::new();
    // Count duplicate target sets (sort labels to normalize order)
    let mut target_set_counts: std::collections::HashMap<Vec<u32>, Vec<u64>> =
        std::collections::HashMap::new();
    for i in 0..n {
        let start = labels_array.offsets()[i] as usize;
        let end = labels_array.offsets()[i + 1] as usize;
        let len = end - start;
        *len_counts.entry(len).or_insert(0usize) += 1;

        let mut labels: Vec<u32> = values
            .values_iter()
            .skip(start)
            .take(len)
            .copied()
            .collect();
        labels.sort();
        let count = counts_array.value(i);
        target_set_counts.entry(labels).or_default().push(count);
    }
    let mut lens: Vec<_> = len_counts.into_iter().collect();
    lens.sort();
    println!("\nLabel length distribution:");
    for (len, count) in &lens {
        println!("  len={}: {} EQCs", len, count);
    }

    let n_unique_sets = target_set_counts.len();
    let n_duplicated = target_set_counts.values().filter(|v| v.len() > 1).count();
    let n_dup_entries: usize = target_set_counts
        .values()
        .filter(|v| v.len() > 1)
        .map(|v| v.len())
        .sum();
    println!(
        "\nDuplicate analysis: {} unique target sets, {} sets appear >1 time ({} extra entries)",
        n_unique_sets,
        n_duplicated,
        n_dup_entries - n_duplicated
    );

    if n_duplicated > 0 {
        println!("Example duplicated target sets:");
        let mut dups: Vec<_> = target_set_counts
            .iter()
            .filter(|(_, v)| v.len() > 1)
            .collect();
        dups.sort_by_key(|(k, _)| k.clone());
        for (targets, counts) in dups.iter().take(10) {
            println!("  targets={:?}, counts={:?}", targets, counts);
        }
    }
}
