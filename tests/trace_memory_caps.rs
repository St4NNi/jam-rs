use jam_rs::trace::memory::{HeapRank, PerContigTopK, retain_per_contig_top_k};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Eq, PartialEq)]
struct Occurrence {
    contig: u32,
    id: u64,
    score: i64,
}

#[test]
fn per_contig_heaps_apply_caps_during_streaming_decode() {
    let occurrences = (0_u64..100)
        .map(|id| Occurrence {
            contig: (id % 3) as u32,
            id,
            score: id as i64,
        })
        .collect::<Vec<_>>();
    let selected = retain_per_contig_top_k(
        occurrences,
        2,
        3,
        |occurrence| occurrence.contig,
        |occurrence| HeapRank::new(occurrence.score, occurrence.id),
    );
    assert_eq!(
        selected.values().map(Vec::len).collect::<Vec<_>>(),
        vec![2, 2, 2]
    );
    assert_eq!(selected[&0][0].id, 99);
    assert_eq!(selected[&0][1].id, 96);
}

#[test]
fn unseen_contigs_are_bounded_deterministically() {
    let mut heaps = PerContigTopK::new(2, 2);
    assert!(heaps.push(10, 10_u64, HeapRank::new(10, 10)));
    assert!(heaps.push(20, 20_u64, HeapRank::new(20, 20)));
    assert!(!heaps.push(30, 30_u64, HeapRank::new(30, 30)));
    let selected = heaps.into_sorted();
    assert_eq!(selected, BTreeMap::from([(10, vec![10]), (20, vec![20])]));
}
