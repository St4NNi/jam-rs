use super::ChainConfig;
use super::Strand;
use super::chain_tests::random_sequence;
use super::{
    AlignmentConfig, AlignmentWorkspace, AnchorConfig, BandConfig, DenseOptimum, JoinContext,
    JoinLimits, JoinScratch, QueryTable, Walk, align_dense_walk, dense_walks, format_cigar,
    holds_walk, walk_alignment,
};

fn walks_of(query: &[u8], window: &[u8], circular: bool) -> (Vec<Walk>, Vec<u8>, QueryTable) {
    let table = QueryTable::new(
        query,
        circular,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let dense = dense_walks(&table, window, 0, BandConfig::default(), 64 << 20).unwrap();
    let (walks, reverse) = (dense.walks, dense.reverse);
    (walks, reverse, table)
}

#[test]
fn imperfect_periodic_pair_walks_every_diagonal_once() {
    // 10 kb of ACGT with one N, against itself: every diagonal a multiple of four apart is one
    // ungapped alignment across both N positions, on both strands.
    let mut query = b"ACGT".repeat(2500);
    query[5000] = b'N';
    let (walks, reverse, table) = walks_of(&query, &query, false);
    let mut lines = std::collections::BTreeSet::new();
    for walk in &walks {
        let diagonal = i64::from(walk.contig) - i64::from(walk.query);
        assert!(
            lines.insert((walk.strand.is_reverse(), diagonal)),
            "one walk per diagonal"
        );
        assert_eq!(diagonal.rem_euclid(4), 0);
        // The whole overlap, but a local alignment leaves out an N in its first or last column.
        let overlap = 10_000 - diagonal.unsigned_abs();
        assert!(overlap - u64::from(walk.length) <= 2, "{walk:?}");
        let oriented = if walk.strand.is_reverse() {
            &reverse
        } else {
            &query
        };
        let alignment = walk_alignment(walk, oriented, 0, &table, BandConfig::default()).unwrap();
        alignment.validate().unwrap();
        assert!(alignment.counts().substitutions <= 2);
        let cigar = format_cigar(&alignment.runs);
        assert!(cigar.ends_with('=') && !cigar.starts_with("1X"), "{cigar}");
        if (walk.strand, diagonal) == (Strand::Forward, 0) {
            assert_eq!(cigar, "10000=");
        }
    }
    let chained = walks.iter().filter(|walk| walk.length >= 40).count();
    assert_eq!(chained, 2 * (2 * (10_000 - 40) / 4 + 1));
}

#[test]
fn circular_walks_cover_the_query_at_most_once() {
    let mut query = b"ACGT".repeat(500);
    query[1000] = b'N';
    let window = [query.clone(), query.clone()].concat();
    let (walks, reverse, table) = walks_of(&query, &window, true);
    assert!(!walks.is_empty());
    assert!(walks.iter().all(|walk| walk.length <= 2000));
    assert!(walks.iter().any(|walk| walk.length == 2000));
    for walk in &walks {
        let oriented = if walk.strand.is_reverse() {
            &reverse
        } else {
            &window
        };
        let alignment = walk_alignment(walk, oriented, 0, &table, BandConfig::default()).unwrap();
        alignment.validate().unwrap();
        assert_eq!(alignment.query_interval.len(), u64::from(walk.length));
        assert!(holds_walk(&alignment, walk, window.len(), 0, 2000));
    }
}

#[test]
fn walks_end_where_the_copy_ends_and_skip_repeated_starts() {
    // A repeat copy inside unrelated sequence: the walk stops at the copy's end, not after the
    // X-drop tail, and starts inside it give no second walk.
    let element = b"ACG".repeat(400);
    let query = [b"GATTACATTGCCAGTTCAGG".to_vec(), element.clone()].concat();
    let window = [random_sequence(3000, 7), element, random_sequence(3000, 8)].concat();
    let (walks, _, _) = walks_of(&query, &window, false);
    let longest = walks
        .iter()
        .filter(|walk| !walk.strand.is_reverse())
        .max_by_key(|walk| walk.length)
        .unwrap();
    assert!((1200..1210).contains(&longest.length), "{longest:?}");
    let diagonal = i64::from(longest.contig) - i64::from(longest.query);
    let same = walks.iter().filter(|walk| {
        !walk.strand.is_reverse() && i64::from(walk.contig) - i64::from(walk.query) == diagonal
    });
    assert_eq!(same.count(), 1);
}

#[test]
fn match_starts_exceeding_the_budget_fail_before_allocation() {
    let query = b"ACGT".repeat(2500);
    let table = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    assert!(dense_walks(&table, &query, 0, BandConfig::default(), 1 << 10).is_err());
}

fn linear(query: &[u8]) -> QueryTable {
    QueryTable::new(
        query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap()
}

#[test]
fn dense_optimum_is_reused_only_for_matching_inputs_and_bands() {
    let query = [vec![b'A'; 100], vec![b'C'; 100]].concat();
    let contig = [vec![b'A'; 100], vec![b'G'; 20], vec![b'C'; 100]].concat();
    let table = linear(&query);
    let mut workspace = AlignmentWorkspace::default();
    let mut saved: Option<DenseOptimum> = None;
    let mut align = |walk: &Walk, config, saved: &mut Option<DenseOptimum>| {
        align_dense_walk(walk, &contig, 0, &table, config, &mut workspace, saved).unwrap()
    };
    let config = BandConfig::default();
    let walk = Walk {
        strand: Strand::Forward,
        contig: 0,
        query: 0,
        length: 200,
    };
    let first = align(&walk, config, &mut saved).unwrap();
    assert_eq!(format_cigar(&first.runs), "100=20D100=");
    assert!(saved.is_some());
    let other = Walk { contig: 4, ..walk };
    assert!(align(&other, config, &mut saved).is_none());
    assert_eq!(Some(first), align(&other, config, &mut None));
    let reverse = Walk {
        strand: Strand::Reverse,
        ..walk
    };
    assert!(align(&reverse, config, &mut saved).is_some());
    let narrow = BandConfig {
        max_drift: 1,
        alignment: AlignmentConfig {
            band_width: 1,
            ..config.alignment
        },
    };
    let shifted = Walk {
        contig: 100,
        ..walk
    };
    assert!(align(&shifted, narrow, &mut saved).is_some());
}

#[test]
fn dense_walk_of_a_long_query_aligns_the_bases_around_it() {
    let query = random_sequence(10_000, 7);
    let mut contig = random_sequence(300, 8);
    contig.extend_from_slice(&query[3_000..5_000]);
    contig.extend(random_sequence(20, 9));
    contig.extend_from_slice(&query[5_000..7_000]);
    contig.extend(random_sequence(300, 10));
    let walk = Walk {
        strand: Strand::Forward,
        contig: 1_300,
        query: 4_000,
        length: 300,
    };
    // The first window, 300 bases and twice the margin on each side, holds only a part of it.
    let config = BandConfig {
        max_drift: 128,
        ..BandConfig::default()
    };
    let table = linear(&query);
    let mut workspace = AlignmentWorkspace::default();
    let found = align_dense_walk(&walk, &contig, 0, &table, config, &mut workspace, &mut None)
        .unwrap()
        .unwrap();
    assert!(
        found.counts().deletions >= 20,
        "{}",
        format_cigar(&found.runs)
    );
    let aligned = found.query_interval;
    assert!(
        aligned.start <= 3_008 && aligned.end >= 6_992,
        "{aligned:?}"
    );
}

#[test]
fn short_blocks_join_across_a_divergent_piece() {
    // Two 30-base blocks of the query with a 6-base piece between them that differs: each block
    // is shorter than the 50 columns a join must reach, but together they reach it.
    let query = random_sequence(400, 31);
    let mut contig = query[100..166].to_vec();
    for base in &mut contig[30..36] {
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    let table = linear(&query);
    let config = BandConfig::default();
    let walks = [
        Walk {
            strand: Strand::Forward,
            contig: 0,
            query: 100,
            length: 30,
        },
        Walk {
            strand: Strand::Forward,
            contig: 36,
            query: 136,
            length: 30,
        },
    ];
    let context = JoinContext {
        table: &table,
        band: config,
        limits: JoinLimits {
            min_columns: 50,
            min_identity: 0.8,
        },
    };
    let mut scratch = JoinScratch::default();
    scratch.prepare(&walks, &context);
    let joined = scratch
        .join_through(&walks, 0, &contig, 1_000, &context)
        .unwrap()
        .unwrap();
    joined.validate().unwrap();
    assert_eq!(joined.query_interval.start, 100);
    assert_eq!(joined.query_interval.end, 166);
    assert_eq!(joined.contig_interval.start, 1_000);
    assert!(JoinScratch::bytes(2, &context) > 0);
}
