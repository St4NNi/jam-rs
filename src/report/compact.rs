//! The compact JSON lines of `jam search`: one line per reported entry at the level of the
//! search, and one summary line per query over all shards.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::{Confirmed, FragmentRecord, Reason, Record, Status, Summary};
use crate::align::Strand;
use crate::search::SearchStage;
use serde::Serialize;

/// A fragment: contig, query start and end, contig start and end, strand and CIGAR.
type FragmentLine<'a> = (&'a str, u64, u64, u64, u64, &'static str, &'a str);

/// A result line. `s` is the stage the entry reached: 1 after sketch, 2 located, 3 confirmed.
#[derive(Serialize)]
struct ResultLine<'a> {
    q: &'a str,
    a: &'a str,
    s: u8,
    /// The containment estimate (stages 1 and 2).
    #[serde(skip_serializing_if = "Option::is_none")]
    c: Option<f64>,
    /// The alignment windows: contig, start and end (stage 2).
    #[serde(skip_serializing_if = "Option::is_none")]
    w: Option<Vec<(&'a str, u64, u64)>>,
    #[serde(flatten, skip_serializing_if = "Option::is_none")]
    aligned: Option<AlignedLine<'a>>,
}

/// The part of a result line that only a confirmed entry has.
#[derive(Serialize)]
struct AlignedLine<'a> {
    cov: f64,
    id: f64,
    bases: u64,
    f: Vec<FragmentLine<'a>>,
    /// Fragments that add no new query base, such as other copies of the same region.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    alt: Vec<FragmentLine<'a>>,
}

/// The result line of `record` in a search up to `level`: possible entries at the sketch
/// level, located entries at the locate level, confirmed entries; `None` for other records.
pub fn result_line(record: &Record, level: SearchStage) -> serde_json::Result<Option<String>> {
    let line = match record {
        Record::Possible(possible) if level == SearchStage::Sketch => ResultLine {
            q: &possible.query,
            a: &possible.entry,
            s: 1,
            c: Some(rounded(possible.containment_estimate, 4)),
            w: None,
            aligned: None,
        },
        Record::Located(located) if level == SearchStage::Locate => ResultLine {
            q: &located.query,
            a: &located.entry,
            s: 2,
            c: Some(rounded(located.containment_estimate, 4)),
            w: Some(
                (located.windows.iter())
                    .map(|window| (window.contig.as_str(), window.start, window.end))
                    .collect(),
            ),
            aligned: None,
        },
        Record::Confirmed(confirmed) => confirmed_line(confirmed),
        _ => return Ok(None),
    };
    serde_json::to_string(&line).map(Some)
}

fn confirmed_line(confirmed: &Confirmed) -> ResultLine<'_> {
    let decision = &confirmed.decision;
    ResultLine {
        q: &confirmed.query,
        a: &confirmed.entry,
        s: 3,
        c: None,
        w: None,
        aligned: Some(AlignedLine {
            cov: rounded(decision.query_coverage, 4),
            id: rounded(decision.identity, 4),
            bases: decision.covered_bases,
            f: fragment_lines(&decision.fragments),
            alt: fragment_lines(&decision.alternatives),
        }),
    }
}

fn fragment_lines(fragments: &[FragmentRecord]) -> Vec<FragmentLine<'_>> {
    fragments.iter().map(fragment_line).collect()
}

/// The fragment as PAF places it: a fragment over the end of a circular query ends past the
/// query length.
fn fragment_line(fragment: &FragmentRecord) -> FragmentLine<'_> {
    let query_start = fragment.query_segments.first().map_or(0, |span| span.start);
    let strand = match fragment.strand {
        Strand::Forward => "+",
        Strand::Reverse => "-",
    };
    (
        &fragment.contig,
        query_start,
        query_start + fragment.query_bases(),
        fragment.contig_interval.start,
        fragment.contig_interval.end,
        strand,
        &fragment.cigar,
    )
}

/// `value` rounded to `places` decimal places, so the lines stay short.
fn rounded(value: f64, places: i32) -> f64 {
    let scale = 10f64.powi(places);
    (value * scale).round() / scale
}

/// The summary of one query over all shards.
#[derive(Clone, Debug, PartialEq)]
pub struct QueryTotals {
    query: String,
    status: Status,
    reason: Option<Reason>,
    error: Option<String>,
    /// Seeded, possible, located and confirmed entries.
    counts: [u64; 4],
    miss: f64,
    /// Seconds of the sketch, locate and align stages.
    seconds: [f64; 3],
}

/// The summary line of a query.
#[derive(Serialize)]
struct SummaryLine<'a> {
    q: &'a str,
    status: Status,
    #[serde(skip_serializing_if = "Option::is_none")]
    reason: Option<Reason>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<&'a str>,
    n: [u64; 4],
    miss: f64,
    t: [f64; 3],
}

impl QueryTotals {
    /// The totals of `query` before any shard reported.
    pub fn new(query: &str) -> Self {
        Self {
            query: query.to_string(),
            status: Status::Complete,
            reason: None,
            error: None,
            counts: [0; 4],
            miss: 0.0,
            seconds: [0.0; 3],
        }
    }

    /// Adds the summary of one shard; the worst status and its reason win.
    pub fn add(&mut self, summary: &Summary) {
        if summary.status > self.status {
            self.status = summary.status;
            self.reason = summary.reason;
            self.error.clone_from(&summary.error);
        }
        let located = summary.locate.entries - summary.locate.dropped_entries;
        let counts = [
            summary.sketch.seeded_entries,
            summary.sketch.possible_entries,
            located,
            summary.align.confirmed_entries,
        ];
        for (total, count) in self.counts.iter_mut().zip(counts) {
            *total += count;
        }
        self.miss = self.miss.max(summary.sketch.miss_upper_bound);
        let times = &summary.times;
        let sketch = times.prepare_query.wall_seconds + times.sketch.wall_seconds;
        let seconds = [sketch, summary.locate.seconds, summary.align.seconds];
        for (total, part) in self.seconds.iter_mut().zip(seconds) {
            *total += part;
        }
    }

    /// The status over all shards so far.
    pub fn status(&self) -> Status {
        self.status
    }

    /// The summary line.
    pub fn line(&self) -> serde_json::Result<String> {
        serde_json::to_string(&SummaryLine {
            q: &self.query,
            status: self.status,
            reason: self.reason,
            error: self.error.as_deref(),
            n: self.counts,
            miss: self.miss,
            t: self.seconds.map(|seconds| rounded(seconds, 3)),
        })
    }

    /// The progress line of `-v`: status, stage counts and seconds.
    pub fn progress(&self) -> String {
        let status = serde_json::to_value(self.status).unwrap_or_default();
        let status = status.as_str().unwrap_or_default().to_string();
        let status = match self.reason.and_then(|r| serde_json::to_value(r).ok()) {
            Some(reason) => format!("{status} ({})", reason.as_str().unwrap_or_default()),
            None => status,
        };
        let [seeded, possible, located, confirmed] = self.counts;
        let [sketch, locate, align] = self.seconds;
        format!(
            "{}: {status}; {seeded} seeded, {possible} possible, {located} located, \
             {confirmed} confirmed; sketch {sketch:.3} s, locate {locate:.3} s, align {align:.3} s",
            self.query
        )
    }
}

#[cfg(test)]
#[path = "compact_tests.rs"]
mod tests;
