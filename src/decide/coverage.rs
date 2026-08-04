//! The query coverage that an entry must reach: a share of the query or a number of bases.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::DecideError;
use std::str::FromStr;

/// The least query coverage that reports an entry.
///
/// Text forms: a share between 0 and 1 such as `0.5`, or a base count such as `800bp`. A share is
/// kept as an exact decimal, so `0.7` never turns into 0.69999.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QueryCoverage {
    /// `numerator / 10^decimal_places` of the query.
    Share { numerator: u64, decimal_places: u32 },
    /// A fixed number of query bases.
    Bases(u64),
}

impl FromStr for QueryCoverage {
    type Err = DecideError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        if let Some(count) = input.strip_suffix("bp") {
            let bases = count
                .parse()
                .map_err(|_| DecideError::InvalidCoverage("not a base count"))?;
            return Ok(Self::Bases(bases));
        }
        let value: f64 = input
            .parse()
            .map_err(|_| DecideError::InvalidCoverage("not a number"))?;
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(DecideError::InvalidCoverage("a share must be from 0 to 1"));
        }
        let decimal = value.to_string();
        let (whole, fraction) = decimal.split_once('.').unwrap_or((&decimal, ""));
        let digits = format!("{whole}{fraction}");
        let digits = digits.trim_start_matches('0').trim_end_matches('0');
        if digits.is_empty() {
            return Ok(Self::Share {
                numerator: 0,
                decimal_places: 0,
            });
        }
        let trailing = fraction.len() - fraction.trim_end_matches('0').len();
        Ok(Self::Share {
            numerator: digits
                .parse()
                .map_err(|_| DecideError::InvalidCoverage("too many digits"))?,
            decimal_places: (fraction.len() - trailing) as u32,
        })
    }
}

impl QueryCoverage {
    /// The number of query bases that must be covered. A base count above the query length is
    /// lowered to the query length.
    pub fn required_bases(self, query_length: u64) -> u64 {
        match self {
            Self::Bases(bases) => bases.min(query_length),
            Self::Share {
                numerator,
                decimal_places,
            } => {
                // A u64 times a u64 fits in u128. A share finer than u128 decimals still needs
                // one base when the query is not empty.
                let product = u128::from(query_length) * u128::from(numerator);
                match 10u128.checked_pow(decimal_places) {
                    Some(denominator) => product.div_ceil(denominator) as u64,
                    None => u64::from(product != 0),
                }
            }
        }
    }

    /// Checks that a share lies from 0.05 to 1, as a large plasmid can show as a few small
    /// pieces, and that a base count is not zero.
    pub fn validate(self) -> Result<(), DecideError> {
        // A length on which every decimal share of up to 17 digits is exact.
        const LENGTH: u64 = 1_000_000_000_000_000_000;
        match self {
            Self::Bases(0) => Err(DecideError::InvalidRule("query coverage is zero bases")),
            Self::Bases(_) => Ok(()),
            Self::Share { .. } if (LENGTH / 20..=LENGTH).contains(&self.required_bases(LENGTH)) => {
                Ok(())
            }
            Self::Share { .. } => Err(DecideError::InvalidRule("query coverage is below 0.05")),
        }
    }

    /// Whether `covered_bases` of a query of `query_length` bases are enough.
    pub fn accepts(self, covered_bases: u64, query_length: u64) -> bool {
        covered_bases >= self.required_bases(query_length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shares_round_up_to_whole_bases() {
        for (text, percent) in [
            ("0", 0),
            (".1", 10),
            (".25", 25),
            (".50", 50),
            (".70", 70),
            (".9", 90),
            ("1", 100),
        ] {
            let coverage: QueryCoverage = text.parse().unwrap();
            for length in [1u64, 99, 100, 101, 50_001, 200_002, u64::MAX] {
                let required = (u128::from(length) * percent).div_ceil(100) as u64;
                assert_eq!(coverage.required_bases(length), required);
                assert!(coverage.accepts(required, length));
                if required > 0 {
                    assert!(!coverage.accepts(required - 1, length));
                }
            }
        }
    }

    #[test]
    fn bad_text_is_rejected() {
        for text in [
            "NaN", "inf", "-inf", "-0.1", "1.01", "bad", "-5bp", "bp", "1.5bp",
        ] {
            assert!(text.parse::<QueryCoverage>().is_err(), "{text}");
        }
    }

    #[test]
    fn tiny_shares_need_one_base() {
        for text in ["5e-324", "1e-100", "1e-20"] {
            let coverage: QueryCoverage = text.parse().unwrap();
            assert_eq!(coverage.required_bases(u64::MAX), 1);
        }
        let coverage: QueryCoverage = "5e-1".parse().unwrap();
        assert_eq!(coverage.required_bases(101), 51);
    }

    #[test]
    fn base_counts_are_capped_by_the_query() {
        let coverage: QueryCoverage = "800bp".parse().unwrap();
        assert_eq!(coverage.required_bases(10_000), 800);
        assert_eq!(coverage.required_bases(500), 500);
        assert!(coverage.accepts(800, 10_000));
        assert!(!coverage.accepts(799, 10_000));
    }
}
