//! A source whose S3 object is checked once per process, and that fails after a change.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::etags::{Start, changed_error};
use super::{ReadAt, ReadCost, S3Object};
use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// A source that reads `object` through other layers (pins, the block cache). Those can answer
/// without a request, so the first read also checks the object with one byte, beside its own
/// work. After any object was found replaced, every read fails: open the shard again.
pub struct Checked {
    inner: Arc<dyn ReadAt>,
    object: Arc<S3Object>,
    changes: u64,
    passed: AtomicBool,
}

impl Checked {
    /// Reads `inner`, which reads `object`.
    pub fn new(inner: Arc<dyn ReadAt>, object: Arc<S3Object>) -> Self {
        let changes = object.changes();
        Self {
            inner,
            object,
            changes,
            passed: AtomicBool::new(false),
        }
    }
}

impl ReadAt for Checked {
    fn read_exact_at(&self, offset: u64, buffer: &mut [u8]) -> io::Result<()> {
        if self.object.changes() != self.changes {
            return Err(changed_error());
        }
        if self.passed.load(Ordering::Acquire) {
            return self.inner.read_exact_at(offset, buffer);
        }
        let latch = match self.object.begin_check() {
            Start::Passed => {
                self.passed.store(true, Ordering::Release);
                return self.inner.read_exact_at(offset, buffer);
            }
            Start::Run(latch) => {
                self.object.run_check(&latch);
                latch
            }
            Start::Wait(latch) => {
                // Without a known ETag the block keys are not known yet.
                if !self.object.is_known() {
                    latch.wait()?;
                }
                latch
            }
        };
        let read = self.inner.read_exact_at(offset, buffer);
        latch.wait()?;
        self.passed.store(true, Ordering::Release);
        read
    }

    fn size(&self) -> u64 {
        self.inner.size()
    }

    fn version(&self) -> io::Result<u64> {
        self.inner.version()
    }

    fn will_need(&self, offset: u64, length: u64) {
        self.inner.will_need(offset, length);
    }

    fn read_cost(&self) -> ReadCost {
        self.inner.read_cost()
    }
}
