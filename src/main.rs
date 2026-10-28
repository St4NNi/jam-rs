//! The jam binary: runs the command line and prints an error in one line.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::process::ExitCode;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

/// Large allocations stay in the usual arenas, which give pages back to the kernel after a
/// delay, so that buffers freed and allocated again for every batch do not fault in their pages
/// each time.
#[cfg(not(target_env = "msvc"))]
#[unsafe(export_name = "_rjem_malloc_conf")]
static MALLOC_CONF: &[u8; 21] = b"oversize_threshold:0\0";

fn main() -> ExitCode {
    match jam_rs::run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("jam: {error}");
            ExitCode::FAILURE
        }
    }
}
