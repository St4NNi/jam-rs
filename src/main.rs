//! The jam binary: runs the command line and prints an error in one line.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use std::process::ExitCode;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() -> ExitCode {
    match jam_rs::run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("jam: {error}");
            ExitCode::FAILURE
        }
    }
}
