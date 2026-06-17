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
