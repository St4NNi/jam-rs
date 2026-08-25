//! Process-local stage profiling for measured Jam batch runs.

use serde::Serialize;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant};

static ACTIVE: OnceLock<Mutex<Option<Arc<Profiler>>>> = OnceLock::new();

fn active_slot() -> &'static Mutex<Option<Arc<Profiler>>> {
    ACTIVE.get_or_init(|| Mutex::new(None))
}

#[derive(Clone, Copy, Debug, Default)]
struct Usage {
    cpu_micros: u64,
    minor_faults: u64,
    major_faults: u64,
}

impl Usage {
    fn for_who(who: libc::c_int) -> Self {
        let mut value = std::mem::MaybeUninit::<libc::rusage>::zeroed();
        // SAFETY: `value` points to writable storage for `getrusage` and is read
        // only when the system call succeeds.
        if unsafe { libc::getrusage(who, value.as_mut_ptr()) } != 0 {
            return Self::default();
        }
        // SAFETY: a successful `getrusage` initialized the complete structure.
        let value = unsafe { value.assume_init() };
        Self {
            cpu_micros: timeval_micros(value.ru_utime)
                .saturating_add(timeval_micros(value.ru_stime)),
            minor_faults: u64::try_from(value.ru_minflt).unwrap_or(0),
            major_faults: u64::try_from(value.ru_majflt).unwrap_or(0),
        }
    }

    fn process() -> Self {
        Self::for_who(libc::RUSAGE_SELF)
    }

    fn thread() -> Self {
        #[cfg(target_os = "linux")]
        {
            Self::for_who(libc::RUSAGE_THREAD)
        }
        #[cfg(not(target_os = "linux"))]
        {
            Self::process()
        }
    }

    fn subtract(self, before: Self) -> Self {
        Self {
            cpu_micros: self.cpu_micros.saturating_sub(before.cpu_micros),
            minor_faults: self.minor_faults.saturating_sub(before.minor_faults),
            major_faults: self.major_faults.saturating_sub(before.major_faults),
        }
    }
}

fn timeval_micros(value: libc::timeval) -> u64 {
    u64::try_from(value.tv_sec)
        .unwrap_or(0)
        .saturating_mul(1_000_000)
        .saturating_add(u64::try_from(value.tv_usec).unwrap_or(0))
}

fn process_read_bytes() -> u64 {
    std::fs::read_to_string("/proc/self/io")
        .ok()
        .and_then(|text| {
            text.lines().find_map(|line| {
                line.strip_prefix("read_bytes:")
                    .and_then(|value| value.trim().parse::<u64>().ok())
            })
        })
        .unwrap_or(0)
}

#[derive(Debug)]
struct StageState {
    calls: u64,
    process_wall: Duration,
    cpu_micros: u64,
    minor_faults: u64,
    major_faults: u64,
    physical_read_bytes: u64,
    active_workers: u64,
    maximum_active_workers: u64,
    last_activity_change: Instant,
    busy_wall: Duration,
    worker_wall: Duration,
}

impl StageState {
    fn new(now: Instant) -> Self {
        Self {
            calls: 0,
            process_wall: Duration::ZERO,
            cpu_micros: 0,
            minor_faults: 0,
            major_faults: 0,
            physical_read_bytes: 0,
            active_workers: 0,
            maximum_active_workers: 0,
            last_activity_change: now,
            busy_wall: Duration::ZERO,
            worker_wall: Duration::ZERO,
        }
    }

    fn update_activity(&mut self, now: Instant) {
        let elapsed = now.saturating_duration_since(self.last_activity_change);
        if self.active_workers != 0 {
            self.busy_wall = self.busy_wall.saturating_add(elapsed);
            self.worker_wall = self.worker_wall.saturating_add(
                elapsed
                    .checked_mul(u32::try_from(self.active_workers).unwrap_or(u32::MAX))
                    .unwrap_or(Duration::MAX),
            );
        }
        self.last_activity_change = now;
    }

    fn enter_worker(&mut self, now: Instant) {
        self.update_activity(now);
        self.active_workers = self.active_workers.saturating_add(1);
        self.maximum_active_workers = self.maximum_active_workers.max(self.active_workers);
    }

    fn leave_worker(&mut self, now: Instant) {
        self.update_activity(now);
        self.active_workers = self.active_workers.saturating_sub(1);
    }
}

#[derive(Debug)]
struct Profiler {
    started: Instant,
    start_usage: Usage,
    start_read_bytes: u64,
    stages: Mutex<BTreeMap<String, StageState>>,
    counters: Mutex<BTreeMap<String, u64>>,
}

impl Profiler {
    fn new() -> Self {
        Self {
            started: Instant::now(),
            start_usage: Usage::process(),
            start_read_bytes: process_read_bytes(),
            stages: Mutex::new(BTreeMap::new()),
            counters: Mutex::new(BTreeMap::new()),
        }
    }

    fn enter_worker(&self, name: &str) {
        let now = Instant::now();
        let mut stages = self
            .stages
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        stages
            .entry(name.to_string())
            .or_insert_with(|| StageState::new(now))
            .enter_worker(now);
    }

    fn leave_worker(&self, name: &str) {
        let now = Instant::now();
        let mut stages = self
            .stages
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if let Some(stage) = stages.get_mut(name) {
            stage.leave_worker(now);
        }
    }

    fn record_thread(&self, name: &str, usage: Usage) {
        let now = Instant::now();
        let mut stages = self
            .stages
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        let stage = stages
            .entry(name.to_string())
            .or_insert_with(|| StageState::new(now));
        stage.calls = stage.calls.saturating_add(1);
        stage.cpu_micros = stage.cpu_micros.saturating_add(usage.cpu_micros);
        stage.minor_faults = stage.minor_faults.saturating_add(usage.minor_faults);
        stage.major_faults = stage.major_faults.saturating_add(usage.major_faults);
    }

    fn record_process(&self, name: &str, wall: Duration, usage: Usage, physical_read_bytes: u64) {
        let now = Instant::now();
        let mut stages = self
            .stages
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        let stage = stages
            .entry(name.to_string())
            .or_insert_with(|| StageState::new(now));
        stage.calls = stage.calls.saturating_add(1);
        stage.process_wall = stage.process_wall.saturating_add(wall);
        stage.cpu_micros = stage.cpu_micros.saturating_add(usage.cpu_micros);
        stage.minor_faults = stage.minor_faults.saturating_add(usage.minor_faults);
        stage.major_faults = stage.major_faults.saturating_add(usage.major_faults);
        stage.physical_read_bytes = stage
            .physical_read_bytes
            .saturating_add(physical_read_bytes);
    }

    fn report(&self) -> StageProfileReport {
        let now = Instant::now();
        let stages = self
            .stages
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        let stages = stages
            .iter()
            .map(|(name, state)| {
                let mut busy_wall = state.busy_wall;
                let mut worker_wall = state.worker_wall;
                if state.active_workers != 0 {
                    let elapsed = now.saturating_duration_since(state.last_activity_change);
                    busy_wall = busy_wall.saturating_add(elapsed);
                    worker_wall = worker_wall.saturating_add(
                        elapsed
                            .checked_mul(u32::try_from(state.active_workers).unwrap_or(u32::MAX))
                            .unwrap_or(Duration::MAX),
                    );
                }
                let reported_wall = if state.process_wall.is_zero() {
                    busy_wall
                } else {
                    state.process_wall
                };
                let average_active_workers = if busy_wall.is_zero() {
                    0.0
                } else {
                    worker_wall.as_secs_f64() / busy_wall.as_secs_f64()
                };
                StageProfile {
                    name: name.clone(),
                    calls: state.calls,
                    wall_micros: micros(reported_wall),
                    cpu_micros: state.cpu_micros,
                    minor_page_faults: state.minor_faults,
                    major_page_faults: state.major_faults,
                    physical_read_bytes: state.physical_read_bytes,
                    average_active_workers,
                    maximum_active_workers: state.maximum_active_workers,
                }
            })
            .collect();
        let usage = Usage::process().subtract(self.start_usage);
        StageProfileReport {
            wall_micros: micros(self.started.elapsed()),
            cpu_micros: usage.cpu_micros,
            minor_page_faults: usage.minor_faults,
            major_page_faults: usage.major_faults,
            physical_read_bytes: process_read_bytes().saturating_sub(self.start_read_bytes),
            stages,
            counters: self
                .counters
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .clone(),
        }
    }
}

fn micros(value: Duration) -> u64 {
    value.as_micros().try_into().unwrap_or(u64::MAX)
}

#[derive(Clone, Debug, Serialize)]
pub struct StageProfile {
    pub name: String,
    pub calls: u64,
    pub wall_micros: u64,
    pub cpu_micros: u64,
    pub minor_page_faults: u64,
    pub major_page_faults: u64,
    pub physical_read_bytes: u64,
    pub average_active_workers: f64,
    pub maximum_active_workers: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct StageProfileReport {
    pub wall_micros: u64,
    pub cpu_micros: u64,
    pub minor_page_faults: u64,
    pub major_page_faults: u64,
    pub physical_read_bytes: u64,
    pub stages: Vec<StageProfile>,
    pub counters: BTreeMap<String, u64>,
}

pub struct ProfileSession {
    profiler: Arc<Profiler>,
}

impl ProfileSession {
    #[must_use]
    pub fn start() -> Self {
        let profiler = Arc::new(Profiler::new());
        *active_slot()
            .lock()
            .unwrap_or_else(|error| error.into_inner()) = Some(Arc::clone(&profiler));
        Self { profiler }
    }

    #[must_use]
    pub fn report(&self) -> StageProfileReport {
        self.profiler.report()
    }
}

impl Drop for ProfileSession {
    fn drop(&mut self) {
        let mut active = active_slot()
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if active
            .as_ref()
            .is_some_and(|profiler| Arc::ptr_eq(profiler, &self.profiler))
        {
            *active = None;
        }
    }
}

fn current() -> Option<Arc<Profiler>> {
    active_slot()
        .lock()
        .unwrap_or_else(|error| error.into_inner())
        .clone()
}

pub struct StageGuard {
    profiler: Option<Arc<Profiler>>,
    name: &'static str,
    usage: Usage,
}

impl Drop for StageGuard {
    fn drop(&mut self) {
        if let Some(profiler) = &self.profiler {
            profiler.record_thread(self.name, Usage::thread().subtract(self.usage));
            profiler.leave_worker(self.name);
        }
    }
}

#[must_use]
pub fn scope(name: &'static str) -> StageGuard {
    let profiler = current();
    if let Some(profiler) = &profiler {
        profiler.enter_worker(name);
    }
    StageGuard {
        profiler,
        name,
        usage: Usage::thread(),
    }
}

pub struct ProcessStageGuard {
    profiler: Option<Arc<Profiler>>,
    name: &'static str,
    started: Instant,
    usage: Usage,
    read_bytes: u64,
}

impl Drop for ProcessStageGuard {
    fn drop(&mut self) {
        if let Some(profiler) = &self.profiler {
            profiler.record_process(
                self.name,
                self.started.elapsed(),
                Usage::process().subtract(self.usage),
                process_read_bytes().saturating_sub(self.read_bytes),
            );
        }
    }
}

#[must_use]
pub fn process_scope(name: &'static str) -> ProcessStageGuard {
    ProcessStageGuard {
        profiler: current(),
        name,
        started: Instant::now(),
        usage: Usage::process(),
        read_bytes: process_read_bytes(),
    }
}

pub struct WorkerGuard {
    profiler: Option<Arc<Profiler>>,
    name: &'static str,
}

impl Drop for WorkerGuard {
    fn drop(&mut self) {
        if let Some(profiler) = &self.profiler {
            profiler.leave_worker(self.name);
        }
    }
}

#[must_use]
pub fn worker(name: &'static str) -> WorkerGuard {
    let profiler = current();
    if let Some(profiler) = &profiler {
        profiler.enter_worker(name);
    }
    WorkerGuard { profiler, name }
}

pub fn add_counter(name: &'static str, value: u64) {
    let Some(profiler) = current() else {
        return;
    };
    let mut counters = profiler
        .counters
        .lock()
        .unwrap_or_else(|error| error.into_inner());
    let counter = counters.entry(name.to_string()).or_default();
    *counter = counter.saturating_add(value);
}
