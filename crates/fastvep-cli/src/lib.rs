//! Library entry point for the fastVEP CLI.
//!
//! Exposes the pipeline + webserver modules so they can be exercised by
//! integration tests in `tests/`. The binary entry point is `src/main.rs`.

mod parallel_records;
pub mod pipeline;
pub mod progress;
pub mod webserver;
