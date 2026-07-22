//! Bounded, ordered parsing for row-oriented supplementary sources.
//!
//! A producer splits decoded input on complete line boundaries. Parser workers
//! process independent batches, while the iterator restores input order before
//! records reach the OSA writer. Both queues are bounded so whole-genome builds
//! use a small, predictable amount of memory.

use anyhow::{anyhow, Context, Result};
use fastvep_sa::common::AnnotationRecord;
use std::collections::{BTreeMap, VecDeque};
use std::io::BufRead;
use std::sync::{mpsc, Arc, Mutex};

const DEFAULT_BATCH_BYTES: usize = 4 * 1024 * 1024;

pub(crate) type BatchParser = Arc<dyn Fn(Vec<u8>) -> Vec<Result<AnnotationRecord>> + Send + Sync>;

struct BatchJob {
    sequence: u64,
    bytes: Vec<u8>,
}

struct ParsedBatch {
    sequence: u64,
    records: Vec<Result<AnnotationRecord>>,
}

/// Iterator over records parsed concurrently but yielded in original order.
pub(crate) struct OrderedParallelRecordIter {
    receiver: mpsc::Receiver<ParsedBatch>,
    pending: BTreeMap<u64, VecDeque<Result<AnnotationRecord>>>,
    current: VecDeque<Result<AnnotationRecord>>,
    next_sequence: u64,
    closed: bool,
}

impl OrderedParallelRecordIter {
    pub(crate) fn new<R>(reader: R, parser: BatchParser, worker_count: usize) -> Self
    where
        R: BufRead + Send + 'static,
    {
        Self::with_batch_size(reader, parser, worker_count, DEFAULT_BATCH_BYTES)
    }

    fn with_batch_size<R>(
        reader: R,
        parser: BatchParser,
        worker_count: usize,
        batch_bytes: usize,
    ) -> Self
    where
        R: BufRead + Send + 'static,
    {
        let worker_count = worker_count.clamp(1, 4);
        let queue_capacity = worker_count * 2;
        let (job_sender, job_receiver) = mpsc::sync_channel::<BatchJob>(queue_capacity);
        let job_receiver = Arc::new(Mutex::new(job_receiver));
        let (result_sender, result_receiver) = mpsc::sync_channel::<ParsedBatch>(queue_capacity);

        let producer_results = result_sender.clone();
        std::thread::spawn(move || {
            let mut reader = reader;
            let mut sequence = 0;
            loop {
                match read_batch(&mut reader, batch_bytes) {
                    Ok(bytes) if bytes.is_empty() => return,
                    Ok(bytes) => {
                        if job_sender.send(BatchJob { sequence, bytes }).is_err() {
                            return;
                        }
                        sequence += 1;
                    }
                    Err(error) => {
                        let _ = producer_results.send(ParsedBatch {
                            sequence,
                            records: vec![Err(error).context("reading supplementary input batch")],
                        });
                        return;
                    }
                }
            }
        });

        for _ in 0..worker_count {
            let jobs = Arc::clone(&job_receiver);
            let results = result_sender.clone();
            let parser = Arc::clone(&parser);
            std::thread::spawn(move || loop {
                let job = match jobs.lock() {
                    Ok(receiver) => receiver.recv(),
                    Err(_) => return,
                };
                let Ok(job) = job else {
                    return;
                };
                let records = parser(job.bytes)
                    .into_iter()
                    .map(|record| {
                        record.with_context(|| {
                            format!("parsing supplementary input batch {}", job.sequence + 1)
                        })
                    })
                    .collect();
                if results
                    .send(ParsedBatch {
                        sequence: job.sequence,
                        records,
                    })
                    .is_err()
                {
                    return;
                }
            });
        }
        drop(result_sender);

        Self {
            receiver: result_receiver,
            pending: BTreeMap::new(),
            current: VecDeque::new(),
            next_sequence: 0,
            closed: false,
        }
    }
}

impl Iterator for OrderedParallelRecordIter {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.current.pop_front() {
                return Some(record);
            }
            if let Some(records) = self.pending.remove(&self.next_sequence) {
                self.current = records;
                self.next_sequence += 1;
                continue;
            }
            if self.closed {
                return if self.pending.is_empty() {
                    None
                } else {
                    self.pending.clear();
                    Some(Err(anyhow!("parallel parser result sequence has a gap")))
                };
            }
            match self.receiver.recv() {
                Ok(batch) => {
                    self.pending.insert(batch.sequence, batch.records.into());
                }
                Err(_) => self.closed = true,
            }
        }
    }
}

fn read_batch(reader: &mut impl BufRead, target_bytes: usize) -> std::io::Result<Vec<u8>> {
    let mut batch = Vec::with_capacity(target_bytes + 64 * 1024);
    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Ok(batch);
        }
        let remaining = target_bytes.saturating_sub(batch.len());
        if remaining < available.len() {
            let search_start = remaining.saturating_sub(1);
            if let Some(relative_newline) = available[search_start..]
                .iter()
                .position(|byte| *byte == b'\n')
            {
                let take = search_start + relative_newline + 1;
                batch.extend_from_slice(&available[..take]);
                reader.consume(take);
                return Ok(batch);
            }
        }
        let take = available.len();
        batch.extend_from_slice(available);
        reader.consume(take);
        if batch.len() >= target_bytes && batch.last() == Some(&b'\n') {
            return Ok(batch);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::time::Duration;

    fn record(position: u32) -> AnnotationRecord {
        AnnotationRecord {
            chrom_idx: 0,
            position,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        }
    }

    #[test]
    fn restores_order_when_later_batch_finishes_first() {
        let parser: BatchParser = Arc::new(|bytes| {
            let value = String::from_utf8(bytes).unwrap();
            if value.starts_with("slow") {
                std::thread::sleep(Duration::from_millis(25));
                vec![Ok(record(1))]
            } else {
                vec![Ok(record(2))]
            }
        });
        let records = OrderedParallelRecordIter::with_batch_size(
            Cursor::new(b"slow\nfast\n".to_vec()),
            parser,
            2,
            5,
        )
        .collect::<Result<Vec<_>>>()
        .unwrap();
        assert_eq!(
            records
                .iter()
                .map(|record| record.position)
                .collect::<Vec<_>>(),
            vec![1, 2]
        );
    }

    #[test]
    fn parser_errors_are_returned_in_input_order() {
        let parser: BatchParser = Arc::new(|bytes| {
            if bytes.starts_with(b"bad") {
                vec![Err(anyhow!("bad row"))]
            } else {
                vec![Ok(record(1))]
            }
        });
        let mut records = OrderedParallelRecordIter::with_batch_size(
            Cursor::new(b"ok__\nbad_\n".to_vec()),
            parser,
            2,
            5,
        );
        assert_eq!(records.next().unwrap().unwrap().position, 1);
        assert!(records
            .next()
            .unwrap()
            .unwrap_err()
            .to_string()
            .contains("batch 2"));
    }
}
