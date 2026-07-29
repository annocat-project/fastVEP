//! Writer for .osa position/allele-level annotation files.
//!
//! Records must be added in chromosome-sorted, position-sorted order.
//! The writer accumulates entries into blocks, compresses them, and writes
//! to the data file while building the index.

use crate::block::{BlockEntry, SaBlock};
use crate::common::{AnnotationRecord, DEFAULT_BLOCK_SIZE, OSA_MAGIC};
use crate::index::{BlockRef, IndexHeader, SaIndex};
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{mpsc, Arc, Mutex};

const FILE_BUFFER_SIZE: usize = 1024 * 1024;

struct BlockJob {
    sequence: u64,
    chrom: String,
    start_pos: u32,
    end_pos: u32,
    block: SaBlock,
}

struct CompressedBlock {
    sequence: u64,
    chrom: String,
    start_pos: u32,
    end_pos: u32,
    compressed: Result<Vec<u8>>,
}

/// Builds an .osa data file and its .osa.idx index file.
pub struct SaWriter {
    index: SaIndex,
    block: SaBlock,
    current_chrom: Option<String>,
    last_key: Option<(u16, u32)>,
    /// Chromosome name -> numeric index mapping.
    chrom_names: Vec<String>,
    data_offset: u64,
}

impl SaWriter {
    pub fn new(header: IndexHeader) -> Self {
        Self {
            index: SaIndex::new(header),
            block: SaBlock::new(DEFAULT_BLOCK_SIZE),
            current_chrom: None,
            last_key: None,
            chrom_names: Vec::new(),
            data_offset: 0,
        }
    }

    /// Build .osa and .osa.idx from an iterator of sorted annotation records.
    ///
    /// Records MUST be sorted by (chrom_idx, position).
    /// `chrom_map` maps chrom_idx -> chromosome name string.
    pub fn write_all<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = AnnotationRecord>,
        chrom_map: &[String],
    ) -> Result<()> {
        self.write_all_results(data_writer, records.map(Ok), chrom_map)
    }

    /// Build .osa and .osa.idx from an iterator that can surface parse errors.
    ///
    /// Records MUST be sorted by (chrom_idx, position).
    pub fn write_all_results<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<()> {
        self.write_all_results_with_workers(
            data_writer,
            records,
            chrom_map,
            compression_worker_count(),
        )
    }

    fn write_all_results_with_workers<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
        worker_count: usize,
    ) -> Result<()> {
        self.chrom_names = chrom_map.to_vec();

        data_writer.write_all(OSA_MAGIC)?;
        data_writer.write_all(&self.index.header.schema_version.to_le_bytes())?;
        self.data_offset = (OSA_MAGIC.len() + 2) as u64;

        let worker_count = worker_count.max(1);
        let queue_capacity = worker_count * 2;
        let (job_sender, job_receiver) = mpsc::sync_channel::<BlockJob>(queue_capacity);
        let job_receiver = Arc::new(Mutex::new(job_receiver));
        let (result_sender, result_receiver) = mpsc::channel::<CompressedBlock>();

        std::thread::scope(|scope| -> Result<()> {
            for _ in 0..worker_count {
                let jobs = Arc::clone(&job_receiver);
                let results = result_sender.clone();
                scope.spawn(move || loop {
                    let job = match jobs.lock() {
                        Ok(receiver) => receiver.recv(),
                        Err(_) => return,
                    };
                    let Ok(job) = job else {
                        return;
                    };
                    let compressed = job.block.compress();
                    if results
                        .send(CompressedBlock {
                            sequence: job.sequence,
                            chrom: job.chrom,
                            start_pos: job.start_pos,
                            end_pos: job.end_pos,
                            compressed,
                        })
                        .is_err()
                    {
                        return;
                    }
                });
            }
            drop(result_sender);

            let mut job_sender = Some(job_sender);
            let mut submitted = 0u64;
            let mut next_to_write = 0u64;
            let mut pending = BTreeMap::new();

            let pipeline_result = (|| -> Result<()> {
                for record in records {
                    let record = record?;
                    if let Some((last_chrom, last_pos)) = self.last_key {
                        if (record.chrom_idx, record.position) < (last_chrom, last_pos) {
                            anyhow::bail!(
                                "SA records are not sorted: previous chrom_idx={}, position={}; current chrom_idx={}, position={}. \
                                 The streaming .osa builder requires input sorted by chromosome (chr1..chr22,X,Y,M) then position \
                                 — sort the source file (e.g. `bcftools sort` / `sort -k1,1 -k2,2n`) and rebuild.",
                                last_chrom,
                                last_pos,
                                record.chrom_idx,
                                record.position
                            );
                        }
                    }
                    self.last_key = Some((record.chrom_idx, record.position));

                    let chrom_name = &chrom_map[record.chrom_idx as usize];
                    if self.current_chrom.as_ref() != Some(chrom_name) {
                        self.submit_block(job_sender.as_ref().unwrap(), submitted)?;
                        if !self.block.is_empty() {
                            unreachable!("submitted block was not cleared");
                        }
                        if self.current_chrom.is_some() {
                            submitted += 1;
                        }
                        self.current_chrom = Some(chrom_name.clone());
                        self.drain_compressed(
                            data_writer,
                            &result_receiver,
                            &mut pending,
                            &mut next_to_write,
                        )?;
                    }

                    let entry = BlockEntry {
                        position: record.position,
                        ref_allele: record.ref_allele,
                        alt_allele: record.alt_allele,
                        json: record.json,
                    };
                    if !self.block.can_add(&entry) {
                        self.submit_block(job_sender.as_ref().unwrap(), submitted)?;
                        submitted += 1;
                        self.drain_compressed(
                            data_writer,
                            &result_receiver,
                            &mut pending,
                            &mut next_to_write,
                        )?;
                    }
                    assert!(self.block.add(entry), "Single entry exceeds block size");
                }

                if !self.block.is_empty() {
                    self.submit_block(job_sender.as_ref().unwrap(), submitted)?;
                    submitted += 1;
                }
                drop(job_sender.take());

                while next_to_write < submitted {
                    let completed = result_receiver
                        .recv()
                        .context("OSA compression worker stopped before completing all blocks")?;
                    pending.insert(completed.sequence, completed);
                    self.write_ready_blocks(data_writer, &mut pending, &mut next_to_write)?;
                }
                Ok(())
            })();

            drop(job_sender.take());
            pipeline_result
        })
    }

    fn submit_block(&mut self, sender: &mpsc::SyncSender<BlockJob>, sequence: u64) -> Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        let block = std::mem::replace(&mut self.block, SaBlock::new(DEFAULT_BLOCK_SIZE));
        sender
            .send(BlockJob {
                sequence,
                chrom: self.current_chrom.as_ref().unwrap().clone(),
                start_pos: block.start_position().unwrap(),
                end_pos: block.end_position().unwrap(),
                block,
            })
            .context("OSA compression workers stopped unexpectedly")
    }

    fn drain_compressed<W: Write>(
        &mut self,
        writer: &mut W,
        receiver: &mpsc::Receiver<CompressedBlock>,
        pending: &mut BTreeMap<u64, CompressedBlock>,
        next_to_write: &mut u64,
    ) -> Result<()> {
        while let Ok(completed) = receiver.try_recv() {
            pending.insert(completed.sequence, completed);
        }
        self.write_ready_blocks(writer, pending, next_to_write)
    }

    fn write_ready_blocks<W: Write>(
        &mut self,
        writer: &mut W,
        pending: &mut BTreeMap<u64, CompressedBlock>,
        next_to_write: &mut u64,
    ) -> Result<()> {
        while let Some(completed) = pending.remove(next_to_write) {
            let compressed = completed.compressed?;
            let compressed_len = compressed.len() as u32;
            writer.write_all(&compressed_len.to_le_bytes())?;
            writer.write_all(&compressed)?;
            self.index.add_block(
                &completed.chrom,
                BlockRef {
                    start_pos: completed.start_pos,
                    end_pos: completed.end_pos,
                    file_offset: self.data_offset,
                    compressed_len,
                },
            );
            self.data_offset += 4 + compressed_len as u64;
            *next_to_write += 1;
        }
        Ok(())
    }

    /// Write the index file.
    pub fn write_index<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.index.write_to(writer)
    }

    /// Convenience: write .osa and .osa.idx to files at the given base path.
    pub fn write_to_files(
        &mut self,
        base_path: &Path,
        records: impl Iterator<Item = AnnotationRecord>,
        chrom_map: &[String],
    ) -> Result<()> {
        let data_path = base_path.with_extension("osa");
        let idx_path = base_path.with_extension("osa.idx");

        let data_file = std::fs::File::create(&data_path).with_context(|| {
            format!(
                "Creating output file {} (does the output directory exist?)",
                data_path.display()
            )
        })?;
        let mut data_writer = BufWriter::with_capacity(FILE_BUFFER_SIZE, data_file);
        self.write_all(&mut data_writer, records, chrom_map)?;
        data_writer.flush()?;

        let idx_file = std::fs::File::create(&idx_path)
            .with_context(|| format!("Creating index file {}", idx_path.display()))?;
        let mut idx_writer = BufWriter::with_capacity(FILE_BUFFER_SIZE, idx_file);
        self.write_index(&mut idx_writer)?;
        idx_writer.flush()?;

        Ok(())
    }

    /// Convenience: write .osa and .osa.idx to files from fallible records.
    pub fn write_results_to_files(
        &mut self,
        base_path: &Path,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<()> {
        let data_path = base_path.with_extension("osa");
        let idx_path = base_path.with_extension("osa.idx");

        let data_file = std::fs::File::create(&data_path).with_context(|| {
            format!(
                "Creating output file {} (does the output directory exist?)",
                data_path.display()
            )
        })?;
        let mut data_writer = BufWriter::with_capacity(FILE_BUFFER_SIZE, data_file);
        self.write_all_results(&mut data_writer, records, chrom_map)?;
        data_writer.flush()?;

        let idx_file = std::fs::File::create(&idx_path)
            .with_context(|| format!("Creating index file {}", idx_path.display()))?;
        let mut idx_writer = BufWriter::with_capacity(FILE_BUFFER_SIZE, idx_file);
        self.write_index(&mut idx_writer)?;
        idx_writer.flush()?;

        Ok(())
    }
}

fn compression_worker_count() -> usize {
    if let Ok(value) = std::env::var("FASTVEP_SA_COMPRESSION_THREADS") {
        if let Ok(count) = value.parse::<usize>() {
            if (1..=8).contains(&count) {
                return count;
            }
        }
    }
    // One compression worker already overlaps zstd with parsing on the caller
    // thread. Additional workers are opt-in because AnnoCat may build several
    // sources concurrently and should own the global CPU budget.
    1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> IndexHeader {
        IndexHeader {
            schema_version: crate::common::SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test".into(),
            version: "test".into(),
            description: "test".into(),
            assembly: "GRCh38".into(),
            match_by_allele: true,
            is_array: false,
            record_list: false,
            is_positional: false,
        }
    }

    fn record(chrom_idx: u16, position: u32) -> AnnotationRecord {
        AnnotationRecord {
            chrom_idx,
            position,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        }
    }

    #[test]
    fn write_all_rejects_unsorted_records() {
        let mut writer = SaWriter::new(header());
        let mut out = Vec::new();
        let err = writer
            .write_all(
                &mut out,
                vec![record(0, 20), record(0, 10)].into_iter(),
                &["1".into()],
            )
            .unwrap_err();

        assert!(err.to_string().contains("SA records are not sorted"));
    }

    #[test]
    fn parallel_compression_is_byte_identical_and_ordered() {
        fn build(worker_count: usize) -> (Vec<u8>, Vec<u8>) {
            let records = (0..10_000).map(|position| {
                Ok(AnnotationRecord {
                    chrom_idx: 0,
                    position,
                    ref_allele: "A".into(),
                    alt_allele: "G".into(),
                    json: format!(
                        r#"{{"score":{},"padding":"{}"}}"#,
                        position,
                        "x".repeat(1024)
                    ),
                })
            });
            let mut writer = SaWriter::new(header());
            let mut data = Vec::new();
            writer
                .write_all_results_with_workers(&mut data, records, &["1".into()], worker_count)
                .unwrap();
            let mut index = Vec::new();
            writer.write_index(&mut index).unwrap();
            (data, index)
        }

        let serial_worker = build(1);
        let parallel_workers = build(3);
        assert_eq!(serial_worker, parallel_workers);
    }
}
