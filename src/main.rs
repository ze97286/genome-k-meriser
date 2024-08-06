use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::sync::{Arc, Mutex};
use rust_htslib::faidx::Reader as Faidx;
use log::{info, error};
use std::time::Instant;

const CHUNK_SIZE: usize = 10_000_000; // 10M bases

/// Reads the reference genome index file (.fai) and extracts the name, length, and offset of each chromosome.
///
/// This function parses a `.fai` index file which contains information about chromosome names, lengths, and offsets.
/// It returns a `HashMap` where the keys are chromosome names and the values are tuples containing
/// the offset and length of each chromosome.
///
/// # Arguments
/// * `fai_path` - Path to the `.fai` file.
///
/// # Returns
/// Returns a `HashMap<String, (usize, usize)>` where:
/// - The key is the chromosome name.
/// - The value is a tuple where the first element is the offset of the chromosome in the FASTA file, 
///   and the second element is the length of the chromosome.
///
/// # Errors
/// Returns an `io::Result<HashMap<String, (usize, usize)>>` which will be an error if the file cannot be read 
/// or if there is a parsing error.
fn read_fai_index(fai_path: &str) -> io::Result<HashMap<String, (usize, usize)>> {
    let mut index = HashMap::new();
    let file = File::open(fai_path)?;
    for line in io::BufReader::new(file).lines() {
        let line = line?;
        let mut parts = line.split('\t');
        let name = parts.next().unwrap().to_string();
        let length = parts.next().unwrap().parse::<usize>().unwrap();
        let offset = parts.next().unwrap().parse::<usize>().unwrap();
        index.insert(name, (offset, length));
    }
    Ok(index)
}

/// Finds k-mers in a given sequence.
///
/// This function scans through a given sequence to find all k-mers of a specified length. It applies an offset 
/// to the positions of the k-mers and returns a `HashMap` where keys are k-mer sequences and values are tuples
/// containing the count of occurrences and a list of positions where the k-mer is found.
///
/// # Arguments
/// * `seq` - The sequence data as a byte slice.
/// * `offset` - The offset to be applied to k-mer positions.
/// * `k` - The length of the k-mers to find.
///
/// # Returns
/// Returns a `HashMap<String, (usize, Vec<usize>)>` where:
/// - The key is the k-mer sequence.
/// - The value is a tuple where the first element is the count of occurrences of the k-mer,
///   and the second element is a vector of positions (with the offset applied) where the k-mer is found.
///
/// # Notes
/// Skips k-mers consisting only of 'N' characters.
fn find_kmers(seq: &[u8], offset: usize, k: usize) -> HashMap<String, (usize, Vec<usize>)> {
    let mut kmers = HashMap::new();
    let binding = seq.to_ascii_uppercase();
    let seq_str = std::str::from_utf8(&binding).expect("Error converting sequence to string");

    for i in 0..=(seq.len() - k) {
        let kmer = &seq_str[i..i + k];
        if kmer.chars().all(|c| c == 'N') {
            continue; // Skip k-mers consisting only of 'N'
        }
        let entry = kmers.entry(kmer.to_string()).or_insert_with(|| (0, Vec::new()));
        entry.0 += 1; // Update count
        entry.1.push(i + offset); // Apply offset here
    }
    kmers
}

/// Processes a single chromosome to find k-mers and saves the results to a TSV file.
///
/// # Arguments
/// * `fasta_path` - Path to the FASTA file containing the genomic sequences.
/// * `fai_index` - An `Arc<Mutex<HashMap<String, (usize, usize)>>>` providing access to the index of the reference genome.
/// * `contig_name` - The name of the chromosome to process.
/// * `k` - The length of k-mers to find.
/// * `output_dir` - The directory where TSV files will be saved.
///
/// # Returns
/// Returns an `io::Result` containing a `HashMap` where keys are k-mer sequences and values are tuples containing
/// the count of occurrences and a list of positions where the k-mer is found.
///
/// # Errors
/// Returns an `io::Result` that will be an error if sequence fetching fails or if results cannot be saved to the file.
fn process_chromosome(
    fasta_path: &str,
    fai_index: Arc<Mutex<HashMap<String, (usize, usize)>>>,
    contig_name: &str,
    k: usize,
    output_dir: &str,
) -> io::Result<()> {
    info!("Starting process_chromosome for {}", contig_name);
    let faidx = Faidx::from_path(fasta_path).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    let mut results = HashMap::new();

    let (start_offset, length) = {
        let binding = fai_index.lock().unwrap();
        binding.get(contig_name).expect("Contig not found in FAI index").clone()
    };

    let mut start = 0; // Start fetching from position 0 within the chunk
    let mut end = CHUNK_SIZE;

    loop {
        // Log the current fetching progress
        let progress = if length > 0 {
            let completed = if start + CHUNK_SIZE > length {
                length
            } else {
                start + CHUNK_SIZE
            };
            let percentage = (completed as f64 / length as f64) * 100.0;
            format!("{:.2}%", percentage)
        } else {
            "N/A".to_string()
        };

        info!("Fetching sequence {} from {} to {} (out of {}), Progress: {}", contig_name, start, end, length, progress);

        // Fetch sequence data
        let seq_data = match faidx.fetch_seq(contig_name, start as usize, end as usize) {
            Ok(data) => data,
            Err(e) => {
                error!("Error fetching sequence {}: {:?}", contig_name, e);
                break;
            }
        };

        if seq_data.is_empty() {
            break;
        }

        // Skip sequences consisting only of 'N'
        if seq_data.iter().all(|&c| c == b'N') {
            info!("Skipping sequence {}: consists only of 'N'", contig_name);
        } else {
            // Find k-mers in the fetched sequence
            let kmers = find_kmers(&seq_data, start, k);
            for (kmer, (count, positions)) in kmers {
                // Accumulate counts and positions for each k-mer
                let entry = results.entry(kmer).or_insert_with(|| (0, Vec::new()));
                entry.0 += count; // Correctly accumulate count
                entry.1.extend(positions);
            }
        }

        start = end - k + 1; // Overlap handling
        end = start + CHUNK_SIZE;

        if seq_data.len() < CHUNK_SIZE {
            break;
        }
    }
    // Save results to TSV file
    if let Err(e) = save_results_to_tsv(output_dir, contig_name, results.clone(), start_offset) {
        eprintln!("Error saving results for {}: {:?}", contig_name, e);
    }
    info!("Finished process_chromosome for {}", contig_name);
    Ok(())
}

/// Saves k-mer results to a TSV file for a given contig.
///
/// # Arguments
/// * `output_dir` - The directory where the TSV file will be saved.
/// * `contig_name` - The name of the contig for which results are saved.
/// * `kmers` - A `HashMap` where keys are k-mer sequences and values are tuples containing
///   the count of occurrences and a list of positions (with offset applied) where the k-mer is found.
/// * `start_offset` - The offset to be applied to positions when saving results.
///
/// # Returns
/// Returns an `io::Result` which will be an error if the file cannot be created or written to.
fn save_results_to_tsv(
    output_dir: &str,
    contig_name: &str,
    kmers: HashMap<String, (usize, Vec<usize>)>,
    start_offset: usize,
) -> io::Result<()> {
    let file_path = format!("{}/{}.tsv", output_dir, contig_name);
    let mut file = File::create(file_path)?;

    writeln!(file, "contig\tkmer\tpositions_chromosome\tpositions_genome\tnum_reps")?;
    for (kmer, (count, positions)) in kmers {
        let positions_chromosome = positions.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(",");
        let positions_genome = positions.iter().map(|p| (p + start_offset).to_string()).collect::<Vec<_>>().join(",");
        writeln!(file, "{}\t{}\t{}\t{}\t{}", contig_name, kmer, positions_chromosome, positions_genome, count)?;
    }
    Ok(())
}

fn main() {
    env_logger::init();

    let fasta_path = env::args().nth(1).expect("Please provide the path to the FASTA file");
    let fai_path = env::args().nth(2).expect("Please provide the path to the FAI file");
    let k: usize = env::args().nth(3).expect("Please provide the k-mer length").parse().expect("Invalid k-mer length");
    let num_threads: usize = env::args().nth(4).expect("Please provide the number of threads").parse().expect("Invalid number of threads");
    let output_dir = env::args().nth(5).expect("Please provide the output directory");

    let start_time = Instant::now();

    // Read the FAI index and wrap it in an Arc<Mutex> for safe sharing between threads
    let fai_index = Arc::new(Mutex::new(read_fai_index(&fai_path).expect("Error reading FAI index")));

    info!("Starting processing for all contigs");

    let fai_index = Arc::clone(&fai_index);
    let contig_names: Vec<String> = fai_index.lock().unwrap().keys().cloned().collect();
    
    // Set up the thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Failed to build thread pool");

    // Process each chromosome in parallel
    contig_names.into_par_iter().for_each(|contig_name| {
        let fasta_path = fasta_path.clone();
        let fai_index = Arc::clone(&fai_index);
        let k = k;
        let output_dir = output_dir.clone();
    
        if let Err(e) = process_chromosome(&fasta_path, fai_index, &contig_name, k, &output_dir) {
            eprintln!("Error processing chromosome {}: {:?}", contig_name, e);
        }
    });

    info!("Processing completed");

    let elapsed_time = start_time.elapsed();
    info!("Total processing time: {:.2?}", elapsed_time);
}
