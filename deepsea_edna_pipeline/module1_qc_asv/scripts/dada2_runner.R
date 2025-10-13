# DADA2 denoising script for ASV generation
# This R script is called from Python to perform DADA2 denoising

# Load required libraries
suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(methods)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript dada2_runner.R <input_dir> <output_dir> <sample_names_file> [additional_params]")
}

input_dir <- args[1]
output_dir <- args[2] 
sample_names_file <- args[3]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read sample names
sample_names <- readLines(sample_names_file)

# Get input file paths
input_files <- file.path(input_dir, paste0(sample_names, ".fastq"))

# Check that all input files exist
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(paste("Missing input files:", paste(missing_files, collapse = ", ")))
}

cat("Processing", length(input_files), "files with DADA2\\n")

# Learn error rates
cat("Learning error rates...\\n")
err <- learnErrors(input_files, multithread = TRUE, verbose = TRUE)

# Plot error rates (if requested)
error_plot_file <- file.path(output_dir, "error_rates.pdf")
pdf(error_plot_file)
plotErrors(err, nominalQ = TRUE)
dev.off()

# Dereplicate sequences
cat("Dereplicating sequences...\\n")
derep <- derepFastq(input_files, verbose = TRUE)
names(derep) <- sample_names

# Run DADA2 denoising algorithm
cat("Running DADA2 denoising algorithm...\\n")
dada_results <- dada(derep, err = err, multithread = TRUE, verbose = TRUE)

# Make sequence table
cat("Creating sequence table...\\n")
seqtab <- makeSequenceTable(dada_results)

cat("Sequence table dimensions:", dim(seqtab), "\\n")
cat("Length distribution of sequences:\\n")
print(table(nchar(getSequences(seqtab))))

# Remove chimeras
cat("Removing chimeras...\\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                   multithread = TRUE, verbose = TRUE)

# Calculate chimera statistics
chimera_rate <- 1 - (sum(seqtab_nochim) / sum(seqtab))
cat("Chimera removal rate:", round(chimera_rate * 100, 2), "%\\n")

# Save results
cat("Saving results...\\n")

# Save sequence table as RDS
saveRDS(seqtab_nochim, file.path(output_dir, "seqtab_nochim.rds"))

# Save ASV sequences as FASTA
asv_seqs <- getSequences(seqtab_nochim)
asv_headers <- paste0("ASV_", seq_along(asv_seqs))
names(asv_seqs) <- asv_headers

writeXStringSet(DNAStringSet(asv_seqs), 
                file.path(output_dir, "asvs.fasta"))

# Save ASV table as TSV
asv_table <- seqtab_nochim
colnames(asv_table) <- asv_headers
write.table(asv_table, file.path(output_dir, "asv_table.tsv"), 
           sep = "\\t", quote = FALSE, row.names = TRUE)

# Generate summary statistics
summary_stats <- data.frame(
  sample = sample_names,
  input_reads = sapply(derep, function(x) sum(x$uniques)),
  denoised_reads = sapply(dada_results, function(x) sum(x$denoised)),
  merged_reads = rowSums(seqtab),
  nonchim_reads = rowSums(seqtab_nochim),
  stringsAsFactors = FALSE
)

# Calculate processing rates
summary_stats$denoising_rate <- summary_stats$denoised_reads / summary_stats$input_reads
summary_stats$chimera_filter_rate <- summary_stats$nonchim_reads / summary_stats$merged_reads

write.table(summary_stats, file.path(output_dir, "processing_summary.tsv"),
           sep = "\\t", quote = FALSE, row.names = FALSE)

# Print final summary
cat("\\n=== DADA2 Processing Summary ===\\n")
cat("Total samples processed:", nrow(summary_stats), "\\n")
cat("Total ASVs generated:", ncol(seqtab_nochim), "\\n")
cat("Total reads retained:", sum(summary_stats$nonchim_reads), "\\n")
cat("Average processing rate:", round(mean(summary_stats$nonchim_reads / summary_stats$input_reads) * 100, 2), "%\\n")

cat("DADA2 processing completed successfully!\\n")