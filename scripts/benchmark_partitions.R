#!/usr/bin/env Rscript
# =============================================================================
# benchmark_partitions.R
#
# Benchmark the original and modified IQ-TREE on multi-partition alignments,
# measuring wall-clock time for ModelFinder (-m TEST) across thread counts.
# Produces a TSV of raw results and diagnostic plots.
#
# Usage:
#   Rscript benchmark_partitions.R [IQTREE_NEW] [IQTREE_ORIG] [SIMDIR] [OUTDIR]
#
# Examples:
#   Rscript benchmark_partitions.R ./iqtree3_new ./iqtree3_orig sim_partitions bench_partitions
# =============================================================================

args        <- commandArgs(trailingOnly = TRUE)
IQTREE_NEW  <- if (length(args) >= 1) args[1] else "iqtree3"
IQTREE_ORIG <- if (length(args) >= 2) args[2] else NULL   # optional
SIMDIR      <- if (length(args) >= 3) args[3] else "sim_partitions"
OUTDIR      <- if (length(args) >= 4) args[4] else "bench_partitions"

THREAD_COUNTS <- c(1, 2, 4, 8, 16)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Verify binaries
check_binary <- function(bin) {
  if (system(paste(bin, "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0)
    stop("IQ-TREE binary not found: ", bin)
}
check_binary(IQTREE_NEW)
if (!is.null(IQTREE_ORIG)) check_binary(IQTREE_ORIG)

# Load scenario metadata written by simulate_partitions.R
meta_path <- file.path(SIMDIR, "scenarios.csv")
if (!file.exists(meta_path))
  stop("Scenario metadata not found: ", meta_path,
       "\nRun simulate_partitions.R first.")
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
meta <- meta[!is.na(meta$phy) & file.exists(meta$phy), ]
if (nrow(meta) == 0) stop("No valid simulation files found in ", SIMDIR)

cat("Modified binary :", IQTREE_NEW, "\n")
cat("Original binary :", if (!is.null(IQTREE_ORIG)) IQTREE_ORIG else "(not provided)", "\n")
cat("Simulation dir  :", SIMDIR, "\n")
cat("Results dir     :", OUTDIR, "\n\n")

# -----------------------------------------------------------------------------
# Helper: time one ModelFinder run, return elapsed seconds or NA on failure
# -----------------------------------------------------------------------------
time_run <- function(iqtree, aln, nex, nt, prefix) {
  cmd <- sprintf(
    '%s -s "%s" -p "%s" -m TEST -nt %d --prefix "%s" --redo -quiet 2>/dev/null',
    iqtree, aln, nex, nt, prefix
  )
  t0     <- proc.time()[["elapsed"]]
  status <- system(cmd)
  t1     <- proc.time()[["elapsed"]]
  if (status == 0) round(t1 - t0, 3) else NA_real_
}

# Build list of (binary_label, binary_path) pairs to benchmark
binaries <- list(modified = IQTREE_NEW)
if (!is.null(IQTREE_ORIG)) binaries$original <- IQTREE_ORIG

# -----------------------------------------------------------------------------
# Main benchmark loop
# -----------------------------------------------------------------------------
rows <- list()

for (i in seq_len(nrow(meta))) {
  s       <- meta$scenario[i]
  phy     <- meta$phy[i]
  n_parts <- meta$n_parts[i]
  nex     <- sub("\\.phy$", "_partition.nex", phy)

  if (!file.exists(nex)) {
    cat("  SKIP", s, "(partition file not found)\n")
    next
  }

  cat("---", s, " (", n_parts, "partitions,", meta$total_len[i], "bp) ---\n")

  for (bname in names(binaries)) {
    bin <- binaries[[bname]]
    for (nt in THREAD_COUNTS) {
      prefix <- file.path(OUTDIR, sprintf("%s_%s_nt%d", s, bname, nt))
      cat(sprintf("  [%-8s] -nt %2d ... ", bname, nt))
      elapsed <- time_run(bin, phy, nex, nt, prefix)
      if (!is.na(elapsed)) cat(elapsed, "s\n") else cat("FAILED\n")
      rows[[length(rows) + 1]] <- data.frame(
        version   = bname,
        scenario  = s,
        n_parts   = n_parts,
        total_len = meta$total_len[i],
        min_len   = meta$min_len[i],
        max_len   = meta$max_len[i],
        num_threads = nt,
        wall_sec  = elapsed,
        stringsAsFactors = FALSE
      )
    }
  }
  cat("\n")
}

results <- do.call(rbind, rows)
tsv_path <- file.path(OUTDIR, "results.tsv")
write.table(results, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Results saved to:", tsv_path, "\n\n")

# -----------------------------------------------------------------------------
# Plots (requires ggplot2)
# -----------------------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  cat("ggplot2 not available — skipping plots.\n")
  quit(save = "no")
}
library(ggplot2)

results_ok <- results[!is.na(results$wall_sec), ]

# -- Plot 1: Wall time vs threads, faceted by scenario, coloured by version --
p1 <- ggplot(results_ok,
             aes(x = num_threads, y = wall_sec,
                 colour = version, shape = version, group = version)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ scenario, scales = "free_y") +
  scale_x_continuous(breaks = THREAD_COUNTS) +
  labs(
    title   = "ModelFinder wall time: modified vs original IQ-TREE (partition data)",
    x       = "Number of threads (-nt)",
    y       = "Wall time (seconds)",
    colour  = "Version",
    shape   = "Version"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "walltime_by_scenario.pdf"), p1, width = 12, height = 8)
ggsave(file.path(OUTDIR, "walltime_by_scenario.png"), p1, width = 12, height = 8, dpi = 150)
cat("Plot 1 saved: walltime_by_scenario.{pdf,png}\n")

# -- Plot 2: Speedup relative to single-thread baseline, per version --
baseline <- results_ok[results_ok$num_threads == 1, c("version", "scenario", "wall_sec")]
names(baseline)[3] <- "wall_1thread"
results_sp <- merge(results_ok, baseline, by = c("version", "scenario"))
results_sp$speedup <- results_sp$wall_1thread / results_sp$wall_sec

p2 <- ggplot(results_sp,
             aes(x = num_threads, y = speedup,
                 colour = version, shape = version, group = version)) +
  geom_line() +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ scenario) +
  scale_x_continuous(breaks = THREAD_COUNTS) +
  labs(
    title   = "Speedup over single-thread baseline (partition data)",
    x       = "Number of threads (-nt)",
    y       = "Speedup (T1 / Tk)",
    colour  = "Version",
    shape   = "Version",
    caption = "Dashed line = ideal linear speedup"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "speedup_by_scenario.pdf"), p2, width = 12, height = 8)
ggsave(file.path(OUTDIR, "speedup_by_scenario.png"), p2, width = 12, height = 8, dpi = 150)
cat("Plot 2 saved: speedup_by_scenario.{pdf,png}\n")

# -- Plot 3: Direct modified vs original comparison (ratio) --
if (!is.null(IQTREE_ORIG)) {
  mod_t  <- results_ok[results_ok$version == "modified",  c("scenario", "num_threads", "wall_sec")]
  orig_t <- results_ok[results_ok$version == "original",  c("scenario", "num_threads", "wall_sec")]
  names(mod_t)[3]  <- "wall_modified"
  names(orig_t)[3] <- "wall_original"
  comp <- merge(mod_t, orig_t, by = c("scenario", "num_threads"))
  comp$ratio <- comp$wall_original / comp$wall_modified  # >1 means modified is faster

  p3 <- ggplot(comp,
               aes(x = num_threads, y = ratio, group = scenario, colour = scenario)) +
    geom_line() +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = THREAD_COUNTS) +
    labs(
      title   = "Speedup of modified over original IQ-TREE (partition data)",
      x       = "Number of threads (-nt)",
      y       = "Time ratio (original / modified)",
      colour  = "Scenario",
      caption = "Ratio > 1: modified is faster; ratio < 1: original is faster"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTDIR, "ratio_modified_vs_original.pdf"), p3, width = 10, height = 6)
  ggsave(file.path(OUTDIR, "ratio_modified_vs_original.png"), p3, width = 10, height = 6, dpi = 150)
  cat("Plot 3 saved: ratio_modified_vs_original.{pdf,png}\n")

  # Print summary table
  cat("\nSummary: original / modified wall time ratio (>1 = modified is faster)\n")
  comp_wide <- reshape(comp[, c("scenario", "num_threads", "ratio")],
                       idvar = "scenario", timevar = "num_threads", direction = "wide")
  names(comp_wide) <- sub("ratio\\.", "nt", names(comp_wide))
  print(comp_wide, row.names = FALSE, digits = 3)
}
