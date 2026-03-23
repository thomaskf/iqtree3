#!/usr/bin/env Rscript
# =============================================================================
# simulate_partitions.R
#
# Simulate multi-partition alignments using IQ-TREE AliSim for benchmarking
# the new thread-scheduling strategy in PartitionFinder.
#
# Each scenario varies the number of partitions and their sizes (nptn * nstate)
# to exercise different regions of the scheduler's behaviour:
#
#   Scenario A – Many uniform light partitions   (16 x 500 bp DNA)
#   Scenario B – Mixed heavy + light partitions  (1 x 8000 bp + 15 x 200 bp DNA)
#   Scenario C – Few heavy partitions            (4 x 4000 bp DNA)
#   Scenario D – Many heavy partitions           (8 x 2000 bp DNA)
#   Scenario E – Very many very light partitions (32 x 100 bp DNA)
#
# The thread-cap formula: m_p = max(1, nptn * nstate / 4000)
# For DNA (4 states):
#   100  bp -> m_p ~ 1      500 bp -> m_p ~ 1
#   1000 bp -> m_p = 1      2000 bp -> m_p = 2
#   4000 bp -> m_p = 4      8000 bp -> m_p = 8
#
# Usage:
#   Rscript simulate_partitions.R [iqtree_bin] [outdir]
# =============================================================================

args      <- commandArgs(trailingOnly = TRUE)
IQTREE    <- if (length(args) >= 1) args[1] else "iqtree3"
OUTDIR    <- if (length(args) >= 2) args[2] else "sim_partitions"
NTAXA     <- 20
SEED      <- 12345
DNA_MODEL <- "GTR+G"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Verify IQ-TREE binary
if (system(paste(IQTREE, "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) {
  stop("IQ-TREE binary not found: ", IQTREE)
}

# -----------------------------------------------------------------------------
# Scenarios: list of named vectors giving per-partition lengths (in bp)
# -----------------------------------------------------------------------------
scenarios <- list(
  A_uniform_light  = rep(500,  16),   # 16 partitions x 500 bp
  B_mixed          = c(8000, rep(200, 15)),  # 1 heavy + 15 light
  C_few_heavy      = rep(4000,  4),   # 4 partitions x 4000 bp
  D_many_heavy     = rep(2000,  8),   # 8 partitions x 2000 bp
  E_many_very_light = rep(100,  32)   # 32 partitions x 100 bp
)

# -----------------------------------------------------------------------------
# Helper: write a NEXUS partition file
# lengths: named integer vector, names = partition names, values = bp lengths
# model:   substitution model string (same for all partitions)
# -----------------------------------------------------------------------------
write_partition_nexus <- function(nexus_path, lengths, model) {
  n   <- length(lengths)
  end <- cumsum(lengths)
  beg <- c(1L, end[-n] + 1L)
  pnames <- paste0("part", seq_len(n))

  lines <- c(
    "#nexus",
    "begin sets;"
  )
  for (i in seq_len(n)) {
    lines <- c(lines, sprintf("  charset %s = %d-%d;", pnames[i], beg[i], end[i]))
  }
  # charpartition line assigns the model to each charset
  cp <- paste(sprintf("%s:%s", model, pnames), collapse = ", ")
  lines <- c(lines, sprintf("  charpartition myparts = %s;", cp))
  lines <- c(lines, "end;")

  writeLines(lines, nexus_path)
}

# -----------------------------------------------------------------------------
# Helper: run AliSim for one scenario
# Returns the output .phy path, or NULL on failure
# -----------------------------------------------------------------------------
run_alisim <- function(scenario_name, lengths, model, ntaxa, seed, outdir, iqtree) {
  prefix    <- file.path(outdir, scenario_name)
  nex_path  <- paste0(prefix, "_partition.nex")
  total_len <- sum(lengths)
  n_parts   <- length(lengths)

  write_partition_nexus(nex_path, lengths, model)

  cmd <- sprintf(
    '%s --alisim "%s" -p "%s" -t "RANDOM{yh/%d}" --seed %d --redo -quiet 2>/dev/null',
    iqtree, prefix, nex_path, ntaxa, seed
  )

  status <- system(cmd)
  phy    <- paste0(prefix, ".phy")

  if (status == 0 && file.exists(phy)) {
    cat(sprintf("  OK   %-25s  %2d partitions  total %6d bp\n",
                scenario_name, n_parts, total_len))
    return(phy)
  } else {
    cat(sprintf("  FAIL %-25s\n", scenario_name))
    return(NULL)
  }
}

# -----------------------------------------------------------------------------
# Run all scenarios
# -----------------------------------------------------------------------------
cat("IQ-TREE binary :", IQTREE, "\n")
cat("Output dir     :", OUTDIR, "\n")
cat("Taxa           :", NTAXA,  "  Seed:", SEED, "\n\n")
cat("=== Simulating partition alignments ===\n")

results <- list()
for (sname in names(scenarios)) {
  phy <- run_alisim(sname, scenarios[[sname]], DNA_MODEL,
                   NTAXA, SEED, OUTDIR, IQTREE)
  results[[sname]] <- list(
    phy       = phy,
    n_parts   = length(scenarios[[sname]]),
    lengths   = scenarios[[sname]],
    total_len = sum(scenarios[[sname]])
  )
}

# Save scenario metadata so the benchmark script can read it
meta <- do.call(rbind, lapply(names(results), function(s) {
  r <- results[[s]]
  data.frame(
    scenario  = s,
    n_parts   = r$n_parts,
    total_len = r$total_len,
    min_len   = min(r$lengths),
    max_len   = max(r$lengths),
    phy       = if (!is.null(r$phy)) r$phy else NA_character_,
    stringsAsFactors = FALSE
  )
}))
meta_path <- file.path(OUTDIR, "scenarios.csv")
write.csv(meta, meta_path, row.names = FALSE)

cat("\nScenario summary:\n")
print(meta[, c("scenario", "n_parts", "total_len", "min_len", "max_len")])
cat("\nMetadata saved to:", meta_path, "\n")
