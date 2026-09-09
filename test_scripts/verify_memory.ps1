# Compare IQ-TREE 3 memory usage against IQ-TREE 2 baseline + threshold.
# When IQ-TREE 2 baseline is 0 (unsupported command), falls back to the
# pre-defined expected value from expected_memory.tsv if a platform column is given.
#
# Args: $IQTree2Log    = IQ-TREE 2 log file
#       $IQTree3Log    = IQ-TREE 3 log file
#       $FallbackColumn = platform column name in expected_memory.tsv for fallback (optional)
param (
    [string] $IQTree2Log     = "time_log_iqtree2.tsv",
    [string] $IQTree3Log     = "time_log_iqtree3.tsv",
    [string] $FallbackColumn = ""
)

$WD = "test_scripts/test_data"
$thresholdFile = Join-Path $WD "expected_memory.tsv"

# Read thresholds (command + diff-threshold only)
$thresholdLines = Get-Content $thresholdFile | Select-Object -Skip 1
$thresholds = foreach ($line in $thresholdLines) {
    $parts = $line -split "`t"
    [PSCustomObject]@{ Command = $parts[0]; Threshold = [double]$parts[1] }
}

# Per-platform threshold column "thr-<platform>" when present, else diff-threshold.
# Peak memory is not comparable across platforms, so one shared allowance is either
# too tight on the noisy ones or meaningless on the quiet ones.
if ($FallbackColumn -ne "") {
    $hdr = (Get-Content $thresholdFile -TotalCount 1) -split "`t"
    $thrIdx = $hdr.IndexOf("thr-$FallbackColumn")
    if ($thrIdx -ge 0) {
        Write-Host "Using per-platform thresholds: thr-$FallbackColumn"
        for ($i = 0; $i -lt [Math]::Min($thresholds.Count, $nLog); $i++) {
            $thresholds[$i].Threshold = [double]($thresholdLines[$i] -split "`t")[$thrIdx]
        }
    } else {
        Write-Host "No thr-$FallbackColumn column; using the shared diff-threshold"
    }
}

# Resolve fallback column index
$fallbackValues = @()
if ($FallbackColumn -ne "") {
    $header = (Get-Content $thresholdFile -TotalCount 1) -split "`t"
    $colIdx = $header.IndexOf($FallbackColumn)
    if ($colIdx -lt 0) {
        Write-Host "WARNING: fallback column '$FallbackColumn' not found in $thresholdFile; skipping fallback"
        $FallbackColumn = ""
    } else {
        $fallbackValues = foreach ($line in $thresholdLines) { [double]($line -split "`t")[$colIdx] }
    }
}

# Read memory (column index 2 = PeakMemory)
$iqtree2Lines = Get-Content $IQTree2Log | Select-Object -Skip 1
$iqtree2Mem = foreach ($line in $iqtree2Lines) { [double]($line -split "`t")[2] }

$iqtree3Lines = Get-Content $IQTree3Log | Select-Object -Skip 1
$iqtree3Mem = foreach ($line in $iqtree3Lines) { [double]($line -split "`t")[2] }

# Column 0 is the command actually executed, kept so a breaching check can be retried.
$iqtree2Cmd = foreach ($line in $iqtree2Lines) { ($line -split "`t")[0] }
$iqtree3Cmd = foreach ($line in $iqtree3Lines) { ($line -split "`t")[0] }
. (Join-Path $PSScriptRoot "remeasure.ps1")

# Reconcile the number of benchmark commands with the number of table rows.
# They are joined POSITIONALLY, so a mismatch means the pairing is wrong.
$nRows = $thresholds.Count
$nLog  = $iqtree3Mem.Count
if ($nLog -ne $nRows) {
    Write-Host "WARNING: the suite ran $nLog commands but memory has $nRows threshold rows."
    if ($nLog -gt $nRows) {
        Write-Host "   Skipping the last $($nLog - $nRows) command(s) - they have no threshold:"
        $iqtree3Cmd[$nRows..($nLog - 1)] | ForEach-Object { Write-Host "     $_" }
    } else {
        Write-Host "   Ignoring the last $($nRows - $nLog) threshold row(s) - no command produced them."
        $thresholds = $thresholds[0..($nLog - 1)]
    }
    Write-Host "   NOTE: rows are matched by POSITION. If the extra command(s) were added in the"
    Write-Host "   middle rather than at the end, every later row is now compared against the"
    Write-Host "   wrong command. Add the missing row(s) to keep the table in step."
}

$failCount = 0

for ($i = 0; $i -lt [Math]::Min($thresholds.Count, $nLog); $i++) {
    $command   = $thresholds[$i].Command
    $threshold = $thresholds[$i].Threshold
    $expected  = $iqtree2Mem[$i]
    $reported  = $iqtree3Mem[$i]

    if ($expected -eq 0) {
        if ($FallbackColumn -ne "" -and $fallbackValues.Count -gt $i) {
            $expected = $fallbackValues[$i]
            Write-Host "ℹ️  ${command}: IQ-TREE 2 baseline unavailable, using pre-defined expected value (${expected}MB)"
        } else {
            Write-Host "⏭ $command skipped (IQ-TREE 2 baseline unavailable, no fallback column provided)"
            continue
        }
    }

    $allowed = $expected + $threshold
    $diff    = $reported - $expected

    # Retry once before failing: re-run this one command for both binaries and
    # re-evaluate. Costs nothing when everything passes.
    if ($reported -gt $allowed -and $iqtree3Cmd.Count -gt $i) {
        Write-Host "↻ $command exceeded (${diff}MB); retrying this command once..."
        $r2 = Measure-Once $iqtree2Cmd[$i]
        $r3 = Measure-Once $iqtree3Cmd[$i]
        if ($r2.Ok -and $r3.Ok) {
            $expected = $r2.Mem
            $reported = $r3.Mem
            $allowed  = $expected + $threshold
            $diff     = $reported - $expected
            Write-Host "   retry: IQ-TREE2 $($r2.Mem)MB, IQ-TREE3 $($r3.Mem)MB, Diff ${diff}MB"
        } else {
            Write-Host "   retry did not produce a usable measurement; keeping the first result"
        }
    }

    if ($reported -gt $allowed) {
        Write-Host "❌ $command exceeded the allowed memory usage."
        Write-Host "   Expected: ${expected}MB, Threshold: ${threshold}MB, IQ-TREE3: ${reported}MB, Diff: ${diff}MB"
        $failCount++
    } else {
        Write-Host "✅ $command passed the memory check."
        Write-Host "   Expected: ${expected}MB, Threshold: ${threshold}MB, IQ-TREE3: ${reported}MB, Diff: ${diff}MB"
    }
}

Write-Host ""

if ($failCount -eq 0) {
    Write-Host "✅ All memory checks passed."
    exit 0
} else {
    Write-Host "❌ $failCount checks failed."
    exit 1
}
