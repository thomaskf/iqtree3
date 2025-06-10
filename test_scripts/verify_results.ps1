param (
    [string]$input_file = "test_scripts/test_data/expect_ans.txt"
)

$fail_count = 0
$line_num = 0

Get-Content $input_file | ForEach-Object {
    $line = $_.Trim()
    $line_num++

    if ($line_num -eq 1 -or $line -eq "") {
        return  # skip header or empty line
    }

    $columns = $line -split "`t"
    if ($columns.Count -lt 4) {
        Write-Host "Skipping malformed line: $line"
        return
    }

    $iqtree_file = $columns[0]
    $field_name = $columns[1]
    $expected_value = [double]$columns[2]
    $threshold = [double]$columns[3]

    if (-not (Test-Path $iqtree_file)) {
        Write-Host "File not found: ${iqtree_file}"
        return
    }

    $actual_line = Select-String -Path $iqtree_file -Pattern ([regex]::Escape($field_name))
    if (-not $actual_line) {
        Write-Host "Field not found in ${iqtree_file}: ${field_name}"
        return
    }

    $match = [regex]::Match($actual_line.Line, '[-]?\d+(\.\d+)?([eE][-+]?\d+)?')
    if (-not $match.Success) {
        Write-Host "No numeric value found in line: $($actual_line.Line)"
        return
    }

    $actual_value = [double]$match.Value
    $abs_diff = [math]::Abs($actual_value - $expected_value)

    if ($abs_diff -le $threshold) {
        Write-Host "PASS: ${iqtree_file} (${field_name})" -ForegroundColor Green
    } else {
        Write-Host "FAIL: ${iqtree_file} (${field_name})" -ForegroundColor Red
        Write-Host "  Expected: ${expected_value}, Found: ${actual_value}, Diff: ${abs_diff}, Threshold: ${threshold}"
        $fail_count++
    }
}

Write-Host ""
if ($fail_count -eq 0) {
    Write-Host "✅ All checks passed." -ForegroundColor Green
} else {
    Write-Host "❌ ${fail_count} checks failed." -ForegroundColor Red
    exit 1
}
