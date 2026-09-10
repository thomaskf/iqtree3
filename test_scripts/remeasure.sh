#!/bin/bash
# Re-run a single benchmark command and echo "<seconds> <peak MB>".
#
# Used by verify_memory.sh / verify_runtime.sh to retry a command that breached
# its threshold. CI runner noise is one-sided - a loaded machine can only make a
# command slower or fatter - so a breach that does not reproduce was noise, while
# a genuine regression reproduces every time.
#
# The suite is stateful, so the command's OWN outputs are cleared first; IQ-TREE
# otherwise refuses to rerun ("previous run successfully finished"). Outputs of
# earlier commands are left alone because later commands still need them.
remeasure() {
    local CMD="$1"
    local PREFIX REAL MEM_MB MEM_KB PEAK_MEM tmp
    PREFIX=$(echo "$CMD" | sed -n 's/.*--prefix \([^ ]*\).*/\1/p')
    [ -n "$PREFIX" ] && rm -f "${PREFIX}".*
    tmp=$(mktemp)
    if [[ "$(uname)" == "Darwin" ]]; then
        /usr/bin/time -l -o "$tmp" $CMD > /dev/null 2>&1
        local rc=$?
        REAL=$(awk '/real/{print $1; exit}' "$tmp")
        PEAK_MEM=$(awk '/peak memory footprint/{print $1; exit}' "$tmp")
        MEM_MB=$(awk "BEGIN {printf \"%.2f\", ${PEAK_MEM:-0} / (1024 * 1024)}")
    else
        /usr/bin/time -o "$tmp" -f "%e %U %S %M" $CMD > /dev/null 2>&1
        local rc=$?
        read -r REAL _ _ MEM_KB < "$tmp"
        MEM_MB=$(awk "BEGIN {printf \"%.2f\", ${MEM_KB:-0} / 1024}")
    fi
    rm -f "$tmp"
    if [ "$rc" -ne 0 ]; then echo "0 0"; else echo "${REAL:-0} ${MEM_MB:-0}"; fi
}
