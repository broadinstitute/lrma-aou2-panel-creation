#!/usr/bin/env bash
#
# Resource-contract tests for GLIMPSE2Phase argument handling and sizing.
#
# GLIMPSE2Phase computes its memory request from phase_threads and phase_kpbwt, so both must
# reach GLIMPSE2 exactly once and exactly as typed. Two things make that non-trivial:
#
#   * callers pass a free-form extra_phase_args string, and the historical default (plus the
#     staged input CSVs) carried --thread and --Kpbwt inside it;
#   * GLIMPSE2 catches boost program_options errors and exits 0, so a duplicate option does
#     not fail loudly -- it surfaces later as a missing output BCF.
#
# Runs locally in a second; no cloud, no Cromwell, no GLIMPSE2 binary required.
#
# Usage: bash test/wdl/test-phase-arg-handling.sh

set -uo pipefail

WDL="$(dirname "$0")/../../wdl/methods/imputation/GLIMPSE2FromPreprocessedPLsJoint.wdl"
pass=0
fail=0

ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass + 1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n     %s\n' "$1" "${2:-}"; fail=$((fail + 1)); }

# ---------------------------------------------------------------------------
# Mirror of the WDL command block. Keep in sync with GLIMPSE2Phase.
# ---------------------------------------------------------------------------
build_cmd() {
    local phase_threads="$1" phase_kpbwt="$2" EXTRA_PHASE_ARGS="$3" OPT

    if [ "$phase_threads" -lt 1 ] || [ "$phase_kpbwt" -lt 1 ]; then
        echo "GUARD-EXIT"
        return 0
    fi

    for OPT in thread Kpbwt; do
        if echo "$EXTRA_PHASE_ARGS" | grep -qE "(^|[[:space:]])--${OPT}([[:space:]]|=|$)"; then
            EXTRA_PHASE_ARGS=$(echo "$EXTRA_PHASE_ARGS" \
                | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]+|=)[^-[:space:]][^[:space:]]*/ /g" \
                | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]|=|$)/ /g")
        fi
    done

    echo "--thread ${phase_threads} --Kpbwt ${phase_kpbwt} ${EXTRA_PHASE_ARGS}"
}

# Count whole-token occurrences, so --threads is not mistaken for --thread.
count_opt() { echo " $1 " | grep -oE "[[:space:]]--$2([[:space:]]|=)" | wc -l | tr -d ' '; }

expect_single() {
    local desc="$1" out
    out=$(build_cmd "$2" "$3" "$4")
    local nt nk
    nt=$(count_opt "$out" thread)
    nk=$(count_opt "$out" Kpbwt)
    if [ "$nt" -eq 1 ] && [ "$nk" -eq 1 ]; then ok "$desc"
    else bad "$desc" "thread=$nt kpbwt=$nk :: $out"; fi
}

expect_guard() {
    local desc="$1" out
    out=$(build_cmd "$2" "$3" "$4")
    [ "$out" = "GUARD-EXIT" ] && ok "$desc" || bad "$desc" "expected GUARD-EXIT, got: $out"
}

expect_contains() {
    local desc="$1" out
    out=$(build_cmd 4 1000 "$2")
    case "$out" in *"$3"*) ok "$desc";; *) bad "$desc" "missing '$3' in: $out";; esac
}

DEFAULT="--impute-reference-only-variants --keep-monomorphic-ref-sites --main 10 --burnin 5 --err-imp 1E-3"
LEGACY="--thread \$(nproc) --impute-reference-only-variants --keep-monomorphic-ref-sites --Kpbwt 1000 --main 10 --burnin 5 --err-imp 1E-3"
STAGED_CSV="--thread 4 --impute-reference-only-variants --keep-monomorphic-ref-sites --Kpbwt 1000 --main 10 --burnin 5 --err-imp 1E-3"

echo
echo "GLIMPSE2Phase receives exactly one --thread and one --Kpbwt"
expect_single "current WDL default"                     4 1000 "$DEFAULT"
expect_single "staged input CSV (--thread 4, --Kpbwt)"  4 1000 "$STAGED_CSV"
expect_single "legacy default (--thread \$(nproc))"     4 1000 "$LEGACY"
expect_single "equals form (--thread=8 --Kpbwt=2000)"   4 1000 "--thread=8 --Kpbwt=2000 --main 10"
expect_single "sibling WDL default (no --Kpbwt)"        4 1000 "--impute-reference-only-variants --keep-monomorphic-ref-sites"
expect_single "bcftools-style --threads left alone"     4 1000 "--threads 4 --main 10"
expect_single "valueless trailing --thread"             4 1000 "--main 10 --thread"
expect_single "reordered (--Kpbwt first, --thread last)" 4 1000 "--Kpbwt 2000 --main 10 --thread 16"
expect_single "empty extra_phase_args"                  4 1000 ""

echo
echo "Invalid sizing parameters fail fast rather than reaching GLIMPSE2"
expect_guard "phase_threads = 0"    0 1000 "$DEFAULT"
expect_guard "phase_kpbwt = 0"      4 0    "$DEFAULT"
expect_guard "phase_threads = -1"  -1 1000 "$DEFAULT"

echo
echo "Unrelated options survive stripping"
for opt in --impute-reference-only-variants --keep-monomorphic-ref-sites "--main 10" "--burnin 5" "--err-imp 1E-3"; do
    expect_contains "preserved: $opt" "$STAGED_CSV" "$opt"
done
expect_contains "preserved: --threads 4 (not --thread)" "--threads 4 --main 10" "--threads 4"

# ---------------------------------------------------------------------------
# Sizing: must reproduce the empirically verified shape at the defaults, and must
# never request more than the N1 limit of 6.5 GB per cpu.
# ---------------------------------------------------------------------------
echo
echo "Memory/CPU sizing"
sizing() {  # kpbwt threads -> "mem cpu"
    python3 -c '
import math, sys
kp, t = int(sys.argv[1]), int(sys.argv[2])
mem = 8 + math.ceil(32.0 * kp * t / 4000.0)
r = math.ceil(mem / 6.5)
cpu = r + (r % 2) if r > t else t
print(mem, cpu)' "$1" "$2"
}

read -r mem cpu <<<"$(sizing 1000 4)"
if [ "$mem" -eq 40 ] && [ "$cpu" -eq 8 ]; then
    ok "defaults reproduce the verified 40 GiB / 8 cpu"
else
    bad "defaults reproduce the verified 40 GiB / 8 cpu" "got ${mem} GiB / ${cpu} cpu"
fi

for combo in "1000 4" "1000 8" "2000 4" "2000 8" "500 4" "4000 4"; do
    set -- $combo
    read -r mem cpu <<<"$(sizing "$1" "$2")"
    within=$(python3 -c "print(1 if $mem/$cpu <= 6.5 else 0)")
    [ "$within" -eq 1 ] \
        && ok "Kpbwt=$1 threads=$2 -> ${mem} GiB / ${cpu} cpu within 6.5 GB/cpu" \
        || bad "Kpbwt=$1 threads=$2" "${mem}/${cpu} exceeds 6.5 GB/cpu"
done

# ---------------------------------------------------------------------------
# The WDL itself must not reintroduce the duplication this guards against.
# ---------------------------------------------------------------------------
echo
echo "WDL source invariants"
if [ -f "$WDL" ]; then
    grep -q 'String extra_phase_args = "--impute-reference-only-variants' "$WDL" \
        && ok "workflow default carries no --thread/--Kpbwt" \
        || bad "workflow default carries no --thread/--Kpbwt" "default string changed"

    if grep -E 'extra_phase_args\s*=\s*"' "$WDL" | grep -qE '\-\-(thread|Kpbwt)'; then
        bad "no default string reintroduces --thread/--Kpbwt" "found one"
    else
        ok "no default string reintroduces --thread/--Kpbwt"
    fi
else
    bad "WDL present at expected path" "$WDL"
fi

echo
echo "-------------------------------------------"
printf 'passed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ] || exit 1
