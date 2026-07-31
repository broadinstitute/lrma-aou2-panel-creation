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
expect_single "legacy default, literal \$(nproc)"       4 1000 "$LEGACY"
# The WDL interpolates extra_phase_args into the script as text, so bash expands $(nproc)
# at assignment; the strip logic therefore sees an already-substituted integer. Cover both.
expect_single "legacy default, pre-expanded by bash"    4 1000 "--thread 8 --impute-reference-only-variants --Kpbwt 1000 --main 10"
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
# Per-shard sizing: mem = 8 + ceil(4 * threads * Kpbwt * L / 1e9), cpu = even(max(threads,
# ceil(mem/6.5))). Must clear the observed OOM boundary, stay under the N1 6.5 GB/cpu limit,
# and stay even.
# ---------------------------------------------------------------------------
echo
echo "Memory/CPU sizing"
sizing() {  # kpbwt threads L -> "mem cpu"
    python3 -c '
import math, sys
kp, t, L = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
mem = 8 + math.ceil(4.0 * t * kp * L / 1e9)
r = math.ceil(mem / 6.5)
u = r if r > t else t
print(mem, u + (u % 2))' "$1" "$2" "$3"
}

# Every shard observed to OOM at 16 GiB must now be sized above 16; the one that passed at 16
# must still be at least what it needed.
while read -r name L outcome floor; do
    read -r mem cpu <<<"$(sizing 1000 4 "$L")"
    if [ "$mem" -ge "$floor" ]; then
        ok "L=$L ($name, $outcome) -> ${mem} GiB >= ${floor}"
    else
        bad "L=$L ($name, $outcome)" "got ${mem} GiB, need >= ${floor}"
    fi
done <<'CASES'
chr22_s0 657784 passed_marginally_at_16 17
chr3_s16 765770 OOMed_at_16 17
chr20_s4 965039 OOMed_at_16 17
chr20_s5 1005104 OOMed_at_16 17
chr11_s9 1122361 OOMed_at_16 17
chr7_s12 1346888 OOMed_at_16 17
CASES

# Ratio limit and even cpu across the parameter space, including odd thread counts.
for combo in "1000 4 400000" "1000 4 1346888" "1000 8 1346888" "2000 4 1346888" \
             "500 4 400000" "1000 5 900000" "1000 1 400000" "1000 3 700000" "1000 4 100000"; do
    set -- $combo
    read -r mem cpu <<<"$(sizing "$1" "$2" "$3")"
    within=$(python3 -c "print(1 if $mem/$cpu <= 6.5 else 0)")
    even=$(python3 -c "print(1 if $cpu % 2 == 0 else 0)")
    if [ "$within" -eq 1 ] && [ "$even" -eq 1 ]; then
        ok "Kpbwt=$1 threads=$2 L=$3 -> ${mem} GiB / ${cpu} cpu (even, <=6.5 GB/cpu)"
    else
        bad "Kpbwt=$1 threads=$2 L=$3" "${mem}/${cpu} within=${within} even=${even}"
    fi
done

# Doubling Kpbwt or threads must double the matrix term, not leave the request unchanged.
read -r m1 _ <<<"$(sizing 1000 4 1000000)"
read -r m2 _ <<<"$(sizing 2000 4 1000000)"
read -r m3 _ <<<"$(sizing 1000 8 1000000)"
[ $((m2 - 8)) -eq $(( (m1 - 8) * 2 )) ] && ok "Kpbwt 1000->2000 doubles the matrix term" \
    || bad "Kpbwt scaling" "$m1 -> $m2"
[ $((m3 - 8)) -eq $(( (m1 - 8) * 2 )) ] && ok "threads 4->8 doubles the matrix term" \
    || bad "thread scaling" "$m1 -> $m3"

# Per-shard sizing must not cost more than the flat worst-case bound it replaces.
python3 - <<'COST'
import math
CPU, MEM, SPOT = 0.0332, 0.00445, 0.30
def size(L, t=4, kp=1000):
    m = 8 + math.ceil(4.0 * t * kp * L / 1e9)
    r = math.ceil(m / 6.5); u = r if r > t else t
    return m, u + (u % 2)
def cost(n, cpu, gb, mins): return n * (cpu*CPU + gb*MEM) * (mins/60) * SPOT
Ls = [400000]*490 + [657784,765770,965039,1005104,1122361,1346888] + [600000]*27
dyn = sum(cost(1, c, m, 35) for m, c in (size(L) for L in Ls))
flat = cost(len(Ls), 8, 40, 35)
print("  \033[32mPASS\033[0m  per-shard $%.1f < flat worst-case $%.1f per batch" % (dyn, flat)
      if dyn < flat else "  \033[31mFAIL\033[0m  per-shard $%.1f >= flat $%.1f" % (dyn, flat))
COST

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
