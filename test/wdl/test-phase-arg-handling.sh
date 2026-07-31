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
                | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]+|=)-?[0-9]+/ /g" \
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
expect_single "negative value (--thread -1)"            4 1000 "--thread -1 --main 10"

# A stripped negative value must not leave "-1" behind as a stray positional argument.
ORPHAN=$(build_cmd 4 1000 "--thread -1 --main 10")
case "$ORPHAN" in *" -1"*) bad "no orphan token from --thread -1" "$ORPHAN";;
                  *) ok "no orphan token from --thread -1";; esac

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

    # build_cmd above mirrors the WDL rather than executing it, so it can pass while the WDL
    # drifts. Assert the two sed patterns it mirrors are actually present, and in particular
    # that the value matcher is not narrowed back to [0-9]+ -- which is the exact regression
    # that let "--thread $(nproc)" through and produced a duplicate option.
    if grep -qF 's/(^|[[:space:]])--${OPT}([[:space:]]+|=)[^-[:space:]][^[:space:]]*/ /g' "$WDL"; then
        ok "WDL still strips non-integer option values"
    else
        bad "WDL still strips non-integer option values" "value-matching sed pattern changed"
    fi
    if grep -qF 's/(^|[[:space:]])--${OPT}([[:space:]]|=|$)/ /g' "$WDL"; then
        ok "WDL still strips valueless leftovers"
    else
        bad "WDL still strips valueless leftovers" "second sed pattern changed"
    fi
    if grep -qE '\-\-\$\{OPT\}\(\[\[:space:\]\]\+\|=\)\[0-9\]\+' "$WDL"; then
        bad "value matcher not narrowed to [0-9]+" "found the regressed pattern"
    else
        ok "value matcher not narrowed to [0-9]+"
    fi

    # The memory expression must keep Float promotion first; reordering overflows Int32.
    if grep -qF 'ceil((((4.0 * phase_threads) * phase_kpbwt) * n_variants)' "$WDL"; then
        ok "memory expression promotes to Float before multiplying"
    else
        bad "memory expression promotes to Float before multiplying" "expression reordered"
    fi

    # Instrumentation must be trap-based and on stderr, or it vanishes exactly when needed.
    if [ "$(grep -c 'trap _instr_report EXIT' "$WDL")" -eq 4 ]; then
        ok "all four measured tasks install an EXIT trap"
    else
        bad "all four measured tasks install an EXIT trap" "found $(grep -c 'trap _instr_report EXIT' "$WDL")"
    fi
    if grep -q 'wall_s=\$((SECONDS-_INSTR_T0))" >&2' "$WDL"; then
        ok "resource lines are written to stderr"
    else
        bad "resource lines are written to stderr" "redirect changed"
    fi
    if grep -qF 's/(^|[[:space:]])--${OPT}([[:space:]]+|=)-?[0-9]+/ /g' "$WDL"; then
        ok "WDL strips signed numeric option values"
    else
        bad "WDL strips signed numeric option values" "signed-value sed pattern missing"
    fi
    if [ "$(grep -c 'Int eff_cpu       = select_first(\[runtime_attr.cpu_cores' "$WDL")" -eq 4 ]; then
        ok "cpu is re-derived from effective memory in all four tasks"
    else
        bad "cpu is re-derived from effective memory in all four tasks" "partial override coupling not enforced"
    fi

    # eff_* must be declared before the command block that references them.
    if python3 - "$WDL" <<'ORDER'
import re,sys
s=open(sys.argv[1]).read(); bad=[]
for m in re.finditer(r'task (\w+) \{', s):
    seg=s[m.end():]; nx=re.search(r'\ntask ',seg); seg=seg[:nx.start()] if nx else seg
    d=seg.find('Float eff_mem_gb'); c=seg.find('command <<<')
    if d>=0 and c>=0 and d>c: bad.append(m.group(1))
sys.exit(1 if bad else 0)
ORDER
    then ok "eff_* declared before the command block in every task"
    else bad "eff_* declared before the command block in every task" "forward reference"; fi

    # df must be POSIX, not GNU-only: -BG does not exist on busybox.
    if grep -q 'df -Pk .' "$WDL" && ! grep -q 'df -P -BG' "$WDL"; then
        ok "df usage is POSIX-portable"
    else
        bad "df usage is POSIX-portable" "GNU-only -BG present"
    fi

    # Container images must be pinned by digest, not by tag: GCR tags are mutable.
    if grep -oE '"[a-z0-9.-]+/[^"]*"' "$WDL" | grep -E 'gcr\.io' | grep -qv '@sha256:'; then
        bad "all container images pinned by digest" "$(grep -oE '"[a-z0-9.-]+/[^"]*"' "$WDL" | grep 'gcr\.io' | grep -v '@sha256:' | head -1)"
    else
        ok "all container images pinned by digest"
    fi

    # The counting task serialises the chromosome; it must not be preemptible.
    if awk '/task CountPanelVariantsPerShard/,/^}/' "$WDL" | grep -qE 'preemptible_tries:\s*0,'; then
        ok "CountPanelVariantsPerShard is non-preemptible"
    else
        bad "CountPanelVariantsPerShard is non-preemptible" "preemptible_tries is not 0"
    fi
else
    bad "WDL present at expected path" "$WDL"
fi

# ---------------------------------------------------------------------------
# Resource reporting: the instrumentation must emit a parseable line, must never be able to
# fail the task it measures, and the harvester must flag both failure directions.
# ---------------------------------------------------------------------------
echo
echo "Resource reporting"
HARVEST="$(dirname "$0")/../../scripts/harvest-resources.py"

# The instrumentation must report even when the measured child is SIGKILLed, because that
# is the case worth measuring. An end-of-script report does not: set -e aborts first.
INSTR_BODY='
_instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk "{printf \"%.2f\", \$1/1073741824}" "$f" 2>/dev/null && return 0; fi; done; printf NA; }
_instr_peak()  { _instr_gib /sys/fs/cgroup/memory.peak /sys/fs/cgroup/memory/memory.max_usage_in_bytes; }
_instr_limit() { _instr_gib /sys/fs/cgroup/memory.max  /sys/fs/cgroup/memory/memory.limit_in_bytes; }
_INSTR_T0=$SECONDS
_INSTR_SAMPLER=""
_instr_report() {
    _rc=$?
    set +e +x
    [ -n "$_INSTR_SAMPLER" ] && kill "$_INSTR_SAMPLER" 2>/dev/null
    echo "[RESOURCE] task=T rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) wall_s=$((SECONDS-_INSTR_T0))" >&2
}
trap _instr_report EXIT
'
KILLED=$(bash -c "set -euxo pipefail; $INSTR_BODY bash -c 'kill -9 \$\$'; echo SHOULD_NOT_RUN" 2>&1 >/dev/null | grep '\[RESOURCE\]')
case "$KILLED" in *"rc=137"*) ok "instrumentation reports after a SIGKILL (rc=137)";;
                  *) bad "instrumentation reports after a SIGKILL (rc=137)" "$KILLED";; esac
case "$KILLED" in *SHOULD_NOT_RUN*) bad "trap does not resurrect the aborted script" "$KILLED";;
                  *) ok "trap does not resurrect the aborted script";; esac

OKRUN=$(bash -c "set -euxo pipefail; $INSTR_BODY true" 2>&1 >/dev/null | grep '\[RESOURCE\]')
case "$OKRUN" in *"rc=0"*peak_rss_gib=*limit_gib=*wall_s=*) ok "instrumentation reports a complete line on success";;
                 *) bad "instrumentation reports a complete line on success" "$OKRUN";; esac

# It must go to stderr: Cromwell delocalizes stderr for FAILED tasks, File outputs do not exist.
STDOUT_ONLY=$(bash -c "set -euxo pipefail; $INSTR_BODY true" 2>/dev/null | grep -c '\[RESOURCE\]' || true)
[ "$STDOUT_ONLY" -eq 0 ] && ok "instrumentation writes to stderr, not stdout" \
    || bad "instrumentation writes to stderr, not stdout" "found $STDOUT_ONLY line(s) on stdout"

# cgroup byte->GiB conversion, and graceful degradation when cgroup is absent
CGDIR=$(mktemp -d); echo 34359738368 > "$CGDIR/peak"
CONV=$(bash -c '_instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk "{printf \"%.2f\", \$1/1073741824}" "$f" 2>/dev/null && return 0; fi; done; printf NA; }; _instr_gib '"$CGDIR"'/peak')
[ "$CONV" = "32.00" ] && ok "cgroup bytes convert to GiB correctly" || bad "cgroup bytes convert to GiB correctly" "got $CONV"
MISS=$(bash -c '_instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk "{printf \"%.2f\", \$1/1073741824}" "$f" 2>/dev/null && return 0; fi; done; printf NA; }; _instr_gib /nope /nada')
[ "$MISS" = "NA" ] && ok "degrades to NA when cgroup is unreadable" || bad "degrades to NA when cgroup is unreadable" "got $MISS"
rm -rf "$CGDIR"

if [ -x "$HARVEST" ]; then
    # under-request must be flagged: a shard that used more than it asked for
    OVER='[RESOURCE] task=GLIMPSE2Phase n_variants=1346888 requested_mem_gib=30 requested_cpu=6 threads=4 kpbwt=1000 peak_rss_gib=31.66 peak_rss_source=cgroup-v2 wall_s=1'
    echo "$OVER" | "$HARVEST" 2>/dev/null | grep -q "OVER REQUEST" \
        && ok "harvester flags an under-sized request" \
        || bad "harvester flags an under-sized request" "no OVER REQUEST warning"
    echo "$OVER" | "$HARVEST" 2>/dev/null | grep -q "bound is NOT a bound" \
        && ok "harvester flags a breached memory bound" \
        || bad "harvester flags a breached memory bound" "no bound warning"

    # over-provisioning must be flagged too, or the data never drives sizing down
    UNDER='[RESOURCE] task=GLIMPSE2Ligate n_shards=24 requested_mem_gib=32 requested_cpu=6 peak_rss_gib=1.86 peak_rss_source=cgroup-v2 wall_s=1'
    echo "$UNDER" | "$HARVEST" 2>/dev/null | grep -q "OVER-PROVISIONED" \
        && ok "harvester flags over-provisioning" \
        || bad "harvester flags over-provisioning" "no OVER-PROVISIONED warning"

    # must tolerate raw logs and missing measurements rather than crashing
    printf 'unrelated log line\n[RESOURCE] task=X peak_rss_gib=NA peak_rss_source=unavailable wall_s=3\n' \
        | "$HARVEST" >/dev/null 2>&1 \
        && ok "harvester tolerates noise and NA measurements" \
        || bad "harvester tolerates noise and NA measurements" "non-zero exit"
else
    bad "harvest-resources.py present and executable" "$HARVEST"
fi

# ---------------------------------------------------------------------------
# Static checks on the real command blocks. Everything above tests a mirror of the WDL's
# shell; this parses the WDL itself, which is the only way to catch a call to a function
# that does not exist in that block.
# ---------------------------------------------------------------------------
echo
echo "WDL command blocks"
CHECKER="$(dirname "$0")/check-command-blocks.py"
if python3 "$CHECKER" "$WDL"; then
    ok "command blocks: defined-before-called, syntax, stderr"
else
    bad "command blocks: defined-before-called, syntax, stderr" "see above"
fi

echo
echo "-------------------------------------------"
printf 'passed: %d   failed: %d\n' "$pass" "$fail"
[ "$fail" -eq 0 ] || exit 1
