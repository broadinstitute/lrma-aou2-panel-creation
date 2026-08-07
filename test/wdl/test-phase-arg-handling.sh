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
PANEL_WDL="$(dirname "$0")/../../wdl/methods/imputation/GLIMPSE2ChunkAndSplitPanel.wdl"
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
# Per-shard sizing: mem = 2 + ceil(8 * threads * Kpbwt * L / 1e9), cpu = even(max(threads,
# ceil(mem/6.5))). Must clear the observed OOM boundary, stay under the N1 6.5 GB/cpu limit,
# and stay even.
# ---------------------------------------------------------------------------
echo
echo "Memory/CPU sizing"
sizing() {  # kpbwt threads L -> "mem cpu"   (constants read from the WDL, never restated)
    python3 -c '
import math, re, sys
kp, t, L = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
m = re.search(r"Int computed_mem_gb\s*=\s*(\d+(?:\.\d+)?)\s*\+\s*"
              r"ceil\(\(\(\((\d+(?:\.\d+)?)\s*\*\s*phase_threads\)", open(sys.argv[4]).read())
CONST, COEFF = (float(m.group(1)), float(m.group(2))) if m else (2.0, 8.0)
f = re.search(r"if computed_mem_gb > (\d+) then computed_mem_gb else \1", open(sys.argv[4]).read())
FLOOR = float(f.group(1)) if f else 0.0
mem = max(FLOOR, CONST + math.ceil(COEFF * t * kp * L / 1e9))
r = math.ceil(mem / 6.5)
u = r if r > t else t
print(int(mem), u + (u % 2))' "$1" "$2" "$3" "$WDL"
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
chr7_s12 1346888 OOMed_at_16 17
CASES

# The floor is the measured number, not the regression: 517 of 523 shards completed at a flat
# 16 GiB, whereas the fit was built from peaks censored at rc=137 and taken at the 2000-state
# default. Small shards must land on the floor, large ones must still be driven by the model.
FLOOR_I=$(grep -oE 'if computed_mem_gb > [0-9]+ then computed_mem_gb else [0-9]+' "$WDL" | grep -oE '[0-9]+' | head -1)
if [ -n "$FLOOR_I" ]; then
    ok "memory floor present in the WDL (${FLOOR_I} GiB)"
    read -r m_small _ <<<"$(sizing 1000 4 50000)"
    [ "$m_small" -eq "$FLOOR_I" ] && ok "small shard L=50000 floored to ${FLOOR_I} GiB" \
        || bad "small shard floored" "got ${m_small}, expected ${FLOOR_I}"
    read -r m_big _ <<<"$(sizing 1000 4 1346888)"
    [ "$m_big" -gt "$FLOOR_I" ] && ok "dense shard L=1346888 -> ${m_big} GiB, model still governs" \
        || bad "dense shard above floor" "got ${m_big}, floor ${FLOOR_I}"
else
    bad "memory floor present in the WDL" "no floor expression found"
fi

# Real measured peaks from the instrumented chr20 run. The request must exceed each peak
# with margin; at <=75% we would be one bad extrapolation from an OOM. This is the assertion
# that would have caught the original 4.0 coefficient, which reached 90% at L=1M.
while read -r L peak; do
    read -r mem cpu <<<"$(sizing 1000 4 "$L")"
    util=$(python3 -c "print(int(100*$peak/$mem))")
    if [ "$util" -le 75 ]; then
        ok "measured L=$L peak=${peak} GiB fits ${mem} GiB (${util}%)"
    else
        bad "measured L=$L peak=${peak} GiB fits ${mem} GiB" "${util}% of request -- too tight"
    fi
done <<'MEASURED'
331325 7.71
338590 7.55
374957 8.45
396083 9.43
379508 8.43
399495 9.44
400379 9.32
503463 11.20
593895 12.57
965039 21.21
1005104 22.44
MEASURED

# Extrapolated to the densest shards in the genome, using the fitted peak model.
while read -r name L; do
    read -r mem cpu <<<"$(sizing 1000 4 "$L")"
    util=$(python3 -c "print(int(100*(0.373+2.1708e-5*$L)/$mem))")
    if [ "$util" -le 75 ]; then
        ok "$name (L=$L) projects to ${util}% of ${mem} GiB"
    else
        bad "$name (L=$L)" "projects to ${util}% of ${mem} GiB"
    fi
done <<'EXTRAP'
chr2_s17 967196
chr7_s12 1346888
EXTRAP

# Ligate: measured 12.39 GiB peak on chr20's 744,316-site seam (the one that took two rc=137
# kills at 12 GiB). Projected genome max is chr7's 1,023,060-site seam at ~13.5 GiB via
# peak = 9.57 + 3.79e-6 * L_isec. The request must clear both with margin.
LIG_MEM=$(awk '/task GLIMPSE2Ligate/,/runtime \{/' "$WDL" | grep -oE 'mem_gb:\s+[0-9]+' | grep -oE '[0-9]+$')
LIG_CPU=$(awk '/task GLIMPSE2Ligate/,/runtime \{/' "$WDL" | grep -oE 'cpu_cores:\s+[0-9]+' | grep -oE '[0-9]+$')
for probe in "12.39 measured_chr20_seam" "13.45 projected_chr7_seam"; do
    set -- $probe
    ratio=$(python3 -c "print(int(100*$1/$LIG_MEM))")
    if [ "$ratio" -le 75 ]; then
        ok "ligate $2 ${1} GiB fits ${LIG_MEM} GiB (${ratio}%)"
    else
        bad "ligate $2 ${1} GiB fits ${LIG_MEM} GiB" "${ratio}% -- too tight"
    fi
done
if python3 -c "import sys; sys.exit(0 if $LIG_MEM/$LIG_CPU <= 6.5 else 1)"; then
    ok "ligate shape ${LIG_MEM}/${LIG_CPU} within the N1 6.5 GB/cpu limit"
else
    bad "ligate shape ${LIG_MEM}/${LIG_CPU} within the N1 limit" "exceeds 6.5"
fi

# Pop: measured worst peak 4.68 GiB across eleven chr20 shards. It scatters ~523 times per
# batch, so this is the task where over-provisioning actually costs something.
POP_MEM=$(awk '/task PopAndMarginalizeCollisions/,/runtime \{/' "$WDL" | grep -oE 'mem_gb:\s+[0-9]+' | grep -oE '[0-9]+$')
POP_CPU=$(awk '/task PopAndMarginalizeCollisions/,/runtime \{/' "$WDL" | grep -oE 'cpu_cores:\s+[0-9]+' | grep -oE '[0-9]+$')
pop_ratio=$(python3 -c "print(int(100*4.68/$POP_MEM))")
if [ "$pop_ratio" -le 75 ]; then
    ok "pop measured 4.68 GiB fits ${POP_MEM} GiB (${pop_ratio}%)"
else
    bad "pop measured 4.68 GiB fits ${POP_MEM} GiB" "${pop_ratio}% -- too tight"
fi
python3 -c "import sys; sys.exit(0 if $POP_MEM/$POP_CPU <= 6.5 else 1)" \
    && ok "pop shape ${POP_MEM}/${POP_CPU} within the N1 6.5 GB/cpu limit" \
    || bad "pop shape ${POP_MEM}/${POP_CPU}" "exceeds 6.5 GB/cpu"

# Threading must not exceed the cpu the task pays for, or the threads contend.
LIG_TH=$(awk '/task GLIMPSE2Ligate/,/^}/' "$WDL" | grep -oE 'Int ligate_threads = [0-9]+' | grep -oE '[0-9]+$' | head -1)
if [ "$LIG_TH" -le "$LIG_CPU" ]; then
    ok "ligate threads ($LIG_TH) <= cpu ($LIG_CPU) at the default"
else
    bad "ligate threads ($LIG_TH) <= cpu ($LIG_CPU)" "oversubscribed"
fi
# ...and under a memory override, which the default check cannot see. The cpu floor must be
# the thread count, not a hardcoded 2, or overriding mem_gb down oversubscribes the threads.
if awk '/task GLIMPSE2Ligate/,/^}/' "$WDL" | grep -q 'eff_ratio_cpu > ligate_threads'; then
    ok "ligate cpu floors on ligate_threads, so overrides cannot oversubscribe"
else
    bad "ligate cpu floors on ligate_threads" "floor is hardcoded; mem_gb override can oversubscribe"
fi

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
# The fixed term comes from the WDL: restating it here made both assertions fail the moment
# it moved from 2 to 4, which is a test tracking a literal rather than the property.
WDL_CONST_I=$(grep -oE 'Int computed_mem_gb[[:space:]]*=[[:space:]]*[0-9]+' "$WDL" | grep -oE '[0-9]+$')
read -r m1 _ <<<"$(sizing 1000 4 1000000)"
read -r m2 _ <<<"$(sizing 2000 4 1000000)"
read -r m3 _ <<<"$(sizing 1000 8 1000000)"
[ $((m2 - WDL_CONST_I)) -eq $(( (m1 - WDL_CONST_I) * 2 )) ] && ok "Kpbwt 1000->2000 doubles the matrix term" \
    || bad "Kpbwt scaling" "$m1 -> $m2"
[ $((m3 - WDL_CONST_I)) -eq $(( (m1 - WDL_CONST_I) * 2 )) ] && ok "threads 4->8 doubles the matrix term" \
    || bad "thread scaling" "$m1 -> $m3"

# Per-shard sizing must not cost more than the flat worst-case bound it replaces.
python3 - "$WDL" <<'COST'
import math, re, sys
CPU, MEM, SPOT = 0.0332, 0.00445, 0.30
# Read the sizing constants out of the WDL rather than restating them: this block validated a
# formula the WDL had already moved off once, and printed a stale cost table while every other
# assertion around it read the live value.
m = re.search(r"Int computed_mem_gb\s*=\s*(\d+(?:\.\d+)?)\s*\+\s*"
              r"ceil\(\(\(\((\d+(?:\.\d+)?)\s*\*\s*phase_threads\)", open(sys.argv[1]).read())
CONST, COEFF = (float(m.group(1)), float(m.group(2))) if m else (2.0, 8.0)
def size(L, t=4, kp=1000):
    m = CONST + math.ceil(COEFF * t * kp * L / 1e9)
    r = math.ceil(m / 6.5); u = r if r > t else t
    return m, u + (u % 2)
def cost(n, cpu, gb, mins): return n * (cpu*CPU + gb*MEM) * (mins/60) * SPOT
Ls = [400000]*491 + [657784,765770,965039,1005104,1346888] + [600000]*27
dyn = sum(cost(1, c, m, 35) for m, c in (size(L) for L in Ls))
flat = cost(len(Ls), 8, 40, 35)   # the worst-case bound this replaced
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
    if grep -qE 'ceil\(\(\(\([0-9]+\.[0-9]+ \* phase_threads\) \* phase_kpbwt\) \* n_variants\)' "$WDL"; then
        ok "memory expression promotes to Float before multiplying"
    else
        bad "memory expression promotes to Float before multiplying" "expression reordered or coefficient changed"
    fi

    # Instrumentation must be trap-based and on stderr, or it vanishes exactly when needed.
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

    # The harvester must derive the sizing constants from the WDL, not hardcode them: they
    # drifted once already (WDL refitted 4.0 -> 8.0 while the harvester still asserted 4.0,
    # which would have reported a bound breach on every shard).
    HARVEST_MODEL=$(python3 -c "
import importlib.util as u, sys
sp=u.spec_from_file_location('h','$(dirname "$0")/../../scripts/harvest-resources.py')
m=u.module_from_spec(sp); sp.loader.exec_module(m)
print('%s %s' % m.wdl_model())" 2>/dev/null)
    WDL_CONST=$(grep -oE 'Int computed_mem_gb\s*=\s*[0-9]+' "$WDL" | grep -oE '[0-9]+$')
    WDL_COEFF=$(grep -oE 'ceil\(\(\(\([0-9]+\.[0-9]+ \* phase_threads' "$WDL" | grep -oE '[0-9]+\.[0-9]+')
    if [ "$HARVEST_MODEL" = "$WDL_CONST.0 $WDL_COEFF" ]; then
        ok "harvester reads the sizing model from the WDL ($HARVEST_MODEL)"
    else
        bad "harvester reads the sizing model from the WDL" "harvester=[$HARVEST_MODEL] wdl=[$WDL_CONST.0 $WDL_COEFF]"
    fi

    # A workflow input passed straight through to a task overrides that task's default, so the
    # two must agree or the task default is dead and the comment beside it lies. A revert
    # changed the ligate task default to 2 and left the workflow default at 4, which shipped a
    # 4-thread ligate sized against a 2-thread measurement. Checked for EVERY passed-through
    # input, not just this one, because the failure is structural rather than specific.
    if python3 - "$WDL" <<'PASSTHRU'
import re, sys
s = open(sys.argv[1]).read()
wf = s[:s.index("\nstruct RuntimeAttr")]
decls = wf[wf.index("input {"):wf.index("\n    Map[String, String] genetic_maps_dict")]
wf_in = dict(re.findall(r"^\s*(?:Int|Float|String|Boolean)\s+(\w+)\s*=\s*([^\n]+?)\s*$", decls, re.M))
passed = set(re.findall(r"^\s*(\w+)\s*=\s*\1\s*,?\s*$", wf, re.M))
bad = []
for name in sorted(passed & set(wf_in)):
    for m in re.finditer(r"task (\w+) \{", s):
        seg = s[m.end():]
        nx = re.search(r"\ntask ", seg)
        seg = seg[:nx.start()] if nx else seg
        t = re.search(r"^\s*(?:Int|Float|String|Boolean)\s+%s\s*=\s*([^\n]+?)\s*$" % name, seg, re.M)
        if t and t.group(1).strip() != wf_in[name].strip():
            bad.append("%s: workflow=%s %s=%s" % (name, wf_in[name], m.group(1), t.group(1)))
if bad:
    print("      " + "; ".join(bad))
sys.exit(1 if bad else 0)
PASSTHRU
    then ok "workflow and task defaults agree for every passed-through input"
    else bad "workflow and task defaults agree for every passed-through input" "the workflow value wins"; fi

    # Structural, so a task added later cannot silently miss either fix. Every task with a
    # runtime block must re-derive cpu from effective memory (or a partial runtime_attr_override
    # breaches the N1 ratio) and must install the EXIT-trap instrumentation (or it reports
    # nothing from the OOMs that matter).
    if python3 - "$WDL" <<'AUDIT'
import re, sys
s = open(sys.argv[1]).read()
bad = []
for m in re.finditer(r"task (\w+) \{", s):
    seg = s[m.end():]
    nx = re.search(r"\ntask ", seg)
    seg = seg[:nx.start()] if nx else seg
    if "runtime {" not in seg:
        continue
    missing = [n for n, ok in (("eff_cpu", "Int eff_cpu" in seg),
                               ("instrumentation", "trap _instr_report" in seg),
                               # boot_disk_gb must stay WIRED, not a silent no-op: an explicit
                               # 0 is provably identical to omitting the attribute (Cromwell's
                               # BootDiskSizeValidation gives 0 + 30 either way), so there is
                               # no reason for the override to be dead.
                               ("boot_disk default", "boot_disk_gb:       0," in seg),
                               ("bootDiskSizeGb wired", "bootDiskSizeGb:" in seg)) if not ok]
    if missing:
        bad.append("%s (%s)" % (m.group(1), ", ".join(missing)))
if bad:
    print("      " + "; ".join(bad))
sys.exit(1 if bad else 0)
AUDIT
    then ok "every task: eff_cpu, instrumentation, boot_disk default and wiring"
    else bad "every task: eff_cpu, instrumentation, boot_disk default and wiring" "see above"; fi

    # A resource edit meant for one task landing on another is invisible in review: both
    # default_attr blocks look alike, and womtool accepts either. This branch downsized
    # GLIMPSE2Chunk from 4/8 to 2/4 while leaving the task it meant to change untouched.
    # Every task this branch did not introduce must keep the resources Sam gave it.
    if python3 - <<'RESDRIFT'
import re, subprocess, sys
P = "wdl/methods/imputation/GLIMPSE2ChunkAndSplitPanel.wdl"
def res(text):
    out, n = {}, None
    for line in text.splitlines():
        m = re.match(r"task (\w+)", line)
        if m: n = m.group(1)
        m = re.search(r"(cpu_cores|mem_gb):\s*(\d+)", line)
        if m and n: out.setdefault(n, {})[m.group(1)] = m.group(2)
    return out
# The merge-base, not sl_aou2_v1 itself: that branch moves, and comparing against its tip
# would report Sam's later changes as our drift.
mb = subprocess.run(["git", "merge-base", "HEAD", "origin/sl_aou2_v1"],
                    capture_output=True, text=True)
if mb.returncode != 0:
    sys.exit(0)   # no base ref (shallow clone); nothing to compare against
base = subprocess.run(["git", "show", mb.stdout.strip() + ":" + P],
                      capture_output=True, text=True)
if base.returncode != 0:
    sys.exit(0)
old, new = res(base.stdout), res(open(P).read())
drift = ["%s %s -> %s" % (k, old[k], new.get(k)) for k in old if old[k] != new.get(k)]
if drift:
    print("      " + "; ".join(drift))
sys.exit(1 if drift else 0)
RESDRIFT
    then ok "pre-existing panel tasks keep their original resources"
    else bad "pre-existing panel tasks keep their original resources" "resources changed"; fi

    # A bare `wait` in a command block that starts the instrumentation sampler blocks
    # forever: `wait` with no arguments waits for every background job, and the sampler loops
    # until the EXIT trap kills it. That deadlocked a non-preemptible task in production.
    if python3 - "$WDL" "$PANEL_WDL" <<'BAREWAIT'
import re, sys
bad = []
for path in sys.argv[1:]:
    s = open(path).read()
    for m in re.finditer(r"task (\w+) \{", s):
        seg = s[m.end():]
        nx = re.search(r"\ntask ", seg)
        seg = seg[:nx.start()] if nx else seg
        if "_INSTR_SAMPLER=$!" in seg and re.search(r"^\s*wait\s*$", seg, re.M):
            bad.append("%s in %s" % (m.group(1), path.split("/")[-1]))
if bad:
    print("      " + "; ".join(bad))
sys.exit(1 if bad else 0)
BAREWAIT
    then ok "no bare wait alongside the instrumentation sampler"
    else bad "no bare wait alongside the instrumentation sampler" "deadlock"; fi

    # The count comes from the panel, not a per-run task: it is a property of the panel and
    # recomputing it once per chromosome per batch produced the same numbers 4,400 times.
    if grep -q 'Array\[Int\] n_variants = chunked_panel\[chromosome\].n_variants' "$WDL"; then
        ok "consumer reads n_variants from the panel"
    else
        bad "consumer reads n_variants from the panel" "field not read"
    fi
    if grep -q 'n_variants: CountPanelVariantsPerShard.n_variants' "$PANEL_WDL"; then
        ok "producer emits n_variants into the panel JSON"
    else
        bad "producer emits n_variants into the panel JSON" "field not emitted"
    fi
    # Both structs must agree or the JSON will not round-trip.
    a=$(awk '/^struct ChunkedPanelChromosome/,/^}/' "$WDL" | grep -oE '^\s+(String|Array\[[A-Za-z]+\])\s+\w+' | tr -s ' ')
    b=$(awk '/^struct ChunkedPanelChromosome/,/^}/' "$PANEL_WDL" | grep -oE '^\s+(String|Array\[[A-Za-z]+\])\s+\w+' | tr -s ' ')
    if [ "$a" = "$b" ]; then
        ok "ChunkedPanelChromosome identical in producer and consumer"
    else
        bad "ChunkedPanelChromosome identical in producer and consumer" "structs differ"
    fi

    # The disk floor must stay a parameter. It is the largest remaining SSD-quota saving and
    # is gated on a measurement; hardcoding it again would put that experiment behind a WDL
    # edit rather than an input override.
    if grep -q 'computed_disk_gb > phase_disk_floor_gb' "$WDL"; then
        ok "phase disk floor is parameterised, not hardcoded"
    else
        bad "phase disk floor is parameterised" "floor inlined again"
    fi

    # The counting task serialises the chromosome; it must not be preemptible.
    if awk '/task CountPanelVariantsPerShard/,/^}/' "$PANEL_WDL" | grep -qE 'preemptible_tries:\s*0,'; then
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
    # Peak derived from the WDL's own constants: a literal stops breaching the bound as soon
    # as the coefficient is refitted upward, and the test then asserts nothing.
    BREACH_PEAK=$(python3 -c "
import re
m = re.search(r'Int computed_mem_gb\s*=\s*(\d+(?:\.\d+)?)\s*\+\s*'
              r'ceil\(\(\(\((\d+(?:\.\d+)?)\s*\*\s*phase_threads\)', open('$WDL').read())
c, k = (float(m.group(1)), float(m.group(2))) if m else (2.0, 8.0)
print('%.2f' % (c + 1.25 * k * 4 * 1000 * 1346888 / 1073741824))")
    OVER="[RESOURCE] task=GLIMPSE2Phase n_variants=1346888 requested_mem_gib=30 requested_cpu=6 threads=4 kpbwt=1000 peak_rss_gib=$BREACH_PEAK peak_rss_source=cgroup-v2 wall_s=1"
    echo "$OVER" | "$HARVEST" 2>/dev/null | grep -q "OVER REQUEST" \
        && ok "harvester flags an under-sized request" \
        || bad "harvester flags an under-sized request" "no OVER REQUEST warning"
    echo "$OVER" | "$HARVEST" 2>/dev/null | grep -q "NOT an upper bound" \
        && ok "harvester flags a breached memory bound" \
        || bad "harvester flags a breached memory bound" "no bound warning"

    # over-provisioning must be flagged too, or the data never drives sizing down
    UNDER='[RESOURCE] task=GLIMPSE2Ligate n_shards=24 requested_mem_gib=24 requested_cpu=4 peak_rss_gib=1.86 peak_rss_source=cgroup-v2 wall_s=1'
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
