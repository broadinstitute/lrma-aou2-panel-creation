# GLIMPSE2 imputation: where the cost actually is

Notes accumulated while fixing the OOMs in `GLIMPSE2FromPreprocessedPLsJoint.wdl`. Recorded so
the rejected options are not re-derived, and so the one change that matters is not lost among
the ones that don't.

Figures are us-central1 `n1-custom` spot (~$0.0332/vCPU-hr, ~$0.00445/GB-hr, spot ≈ 30% of
on-demand), a 250-sample batch, 523 phase shards across 22 chromosomes.

## Summary

**Compute is already at the floor. The gap is preemption churn, and the fix is a backend
config change, not a WDL change.**

| | per batch | per sample |
|---|---|---|
| base compute | ~$27 | ~11¢ |
| × 1.97 full-run-equivalents per success | ~$53 | ~21.5¢ |

11¢ matches Sam's Terra figure, so the base is not where the money is. The 1.97 multiplier is,
and it comes from one thing: preempted shards restart from zero.

## The lever worth pulling

`checkpointingInterval` — `batchConfiguration.batchAttributes`, backend config, not settable
from a WDL.

Cromwell uploads an empty checkpoint placeholder at t=0, then syncs every **600 s**. Median
preemption is at **7 minutes**. So most preempted shards are killed before their first real
sync and resume from nothing — every retry inspected showed `resumed=no`.

Halving the interval should let most preempted shards resume rather than restart, pulling the
multiplier from ~1.97 toward ~1.3 and per-sample cost to roughly **14¢**. That is larger than
every WDL-side change available, combined.

**Action: ask VWB to lower `checkpointingInterval`.** One email.

Related margin note: spot stops paying for itself around 85% preemption. The observed rate is
71%. That is less headroom than it appears, and it is why the zone widening
(`zones` restored to Cromwell's stock four-zone default) matters beyond cost.

## The one WDL lever left, untested

`cpuPlatform` selects the machine family — `GcpBatchMachineConstraints.scala`:

```scala
case Some(CpuPlatformAMDRomeValue)          => N2DCustomMachineType
case Some(CpuPlatformIntelCascadeLakeValue) => N2CustomMachineType
case Some(CpuPlatformIntelIceLakeValue)     => N2CustomMachineType
// unset                                    => N1CustomMachineType
```

Unset means **N1**, the oldest and most expensive family. N2D is roughly 14% cheaper per vCPU
and per GB: about $3/batch, ~$600 across the 50k, for one runtime line. N2/N2D also permit
8.0 GB/cpu against N1's 6.5, which would drop ligate from 6 cpu to 4 — though it rarely helps
phase, where `phase_threads = 4` is usually the binding floor rather than the ratio.

**Not wired in, deliberately.** `setMinCpuPlatform` *restricts* placement to hosts with that
CPU, so it can reduce spot availability. At 71% preemption, losing placement breadth would
cost more than 14% saves. This needs one chromosome and a preemption-rate comparison before
adoption, not a default change.

To test: add `cpuPlatform: "AMD Rome"` to a task's `runtime` block, run one chromosome, and
compare preemption rate and cost against the same chromosome on N1.

## Rejected, with reasons

**Fewer phase threads.** `phase_threads = 2` is ~11% cheaper on paper ($0.031 vs $0.035 per
shard) but stretches the shard from ~35 to ~52 minutes. Against a 7-minute median preemption
the longer task is far likelier to die, and the paper saving evaporates. 4 is near-optimal.
8 threads is worse on both axes ($0.045, and the memory grows with the thread count).

**Larger batches.** Memory is genuinely batch-size independent — `Alpha` sizes on `L × Kpbwt`,
not on sample count — so 1000 samples would amortize the ~6.2 GB panel bin over 4× as many
samples. But localization is only ~6.4 of 35 minutes, capping the gain near 15%, and a 2-hour
task at this preemption rate would rarely finish. This is the mechanism behind Sam's
observation that "1k batches probably too big."

**Pre-splitting the PL VCF per chromosome.** Saves ~17 s of localization per shard, about
**$0.15/batch**. In-region GCS egress is free, so this is a reliability and tidiness win, not
a cost one. Still worth doing eventually — it removes ~1.4 TB of redundant transfer per batch
— but it will not move the bill.

**Dropping the reheader step.** ~2.9 of every 35 minutes, ~$1.5/batch, ~$300 across the 50k.
The stated reason for it is that GLIMPSE2's single-contig header breaks
`bcftools concat --naive` — but nothing concatenates phased shards; ligate consumes them
directly. Possibly removable, needs testing, and the saving is small.

**Lower `--Kpbwt` or fewer MCMC iterations.** Both are real compute reductions and both trade
accuracy. Sam has already described the current settings (Kpbwt 1000, 5 burn-in + 10 main) as
"a bit aggressive" relative to the defaults. Not knobs to turn without him.

## Still open

- **chr11 shards 7 and 8.** Never completed. shard-7 reached attempt-10 with no attempt-11,
  despite the backend source promising a non-preemptible fallback after `maxPreemption`.
  Unexplained, and an unexplained failure mode is worse at 200 batches than a known one.
- **Ligate's 32 GiB.** An empirical cap. Three mechanisms have been proposed and all three
  refuted (see the comment in the WDL). The instrumentation added to that task reports peak
  RSS from an EXIT trap, so the next run of a seam that used to die at 12 GiB settles it.
- **Whether task-level `zones` overrides the backend's `allowedLocations`.** Read
  `allocationPolicy.location` after one chromosome.
