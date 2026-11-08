# E14 execution note, 2026-11-07, about 00:45 (schedule change only)

**Reason:** the user needed E14 finished by 2026-11-08. On the sequential schedule in `run_steps3_4.sh`,
step 4 would only start after the N_DS subsample finished (about 02:30 on 2026-11-08), so E14 would
have ended around 2026-11-09.

**What changed:** step 4 (Pangolin and MMSplice on T-E13, T-Vex and T-MFASS) was started **while** the
step-3 N_DS subsample was still running, using `run_step4_parallel.sh`.
- Same scripts, inputs and output paths as step 4 of `run_steps3_4.sh`.
- The three sets run concurrently.
- Thread caps: Pangolin `OMP_NUM_THREADS=2`; MMSplice `OMP_NUM_THREADS=1`.

**What did not change:**
- The pre-registration, the scripts, the inputs and the seeds.
- The N_DS subsample: it was not stopped or restarted.
- "No outcome is examined until all of steps 2–4 are complete": still holds.
- E14c is still run once, after the code-path test.

**Duplicate step 4:** `guard_old_step4.sh` waits until `run_steps3_4.sh` logs "step4 start", which only
happens after the N_DS subsample has exited. It then stops that script's bash process and its own child
processes, so step 4 is not computed twice. Nothing else is stopped.

**Possible numerical effect:** limiting CPU thread counts can change floating-point summation order, so
results may differ at about the 1e-6 level from a single-job run. This has no bearing on any
pre-registered decision.

## Addendum, 2026-11-07 06:45

- **Standby gap:** the laptop was in Modern Standby from 01:06 to 06:39 (idle timeout, on battery). The
  processes were suspended, not killed, and resumed afterwards. No outputs were affected.
- **Keep-awake:** `keep_awake_until_e14_done.ps1` holds `ES_SYSTEM_REQUIRED` while any E14 scoring
  process is alive. No settings are changed.
- **Priority:** the Pangolin and MMSplice processes were set to BelowNormal so that the N_DS subsample,
  the slowest remaining step, gets CPU first. This is scheduling only; results are unaffected.

## Addendum, 2026-11-08 04:05: completion of scoring
- Parallel step 4 finished 02:33:26 (all exit 0). Row counts match between Pangolin and MMSplice: e13test 1,573; vexseq 1,960; mfass 3,000 (plus header).
- DeltaSplice native subsample finished 200/200 (`ds_native_fresh.tsv`, 1,001 lines); no Traceback in `log_ds_fresh.txt`.
- The original chain reached "step4 start" at 04:04:49. The guard killed its Pangolin e13test processes at 04:04:56 (7 s). No score or track file was modified: `pangolin_e13test/scores.tsv` and `tracks.npz` keep their 2026-11-07 21:03 timestamps. The only side effect is that the duplicate run truncated and overwrote `log_pangolin_e13test.txt`, so the original e13test Pangolin log is not preserved. Its completion is recorded in the parallel-run task output ("pangolin e13test done exit=0").
- Keep-awake script exited once no scoring process remained.
