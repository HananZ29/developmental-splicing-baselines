# E14 lock verification, 2026-11-07

`sha256sum -c PREREG_E14.sha256` reports one mismatch: `PREREGISTRATION_E14.md`. All other locked
files match.

**Cause:**
- The lock was written on 2026-11-05 at 01:27:32.
- At 01:31:35 the deviations log was edited once, replacing `*(empty)*` with entries C-E14-1 (G-X
  criterion, Holm-adjusted CI, N_DS-missing rule) and C-E14-2 (MMSplice input conventions). Both
  entries are dated before any scoring output existed.
- This is the only edit to the file after the freeze (session transcript).

**Check:**
- Undoing exactly that edit on the current file gives SHA-256 `22946e1c…d4899b`.
- That is identical to the locked hash.
- So the frozen protocol text is byte-for-byte unchanged. The only post-freeze content is the dated
  deviations log, which is what the log is for.

| Version | SHA-256 of `PREREGISTRATION_E14.md` |
|---|---|
| Frozen (locked), 2026-11-05 01:27 | `22946e1c9989a3907c1f096c08ab5546f93bd3c5b532c0db0e6c33599bd4899b` |
| Current (frozen text + deviations log C-E14-1/2), 2026-11-05 01:31 | `1376e79c4a910f6d69199fc568e437561086609c4b51732575e5e797d91f6d07` |

`PREREG_E14.sha256` is left unchanged, to preserve the original lock. This mismatch must be disclosed
next to the E14 results.
