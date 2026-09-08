# Session S04 — M3 Lambda Sweep 2D Complete
**Date:** 2026-04-19
**Session:** S04

---

## Summary

M3 Lambda Sweep 2D fully implemented, tested, and run to completion.

### Files built

**`layer2/lambda_sweep_2D.py`** — Full sweep driver with:
- SQLite results.db logging (`_ensure_db`, `log_result`, `query_results`, `conjecture_35_verdict`)
- `run_ccf_sweep()` — CCF profiles, all 6 λ values, EXP-L2-R1-CCF-001…006
- `run_boussinesq_sweep()` — 3 Boussinesq profiles, EXP-L2-R1-BOUS-001…003
- `run_adversarial()` — adversarial min (λ→0), EXP-L2-R1-ADV-001
- `run_smoke_test()` — fast N=32 validation for CI
- `print_results_table()`, `--query` CLI flag, `--profile all/ccf/bous/adv`

**`layer2/test_lambda_sweep_2D.py`** — 35 tests across 7 classes:
- TestDatabase (6): file creation, schema, insert, upsert, filter, empty
- TestConjecture35Verdict (3): incomplete/pass/fail logic
- TestCCFSweep (10): all 6 λ pass, A_min > 0, exp IDs, critical λ=0.4703, t_final
- TestBoussinesqSweep (5): 3 profiles, all pass, DB rows
- TestAdversarial (4): λ=0.01 pass, exp_id, DB logging
- TestSmokeTest (3): pass/n_runs/all_pass
- TestM3EndToEnd (3): 9 rows in DB, Conjecture 3.5 PASS verdict, global A_min > 0

---

## Experiment Results

Full sweep: `python3 layer2/lambda_sweep_2D.py --profile all --N 64 --t-end 0.5 --verbose`

| exp_id | λ | verdict | A_min | ω_max |
|---|---|---|---|---|
| EXP-L2-R1-CCF-001 | 1.2000 | ✅ PASS | 3.0801e-05 | 9.89e-01 |
| EXP-L2-R1-CCF-002 | 0.9000 | ✅ PASS | 3.0801e-05 | 9.89e-01 |
| EXP-L2-R1-CCF-003 | 0.6057 | ✅ PASS | 3.0801e-05 | 1.022 |
| EXP-L2-R1-CCF-004 | 0.4703 | ✅ PASS | 3.0801e-05 | 1.037 |
| EXP-L2-R1-CCF-005 | 0.3000 | ✅ PASS | 3.0801e-05 | 1.056 |
| EXP-L2-R1-CCF-006 | 0.1000 | ✅ PASS | 3.0801e-05 | 1.078 |
| EXP-L2-R1-BOUS-001 | 1.9206 | ✅ PASS | 3.0801e-05 | 9.87e-01 |
| EXP-L2-R1-BOUS-002 | 1.3991 | ✅ PASS | 3.0801e-05 | 9.87e-01 |
| EXP-L2-R1-BOUS-003 | 1.1843 | ✅ PASS | 3.0801e-05 | 9.87e-01 |
| EXP-L2-R1-ADV-001 | 0.0100 | ✅ PASS | 3.1740e-05 | 9.999e-01 |

**Conjecture 3.5: PASS — λ_c not found. All λ suppressed.**

---

## Test Results

| Suite | Tests | Status |
|---|---|---|
| layer1/test_props.py | 71 | ✅ PASS |
| layer2/test_self_similar_IC.py | 53 | ✅ PASS |
| layer2/test_M3.py | 64 | ✅ PASS |
| layer2/test_lambda_sweep_2D.py | 35 | ✅ PASS |
| **TOTAL** | **223** | **✅ ALL PASS** |

---

## Next

1. Build `layer2/route2_2D.py` — Route 2: exact Prize NS + auxiliary θ (M2)
2. 4 ε_param runs → EXP-L2-R2-001…004
3. Build `layer2/mu_limit_2D.py` — μ(T)→ν limit (M4)
