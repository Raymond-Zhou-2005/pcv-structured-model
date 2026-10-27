# Numerical portability audit

This audit supplements the frozen manuscript-linked `v1.0.0` release. It does
not alter the scientific solver, numerical acceptance tolerances, inputs,
manuscript results or frozen publication figures.

## Recorded runs

- The fresh reference run on Windows x86-64 / CPython 3.12.14 passed 811/811
  checks. Its source/table comparisons are in `../publication_validation.json`.
- [Initial Linux run](https://github.com/Raymond-Zhou-2005/pcv-structured-model/actions/runs/40555237696):
  all model analyses and figures completed, followed by 27 failed strict
  reference-column comparisons (`rtol=1e-12`, `atol=1e-7`). The largest reported
  full-population monetary residual was approximately 0.000380 scenario USD.
  `initial_linux_failures.json` preserves the failure list. The original runner's
  detailed output files were not archived, so individual failing rows and an
  exact root-cause attribution cannot be recovered from that log.
- [Diagnostic Linux rerun](https://github.com/Raymond-Zhou-2005/pcv-structured-model/actions/runs/40557076944):
  Ubuntu / CPython 3.12.15 passed the unchanged 811/811 checks. Only failure
  reporting and artifact retention had changed. All 144 input/economic-parameter
  records matched Windows semantically, and all state arrays agreed within the
  original tolerances (maximum absolute difference 2.91e-11). See
  `linux_pass_verification.json` and `cross_platform_summary.json`.

A passing rerun does not erase the first failure or certify every hardware and
library combination. Windows is the fully tested reference environment. The
Linux workflow is available by manual dispatch for further investigation;
failures remain failures and are not converted into passing checks.

## Independent numerical probes

`fixed_input_month_oracle.py` compares the retained RK45 monthly solve with
the frozen-coefficient exponential/exposure formula, under identical inputs,
for no vaccination and the three full111 scenarios at months 0, 24, 60 and 119.
Default and tighter RK45 tolerances are retained in the CSV. Maximum default
single-month death discrepancy was 1.28e-13, corresponding to approximately
5.26e-9 USD in that month's discounted death valuation. Tighter tolerances did
not improve this probe. It is not a bound on the full nonlinear trajectory.

The same script audits probability/hazard conversions against 80-digit Decimal
calculations conditional on the same binary64 input. The inherited `1-exp(-x)`
and `-log(1-p)` forms lose relative precision at small risks; `expm1`/`log1p`
reduce that loss. The current release retains the inherited forms. These tests
identify a numerical sensitivity, without proving the cause of the earlier
Linux failure or validating epidemiological inputs. Frozen coefficients and
all probe rows are included.

From the repository root, the optional probe can be rerun separately:

```sh
python checks/portability/fixed_input_month_oracle.py --repo . --output checks/logs/month_oracle
```

Future portability work should compare per-month primitive flows and stocks,
trace their fixed economic weights, and propagate error budgets through
contrasts and ratios. Simply increasing one global tolerance would conceal the
conditioning of small difference estimates. None of that future work is
represented as completed here.
