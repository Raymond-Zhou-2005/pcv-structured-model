# Data definitions and estimands

## CFR translation diagnostic in the local revision

Case and fatal-event rewards share the saved monthly alive exposure. The
inherited probability transformation uses lambda = -log(1-p) and
mu = -log(1-p*f), so their ratio is an effective fatal-event coefficient,
mu/lambda. For pneumonia, f is the severe fraction times the severe-case CFR;
the effective CFR among severe episodes therefore uses severity*lambda as its
denominator. Meningitis sequela sources retain (1-CFR)*episode rewards as a
separate coefficient convention.

`checks/revision_20261103/cfr_translation_audit.py` reconstructs these quantities
from the 144 saved trajectories and reports their differences without changing
the population trajectories, sequela stocks, QALYs or costs. It is a fixed-exposure
reward diagnostic; changing the underlying death hazards would require a full
alternative simulation.

## Current open-population application

`research/inputs/Final/outputs/runs/<scenario>/` contains a complete model run.
There are 144 configurations: 37 core cases, 27 schedule cases and 80
perturbation/comparator configurations. The 60 vaccinated perturbations appear
in `sensitivity.csv`; their 20 matched no-vaccination comparators remain in run
directories.

| File | Definition |
|---|---|
| `parameters.json` | Exact epidemiological/programme constants for this run |
| `economic_parameters.json` | Exact valuation constants and disclosed proxy conversions |
| `states.npz` | `alive_after_births_and_doses`: 120 recorded months × 5 alive dose states × 60 attained-age bins; entries are bin masses, not point densities |
| `monthly.csv` | Expected endpoint rewards and alive/sequela quantities at the declared monthly timing |
| `doses.csv` | Actual modelled first, second and booster administrations |
| `herd.csv` | Weighted dose proxy C, normalised r, target and realised endpoint multipliers H |
| `balance.csv` | Population/dose conservation checks |
| `yearly_unrounded.csv` | Annual aggregation of unrounded event outputs |
| `economic_details.csv` | Calendar-discounted event, incident-sequela and stock valuation ledger |

Check the actual CSV headers as the canonical field list. Historical files may
contain additional bookkeeping columns. The five alive programme states are
non-programme, enrolled/no dose, one dose, two doses and booster. These are
different from the constructed CTMC verification states, whose state 4 is
absorbing.

The full factorial uses binary D/W/G: delayed/incomplete delivery, attained-age
waning and additional state-dependent herd build-up. `000` retains annual
target-derived indirect protection. `G=0` does not disable herd; matched
`herd_off` runs set H=1. Mechanism changes and paired indirect contributions
therefore have different estimands.

Main gains are differences from the matched no-vaccination run, per 100,000
gross recorded-window births. This denominator is neither one fully followed
birth cohort nor the count of all exposed children. Cases and deaths are
expected event rewards; reported deaths are the pneumococcal endpoint rewards,
while neonatal/background mortality also changes exposure. QALYs and costs use
the declared mixed horizon and scenario monetary convention. A negative ICER
has a dominance label; negative ratios are not a product-ranking scale.

## Numerical and stochastic checks

`research/outputs/refinement.csv` uses k=1,2,4,8,16 substeps per month. The finest
computed grid is a finite reference, not an exact solution; a zero difference
at k=16 is definitional. Dose opportunities are applied once at the scheduled
integer age, not once per substep.

`research/outputs/aggregation_states.npz` contains the two constructed
50,000-person populations with identical age and dose marginals but different
age–dose associations. The matching CSV gives one frozen-coefficient month of
expected disease rewards and the initial feedback indices.

`revision/outputs/ctmc_variability.csv` contains 20 time–state summaries.
`ctmc_absorbing_event.csv` contains four absorbing-event probabilities and exact
binomial count prediction intervals for N=100,000. `ctmc_absorbing_replicates.csv`
contains the 4,000 independently drawn terminal absorbing counts. The seed rule
is 20261004 + int(month × 100), with 1,000 draws at each time. The four time
samples are independent; they are not linked trajectories or parameter draws.

`research/outputs/cohort/` retains the separate seed-matched closed-cohort
microsimulation audit. Its population, protection convention and economics are
not interchangeable with the current open-calendar factorial. Historical
200-draw PSA claims lack recovered original draws/seed and are not reproduced
as current calendar uncertainty evidence.

## Sources and evidence status

`research/sources/parameter_sources.csv` provides the exact input, source role,
transformation and unresolved transfer. `input_comparison.csv` is an internal
arithmetic consistency comparison, not an independent primary-data validation.
`productivity_age_weights.csv` reconstructs all 46 stored contributions.
Bibliographic identity metadata do not verify every computational number.
