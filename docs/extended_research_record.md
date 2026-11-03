# Extended research record: mathematical methods, application outputs and evidence

> Repository research record, not a separate journal supplement. The original S1–S15 identifiers and displayed table precision are retained for traceability. Core model definitions and arguments are given in the main manuscript; this record preserves the full extended evidence, historical comparisons and qualifications. Full-precision numerical sources remain in `research/` and `revision/`.
>
> Reference keys link to `reference_key_map.csv`, which records bibliographic identity, not verification of every transferred parameter. Parameter-specific evidence status remains in `research/sources/parameter_sources.csv`. Historical cost years, estimands, unrecovered PSA draws and unverified proxy inputs retain their original limitations.

Abbreviations: PCV, pneumococcal conjugate vaccine; AOM, acute otitis media; IPD, invasive pneumococcal disease; QALY, quality-adjusted life year; VE, vaccine effectiveness; CFR, case-fatality ratio.

## S1. Characteristic formula and local nonlinear well-posedness

For prescribed feedback $`h(t)`$, write $`U_h(t,s;a_0)`$ for the fundamental solution along $`a=a_0+t-s`$, including the ordered interface maps. Initial-population characteristics give $`A(t,a)=U_h(t,0;a-t)A_0(a-t)`$ for $`a>t`$; birth characteristics give $`A(t,a)=U_h(t,t-a;0)b(t-a)`$ for $`a<t`$. At interface ages the appropriate one-sided trace and map are used. These formulas also define measure-valued cohort propagation.

For a bounded weight $`w`$ and state weights $`\rho`$, let $`D(A)=\sum_g\int wA_g`$ and $`N_\rho(A)=\sum_g\rho_g\int wA_g`$. On a bounded positive set where $`D(A),D(B)\ge\delta`$, both numerator and denominator are bounded linear functionals. Splitting the ratio difference gives
``` math
\left|\frac{N_\rho(A)}{D(A)}-\frac{N_\rho(B)}{D(B)}\right|
 \le\frac{|N_\rho(A)-N_\rho(B)|}{\delta}
    +\frac{|N_\rho(B)|\,|D(A)-D(B)|}{\delta^2}
 \le L_C\|A-B\|_1.
```
Clipping is Lipschitz. Thus the nonlocal multiplier is locally Lipschitz on this set. For bounded locally Lipschitz rates, the characteristic integral operator obeys a difference estimate of the form $`\|\mathcal T A-\mathcal T B\|_{C([0,T];L^1)}\le TL_F\|A-B\|_{C([0,T];L^1)}`$ on a suitably bounded invariant neighbourhood. For sufficiently small $`T`$ it is a contraction. Fixed birth data and bounded linear interface maps preserve the estimate; successive characteristic strips extend the solution while the denominator bound and coefficient bounds persist. The proof establishes local uniqueness under these assumptions. It does not guarantee global continuation across vanishing weighted population, arbitrary discontinuous feedback or newly specified policies.

### Revision: explicit local neighbourhood and density counterexample

The current manuscript makes the local fixed-point construction explicit.
Assume a bounded strongly continuous positive transport evolution with the fixed
dose interfaces, and initial weighted mass D0 > 0. The reaction-free evolution
z(t) has weighted mass at least 3D0/4 for sufficiently short time. A closed
non-negative path tube of radius R < D0/(4 ||w|| infinity) around z then has
weighted mass at least D0/2. Bounded reaction coefficients give a self-map bound
C0 T, and variation of constants gives the contraction bound C1 T exp(C1 T).
Taking time small enough for both bounds proves the stated local result.
The manuscript gives the ratio inequality and continuation conditions explicitly.

The aggregation proposition now uses normalised indicator densities on two
disjoint intervals of positive length. Swapping their assignments between two
states preserves both marginals. A nonzero interval-averaged cross-difference
gives a nonzero reward difference. This places the diagnostic construction in
the same L1 space as the main density formulation. The saved numerical figure
continues to represent monthly bin masses.

## S2. Why a finite stochastic interacting model is not exactly its mean equation

For $`X=N^{-1}\operatorname{Binomial}(N,p)`$, $`\mathbb E[X^2]=p^2+p(1-p)/N`$. With $`N=100`$ and $`p=0.5`$, this is 0.2525, whereas $`(\mathbb E[X])^2=0.25`$. Consequently, replacing a nonlinear random feedback by the feedback at the mean can discard covariance terms. The continuous-time chain verification in the main text uses externally specified rates; the PCV herd functional is deterministic population feedback. No exact finite-particle expectation equivalence is asserted between them.

The five-state verification generator is
``` math
Q=\begin{pmatrix}
-.21&.20&0&0&.01\\
0&-.162&.15&0&.012\\
0&0&-.098&.08&.018\\
0&0&0&-.022&.022\\
0&0&0&0&0
\end{pmatrix}.
```
The initial probability is $`(1,0,0,0,0)^\top`$. Forward integration uses RK45 with relative tolerance $`10^{-11}`$ and absolute tolerance $`10^{-13}`$. Terminal multinomial samples use 1,000 independent replicates of 100,000 individuals at each observation time; the seed is 20261004 plus the integer number of hundredths of a month.

#### Terminal variability and model-implied event probabilities.

For the fixed verification generator, the count vector of $`N=100{,}000`$ independent individuals is multinomial. Its marginal state-$`i`$ count is binomial, with mean $`Np_i(t)`$, variance $`Np_i(t)\{1-p_i(t)\}`$ and cross-state covariance $`-Np_i(t)p_j(t)`$. Since state 4 is absorbing and initially empty, $`p_4(t)`$ is the probability of absorption by $`t`$ under this constructed process. The forward ODE and matrix exponential give the same terminal probabilities. The Markov interpretation identifies the probability law around the expected counts; this is distinct from a deterministic mean-count report.

The original seeds regenerate 1,000 terminal samples at each observation time. Source data retain all 4,000 absorbing counts and all 20 time–state summaries. Count standard deviations are compared with the binomial expression, and the tabulated central 95% prediction intervals use the exact binomial quantiles. These are intervals for random terminal counts conditional on the stated generator and population size, not confidence intervals for the mean or uncertainty in the generator. Sampling times are independent; no linked trajectories or path-dependent endpoint risks are inferred. The verification states and rates are constructed and are not PCV clinical estimates.

<div id="tab:ctmc_variability">

| Month | Probability | Mean count | Exact SD |  MC SD | 95% count interval |
|------:|------------:|-----------:|---------:|-------:|-------------------:|
|  0.25 |    0.002510 |     250.96 |    15.82 |  16.17 |       \[220, 282\] |
|     1 |    0.010163 |   1,016.30 |    31.72 |  31.47 |     \[955, 1,079\] |
|     3 |    0.031623 |   3,162.28 |    55.34 |  53.69 |   \[3,054, 3,271\] |
|    12 |    0.146056 |  14,605.58 |   111.68 | 114.06 | \[14,387, 14,825\] |

Terminal absorbing-event distribution in the constructed verification chain. Counts concern 100,000 independent individuals. SD denotes standard deviation; the Monte Carlo (MC) estimate uses 1,000 seeded terminal samples at each observation time. The central 95% binomial interval is a prediction interval for a random count with fixed inputs.

</div>

## S3. Transport discretisation, rewards and accounting

Cell arrays store masses, not point densities. Births enter at the left boundary; the reaction is exact for the frozen coefficient set, then survivor masses shift by one cell at Courant number one. Conditional dose maps operate once at each integer-age visit. Refined birth counts are $`B_m/k`$ at each substep, and attained ages determine piecewise-constant risk and waning profiles. This preserves programme opportunities while specifying a within-month interpolation. The finest grid is a numerical reference only. The manufactured test compares reconstructed cell masses with exact integrals over those cells, including the moving front.

Sequela stocks are aggregated over ages. Each reaction step scales existing stocks by overall alive survival and adds fractions of meningitis survivor events. This is an aggregate recorded-window accounting approximation. Annual economics use summed unrounded events and year-end stocks. One-time hearing, epilepsy, hemiplegia and implant costs use actual new sequela survivor flows summed across the year; net annual stock increments are retained only as historical diagnostics. Chronic losses and recurring costs stop at the recorded-window end, whereas deaths use remaining-life valuation. This mixed horizon is not a complete lifetime economic analysis. Acute QALYs are events multiplied by duration/365 and disutility. Death QALYs use remaining life from assumed mean death age 2.5, valued relative to the event and then discounted at calendar-year midpoint. Actual dose counts are costed at time $`y-1+0.35`$; disease and other annual quantities at $`y-0.5`$. The input price alignment uses a retained 3% proxy; PCV13 price 4.22 retains its later source year. The reported currency is scenario USD, with exact conversions and unresolved source status recorded below.

With a remaining-life weight $`W(a)`$ already discounted from birth, translating to an event-relative quantity requires multiplying by $`(1+d)^a`$ before applying calendar discount $`(1+d)^{-t}`$. Omitting that conversion treats cohort age as calendar time. Earlier independent timing-economic outputs affected by this error are not included in current estimates. This correction is an accounting requirement and not evidence for a biological mechanism.

## S4. Research coverage and evidence hierarchy

<div id="tab:coverage">

| Research strand | Treatment | Location |
|:---|:---|:---|
| Research strand | Treatment | Location |
| Piecewise-constant cohort propagation | Exact exponential reaction; matched ODE comparison | Main Methods; S2–S3 |
| Accounting, parameter and structural tiers | Separate output conventions, transferred inputs and model mechanisms | Main evidence status; S5 |
| Monthly stochastic / deterministic cohort | Retained 1,000-replicate audit rechecked; historical comparisons | Main verification; S6 |
| AOM denominator reconciliation | All-cause versus pneumococcal-attributed events | Main AOM; S6 and S10 |
| Age–time cohort limit | Characteristic derivation; earlier simplified limit retained as historical | Main propositions; S1 |
| Open calendar population | Warm-up, continuing births and age-out | Main transport; Fig 3 |
| Endpoint-specific VE and SII bridge | Three profiles; proxy and serotype-weight limitations | Main application; S5 |
| Programme dose states and waning | Five states, finite windows, attained-age approximation | Main application; S5 |
| Indirect protection | C, r, target and actual H; paired off counterfactual | Main functional; S7 |
| Health-economic refinement | Unrounded outputs, calendar discount and actual dose costs | Main outputs; S3, S7 |
| Factorial and interaction analysis | Eight configurations per product; Shapley closure | Main Fig 4; S7 |
| Schedule and completion recovery | Nine rules per product with common W/G | Main Fig 5; S8 |
| Deterministic perturbations | 20 settings and matched NV; perspectives | Main Fig 6; S9 |
| Late-period outcomes | Year 10 endpoint summaries; no fitted equilibrium claim | S11 |
| Numerical and aggregation experiments | Matrix/ODE, marginals counterexample, refinement | Main Figs 1–2; S12 |
| Probabilistic screening code | Historical 200-draw summary; no recoverable draws or calendar PSA | Discussion; audit archive |
| Historical reports and manuscript versions | Archival lineage, not independent datasets | Source inventory and audit archive |

Coverage of the project research strands and their evidential roles.

</div>

## S5. Complete application input ledger and source status

The exact constants below are computational inputs. A verified reference identity does not establish that it supports the exact number, population or conversion. Births and mortality are frozen demographic scenario inputs; source datasets are potential calibration targets [UNWPP,WorldBankMortality](reference_key_map.csv). Transferred incidence and sequela sources provide context [Wahl2018,PhilippinesHTA,Dilok2019,Edmond2010](reference_key_map.csv). Hospital AOM detection motivates alternatives rather than verifying the attribution proxy [Do2025](reference_key_map.csv). Dose/carriage/waning parameters are structural. Utility and delivery references motivate valuation but do not validate local costs [Oh1996,Oostenbrink2002,Bennett2000,vanHoek2012,Mai2022](reference_key_map.csv). No independent inspection of patient-level hospital cost or survey data was possible; their retained numerical aggregates are unverified illustrative constants.

<div id="tab:epi_ledger">

| Parameter | Value | Evidence / interpretation |
|:---|:---|:---|
| Parameter | Value | Evidence / interpretation |
| birth_cohort | 1384610 | Structural / illustrative; not independently calibrated |
| u5_mortality_rate | 20 | Frozen all-cause inputs; reconciled by declared allocation |
| neonatal_mortality_rate | 10.3 | Frozen all-cause inputs; reconciled by declared allocation |
| inc_pneumonia | 1771 | Transferred / inherited; exact national estimate unverified |
| inc_meningitis | 12 | Transferred / inherited; exact national estimate unverified |
| inc_sepsis | 124 | Transferred / inherited; exact national estimate unverified |
| inc_aom_all | 9570 | Transferred / inherited; exact national estimate unverified |
| prop_pneumo_severe | 0.389 | Structural / illustrative; not independently calibrated |
| prop_aom_spn | 0.303 | Unverified attribution proxy; primary detection denominators differ |
| cfr_pneumo_severe | 0.014241 | Transferred / inherited; exact national estimate unverified |
| cfr_meningitis | 0.032258 | Transferred / inherited; exact national estimate unverified |
| cfr_sepsis | 0.08 | Transferred / inherited; exact national estimate unverified |
| prob_cog | 0.016 | Transferred / inherited; exact national estimate unverified |
| prob_epi | 0.022 | Transferred / inherited; exact national estimate unverified |
| prob_hemi | 0.032 | Transferred / inherited; exact national estimate unverified |
| prob_hear | 0.046 | Transferred / inherited; exact national estimate unverified |
| prob_cochlear_implant | 0.025 | Transferred / inherited; exact national estimate unverified |
| prob_special_education | 0.2574 | Transferred / inherited; exact national estimate unverified |
| vaccine_coverage | 0.9 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-GSK / 0 | 0.59 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-GSK / 1 | 0.49 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-GSK / 2 | 0.39 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-GSK / 3 | 0.32 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-GSK / 4 | 0.28 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-GSK / 5 | 0.26 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV13 / 0 | 0.54 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV13 / 1 | 0.38 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV13 / 2 | 0.32 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV13 / 3 | 0.28643333 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV13 / 4 | 0.25286667 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV13 / 5 | 0.2193 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-SII / 0 | 0.2211 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-SII / 1 | 0.2211 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-SII / 2 | 0.2211 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-SII / 3 | 0.2211 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-SII / 4 | 0.2211 | Structural / illustrative; not independently calibrated |
| vaccine_irrs / PCV10-SII / 5 | 0.2211 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 0 | 0.08 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 1 | 0.074 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 2 | 0.078 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 3 | 0.082 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 4 | 0.086 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 5 | 0.091 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 6 | 0.094 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 7 | 0.092 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 8 | 0.087 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 9 | 0.083 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 10 | 0.078 | Structural / illustrative; not independently calibrated |
| birth_month_shares / 11 | 0.075 | Structural / illustrative; not independently calibrated |
| warmup_months | 60 | Structural / illustrative; not independently calibrated |
| record_months | 120 | Structural / illustrative; not independently calibrated |
| age_bins | 60 | Structural / illustrative; not independently calibrated |
| vaccine_schedule_months / 0 | 2 | Structural / illustrative; not independently calibrated |
| vaccine_schedule_months / 1 | 4 | Structural / illustrative; not independently calibrated |
| vaccine_schedule_months / 2 | 12 | Structural / illustrative; not independently calibrated |
| direct_protection_start_month | 4 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV10-GSK / ipd | 0.5 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV10-GSK / pneumonia | 0.22 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV10-GSK / aom | 0.53 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV13 / ipd | 0.602 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV13 / pneumonia | 0.19 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV13 / aom | 0.46 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV10-SII / ipd | 0.49383048 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV10-SII / pneumonia | 0.15586012 | Structural / illustrative; not independently calibrated |
| vaccine_direct_effectiveness / PCV10-SII / aom | 0.37734555 | Structural / illustrative; not independently calibrated |
| pcv10_sii_relative_ve_assumption | central_weighted_gmc | Retained summary multiplier; original row-level derivation unavailable |
| pcv10_sii_relative_ve_multiplier | 0.82031641 | Retained summary multiplier; original row-level derivation unavailable |
| routine_dose1_transition_probs / 2 | 0.65 | Structural / illustrative; not independently calibrated |
| routine_dose1_transition_probs / 3 | 0.57 | Structural / illustrative; not independently calibrated |
| routine_dose1_transition_probs / 4 | 0.67 | Structural / illustrative; not independently calibrated |
| routine_dose1_transition_probs / 5 | 0.6 | Structural / illustrative; not independently calibrated |
| routine_dose2_transition_probs / 4 | 0.7 | Structural / illustrative; not independently calibrated |
| routine_dose2_transition_probs / 5 | 0.45 | Structural / illustrative; not independently calibrated |
| routine_dose2_transition_probs / 6 | 0.3 | Structural / illustrative; not independently calibrated |
| routine_dose2_transition_probs / 7 | 0.2 | Structural / illustrative; not independently calibrated |
| routine_dose2_transition_probs / 8 | 0.1 | Structural / illustrative; not independently calibrated |
| routine_booster_transition_probs / 12 | 0.55 | Structural / illustrative; not independently calibrated |
| routine_booster_transition_probs / 13 | 0.35 | Structural / illustrative; not independently calibrated |
| routine_booster_transition_probs / 14 | 0.25 | Structural / illustrative; not independently calibrated |
| routine_booster_transition_probs / 15 | 0.15 | Structural / illustrative; not independently calibrated |
| routine_booster_transition_probs / 16 | 0.1 | Structural / illustrative; not independently calibrated |
| catchup_enabled | 0 | Structural / illustrative; not independently calibrated |
| catchup_campaign_months | 12 | Structural / illustrative; not independently calibrated |
| catchup_7_11_to_dose1_prob | 0.12 | Structural / illustrative; not independently calibrated |
| catchup_12_23_to_dose1_prob | 0.08 | Structural / illustrative; not independently calibrated |
| dose1_relative_ve | 0.45 | Structural / illustrative; not independently calibrated |
| dose2_relative_ve | 1 | Structural / illustrative; not independently calibrated |
| booster_relative_ve | 1.1 | Structural / illustrative; not independently calibrated |
| max_direct_ve | 0.95 | Structural / illustrative; not independently calibrated |
| primary_waning_by_age / 24 | 0.85 | Structural / illustrative; not independently calibrated |
| primary_waning_by_age / 36 | 0.7 | Structural / illustrative; not independently calibrated |
| booster_waning_by_age / 24 | 0.95 | Structural / illustrative; not independently calibrated |
| booster_waning_by_age / 36 | 0.8 | Structural / illustrative; not independently calibrated |
| booster_waning_by_age / 48 | 0.7 | Structural / illustrative; not independently calibrated |
| carriage_reduction_by_group / Program_Dose0 | 0 | Structural / illustrative; not independently calibrated |
| carriage_reduction_by_group / Program_Dose1 | 0.15 | Structural / illustrative; not independently calibrated |
| carriage_reduction_by_group / Program_Dose2 | 0.34 | Structural / illustrative; not independently calibrated |
| carriage_reduction_by_group / Program_Booster | 0.4 | Structural / illustrative; not independently calibrated |
| carriage_age_weights / 0_5 | 0 | Structural / illustrative; not independently calibrated |
| carriage_age_weights / 6_17 | 1 | Structural / illustrative; not independently calibrated |
| carriage_age_weights / 18_23 | 0.8 | Structural / illustrative; not independently calibrated |
| carriage_age_weights / 24_59 | 0.35 | Structural / illustrative; not independently calibrated |
| mature_effective_dose_share | 0.88 | Structural / illustrative; not independently calibrated |
| carriage_to_herd_fraction | 1 | Structural / illustrative; not independently calibrated |
| mortality_accounting | cause_deleted | Structural / illustrative; not independently calibrated |

Frozen epidemiological and programme constants. Probability and rate units follow the main Methods.

</div>

<div id="tab:econ_ledger">

| Parameter | Value | Status |
|:---|:---|:---|
| Parameter | Value | Status |
| birth_cohort | 1384610 | Illustrative / inherited; see S3 |
| coverage_rate / pcv10 | 0.9 | Illustrative / inherited; see S3 |
| coverage_rate / pcv10sii | 0.9 | Illustrative / inherited; see S3 |
| coverage_rate / pcv13 | 0.9 | Illustrative / inherited; see S3 |
| doses_per_child | 3 | Illustrative / inherited; see S3 |
| costs_medical / inpatient_pneumonia | 214.65 | Illustrative / inherited; see S3 |
| costs_medical / outpatient_pneumonia | 15.6214 | Illustrative / inherited; see S3 |
| costs_medical / inpatient_meningitis | 557.767 | Illustrative / inherited; see S3 |
| costs_medical / inpatient_sepsis | 572.427 | Illustrative / inherited; see S3 |
| costs_medical / outpatient_aom | 12.9223 | Illustrative / inherited; see S3 |
| costs_nonmedical / inpatient_pneumonia | 134.16 | Illustrative / inherited; see S3 |
| costs_nonmedical / outpatient_pneumonia | 17.8905 | Illustrative / inherited; see S3 |
| costs_nonmedical / inpatient_meningitis | 280.846 | Illustrative / inherited; see S3 |
| costs_nonmedical / inpatient_sepsis | 287.718 | Illustrative / inherited; see S3 |
| costs_nonmedical / outpatient_aom | 19.9925 | Illustrative / inherited; see S3 |
| costs_caregiver / inpatient_pneumonia | 173.673 | Illustrative / inherited; see S3 |
| costs_caregiver / outpatient_pneumonia | 117.344 | Illustrative / inherited; see S3 |
| costs_caregiver / inpatient_meningitis | 280.724 | Illustrative / inherited; see S3 |
| costs_caregiver / inpatient_sepsis | 182.439 | Illustrative / inherited; see S3 |
| costs_caregiver / outpatient_aom | 131.134 | Illustrative / inherited; see S3 |
| costs_sequelae / cognitive_annual | 71.2233 | Illustrative / inherited; see S3 |
| costs_sequelae / hearing_case | 990.291 | Illustrative / inherited; see S3 |
| costs_sequelae / epilepsy_case | 319.476 | Illustrative / inherited; see S3 |
| costs_sequelae / hemiplegia_case | 969.393 | Illustrative / inherited; see S3 |
| costs_sequelae / cochlear_implant_case | 24106.2 | Illustrative / inherited; see S3 |
| costs_sequelae / special_education_annual | 2384.88 | Illustrative / inherited; see S3 |
| vaccine_prices / pcv10 | 1.94175 | Illustrative / inherited; see S3 |
| vaccine_prices / pcv13 | 4.22 | Illustrative / inherited; see S3 |
| vaccine_prices / pcv10sii | 1.94175 | Illustrative / inherited; see S3 |
| service_delivery_per_dose | 2.038835 | Illustrative / inherited; see S3 |
| wastage_rate | 0.05 | Illustrative / inherited; see S3 |
| ae_incidence_per_100k_doses | 4.01 | Illustrative / inherited; see S3 |
| ae_cost_per_case | 4800.355 | Illustrative / inherited; see S3 |
| utilities / inpatient_pneumonia | 0.9941 | Illustrative / inherited; see S3 |
| utilities / outpatient_pneumonia | 0.9963 | Illustrative / inherited; see S3 |
| utilities / meningitis | 0.9768 | Illustrative / inherited; see S3 |
| utilities / sepsis | 0.9921 | Illustrative / inherited; see S3 |
| utilities / aom | 0.995 | Illustrative / inherited; see S3 |
| utilities / cognitive | 0.62 | Illustrative / inherited; see S3 |
| utilities / hearing | 0.91 | Illustrative / inherited; see S3 |
| utilities / epilepsy | 0.83 | Illustrative / inherited; see S3 |
| utilities / hemiplegia | 0.3903 | Illustrative / inherited; see S3 |
| life_expectancy | 74 | Illustrative / inherited; see S3 |
| discount_rate | 0.03 | Illustrative / inherited; see S3 |
| discount_rate_align | 0.03 | Illustrative / inherited; see S3 |
| prod_loss_per_death_age0 | 41158.175 | Illustrative / inherited; see S3 |

Economic constants. Scenario USD values use the disclosed conversions and source-year exceptions, utilities are dimensionless, durations are days, and the discount rate is annual.

</div>

## S6. Closed-cohort checks and AOM denominator convention

The cohort count model and deterministic expectation implementation share a discrete monthly event kernel. In the saved replicate summaries, the death reward is the conditional expectation obtained from sampled disease events multiplied by endpoint fatality probabilities, whereas actual sampled deaths update the alive population. The variation in the saved death reward must therefore not be interpreted as a distribution of integer fatalities. Its 2024-valued economics and mature-effect proxies differ from the current calendar ledger; its results are not pooled with calendar outcomes. The historical stored 1,000-replicate means, standard deviations and pointwise Monte Carlo intervals were independently reconstructed from replicate-level summaries. The seeded four-configuration audit was generated for the preceding methods version; this revision independently reconstructs its means and intervals from the stored replicates. AOM conventions separately count only the operational attributed component or retain the unaffected all-cause component. The difference in total cases is an outcome-definition difference; it does not validate clinical vaccine effects or constitute exact replication of an unavailable external implementation.

<div id="tab:cohort">

| AOM mode / strategy | Metric | Deterministic | MC mean | MC lower | MC upper |
|:---|:---|---:|---:|---:|---:|
| AOM mode / strategy | Metric | Deterministic | MC mean | MC lower | MC upper |
| local_spn_only / no_vaccination | total_cases | 326,304.82 | 326,312.01 | 326,275.56 | 326,348.45 |
| local_spn_only / no_vaccination | total_deaths | 1,369.80 | 1,370.14 | 1,369.67 | 1,370.61 |
| local_spn_only / no_vaccination | pv_total_cost | 134,498,639.16 | 134,519,031.64 | 134,491,351.50 | 134,546,711.77 |
| local_spn_only / no_vaccination | pv_total_qaly_loss | 39,426.65 | 39,436.22 | 39,422.70 | 39,449.75 |
| local_spn_only / mature_pcv13 | total_cases | 71,593.98 | 71,601.05 | 71,584.33 | 71,617.76 |
| local_spn_only / mature_pcv13 | total_deaths | 300.51 | 300.65 | 300.43 | 300.87 |
| local_spn_only / mature_pcv13 | pv_total_cost | 54,883,578.18 | 54,892,588.31 | 54,879,940.68 | 54,905,235.95 |
| local_spn_only / mature_pcv13 | pv_total_qaly_loss | 8,649.48 | 8,653.67 | 8,647.39 | 8,659.95 |
| reference_aligned / no_vaccination | total_cases | 778,232.73 | 778,235.72 | 778,181.32 | 778,290.11 |
| reference_aligned / no_vaccination | total_deaths | 1,369.80 | 1,369.54 | 1,369.07 | 1,370.02 |
| reference_aligned / no_vaccination | pv_total_cost | 207,333,316.64 | 207,319,740.52 | 207,290,837.12 | 207,348,643.92 |
| reference_aligned / no_vaccination | pv_total_qaly_loss | 39,466.87 | 39,459.49 | 39,445.89 | 39,473.10 |
| reference_aligned / mature_pcv13 | total_cases | 523,778.58 | 523,789.89 | 523,744.56 | 523,835.22 |
| reference_aligned / mature_pcv13 | total_deaths | 300.51 | 300.62 | 300.39 | 300.84 |
| reference_aligned / mature_pcv13 | pv_total_cost | 127,759,000.41 | 127,767,790.91 | 127,752,956.58 | 127,782,625.23 |
| reference_aligned / mature_pcv13 | pv_total_qaly_loss | 8,689.72 | 8,692.97 | 8,686.50 | 8,699.44 |

Retained seeded closed-cohort expectation audit: 1,000 replications for each mode/strategy configuration. Intervals are for simulated means.

</div>

Historical cross-model summaries below preserve population and accounting boundaries. They were not rerun as comparable factorial configurations; the old two-state comparison also changed birth and dose-cost conventions. They are retained to show the transition from cohort to calendar estimands.

<div id="tab:supp_baseline">

| **Model** | **Estimand** | **Cases** | **Deaths** | **PV cost** | **PV QALY loss** | **Births** |
|:---|:---|:---|:---|:---|:---|:---|
| **Model** | **Estimand** | **Cases** | **Deaths** | **PV cost** | **PV QALY loss** | **Births** |
| Closed-cohort ODE | One birth cohort, 5 years | 325,805 | 1,369 | 129.90 | 39,411.4 | 1,384,610 |
| Closed-cohort Markov | One birth cohort, 5 years | 325,374 | 1,367 | 129.73 | 39,357.4 | 1,384,610 |
| Calendar PDE, on-schedule protection | Under-five calendar population, 10 years | 3,258,100 | 13,690 | 1,197.87 | 353,979.6 | 13,846,100 |
| Calendar PDE, delayed doses and waning | Under-five calendar population, 10 years | 3,258,100 | 13,690 | 1,197.87 | 353,979.6 | 13,846,100 |

Baseline model totals by structure. Present-value costs are reported in million 2023 US dollars.

</div>

<div id="tab:supp_scenario_totals">

| **Model** | **Scenario** | **Vaccine** | **Cases** | **Deaths** | **Doses** | **PV cost** | **PV QALY loss** |
|:---|:---|:---|:---|:---|:---|:---|:---|
| **Model** | **Scenario** | **Vaccine** | **Cases** | **Deaths** | **Doses** | **PV cost** | **PV QALY loss** |
| Closed-cohort ODE | Dynamic rollout IRR | PCV10-GSK | 145,980 | 695 | – | 79.51 | 20,155.7 |
| Closed-cohort ODE | Immediate mature IRR | PCV10-GSK | 84,686 | 358 | – | 50.11 | 10,303.4 |
| Closed-cohort ODE | Dynamic rollout IRR | PCV13 | 126,159 | 603 | – | 79.98 | 17,509.3 |
| Closed-cohort ODE | Immediate mature IRR | PCV13 | 71,430 | 301 | – | 53.65 | 8,665.7 |
| Closed-cohort ODE | Immediate mature IRR | PCV10-SII | 72,018 | 303 | – | 44.99 | 8,724.0 |
| Closed-cohort Markov | Dynamic rollout IRR | PCV10-GSK | 145,868 | 694 | – | 79.44 | 20,127.4 |
| Closed-cohort Markov | Immediate mature IRR | PCV10-GSK | 84,616 | 353 | – | 49.88 | 10,165.7 |
| Closed-cohort Markov | Dynamic rollout IRR | PCV13 | 126,046 | 603 | – | 79.96 | 17,508.2 |
| Closed-cohort Markov | Immediate mature IRR | PCV13 | 71,370 | 299 | – | 53.56 | 8,606.7 |
| Closed-cohort Markov | Immediate mature IRR | PCV10-SII | 71,950 | 303 | – | 44.97 | 8,721.5 |
| Calendar PDE, on-schedule protection | Rollout-to-stable programme | PCV10-GSK | 1,095,087 | 4,572 | 35,149,439 | 539.31 | 120,644.9 |
| Calendar PDE, on-schedule protection | Rollout-to-stable programme | PCV13 | 934,879 | 3,912 | 35,149,439 | 552.49 | 103,212.1 |
| Calendar PDE, on-schedule protection | Rollout-to-stable programme | PCV10-SII | 720,136 | 3,037 | 35,149,439 | 394.26 | 78,240.8 |
| Calendar PDE, delayed doses and waning | Rollout-to-stable programme | PCV10-GSK | 1,448,251 | 5,843 | 30,046,011 | 650.46 | 157,022.8 |
| Calendar PDE, delayed doses and waning | Rollout-to-stable programme | PCV13 | 1,299,799 | 5,330 | 30,046,281 | 664.06 | 143,817.1 |
| Calendar PDE, delayed doses and waning | Rollout-to-stable programme | PCV10-SII | 1,194,826 | 4,916 | 30,046,573 | 562.71 | 132,411.1 |

Vaccine-scenario total outcomes by model structure. Present-value costs are reported in million 2023 US dollars.

</div>

<div id="tab:supp_full_qaly">

| **Model** | **Scenario** | **Vaccine** | **Cases averted** | **Deaths averted** | **QALYs gained** | **Cost diff.** | **USD / QALY** |
|:---|:---|:---|:---|:---|:---|:---|:---|
| **Model** | **Scenario** | **Vaccine** | **Cases averted** | **Deaths averted** | **QALYs gained** | **Cost diff.** | **USD / QALY** |
| Closed-cohort ODE | Dynamic rollout IRR | PCV10-GSK | 12,987.4 | 48.7 | 1,390.7 | -3.64 | -2,617 |
| Closed-cohort ODE | Immediate mature IRR | PCV10-GSK | 17,414.2 | 73.0 | 2,102.2 | -5.76 | -2,741 |
| Closed-cohort ODE | Dynamic rollout IRR | PCV13 | 14,418.9 | 55.3 | 1,581.8 | -3.61 | -2,279 |
| Closed-cohort ODE | Immediate mature IRR | PCV13 | 18,371.6 | 77.1 | 2,220.5 | -5.51 | -2,480 |
| Closed-cohort ODE | Immediate mature IRR | PCV10-SII | 18,329.1 | 77.0 | 2,216.3 | -6.13 | -2,767 |
| Closed-cohort Markov | Dynamic rollout IRR | PCV10-GSK | 12,964.3 | 48.6 | 1,388.8 | -3.63 | -2,615 |
| Closed-cohort Markov | Immediate mature IRR | PCV10-GSK | 17,388.1 | 73.2 | 2,108.3 | -5.77 | -2,735 |
| Closed-cohort Markov | Dynamic rollout IRR | PCV13 | 14,396.0 | 55.2 | 1,578.0 | -3.59 | -2,278 |
| Closed-cohort Markov | Immediate mature IRR | PCV13 | 18,344.8 | 77.1 | 2,220.9 | -5.50 | -2,477 |
| Closed-cohort Markov | Immediate mature IRR | PCV10-SII | 18,303.0 | 76.8 | 2,212.6 | -6.12 | -2,767 |
| Calendar PDE, on-schedule protection | Rollout-to-stable programme | PCV10-GSK | 15,621.8 | 65.9 | 1,685.2 | -4.76 | -2,822 |
| Calendar PDE, on-schedule protection | Rollout-to-stable programme | PCV13 | 16,778.9 | 70.6 | 1,811.1 | -4.66 | -2,574 |
| Calendar PDE, on-schedule protection | Rollout-to-stable programme | PCV10-SII | 18,329.8 | 76.9 | 1,991.5 | -5.80 | -2,914 |
| Calendar PDE, delayed doses and waning | Rollout-to-stable programme | PCV10-GSK | 13,071.2 | 56.7 | 1,422.5 | -3.95 | -2,779 |
| Calendar PDE, delayed doses and waning | Rollout-to-stable programme | PCV13 | 14,143.3 | 60.4 | 1,517.8 | -3.86 | -2,540 |
| Calendar PDE, delayed doses and waning | Rollout-to-stable programme | PCV10-SII | 14,901.5 | 63.4 | 1,600.2 | -4.59 | -2,867 |

Complete incremental health-economic outcomes standardised per 100,000 births. Cost differences are reported in million 2023 US dollars.

</div>

<div id="tab:supp_refinement">

| **Contrast** | **Vaccine** | **Scenario** | **Cases change** | **Deaths change** | **QALY change** | **Cost change, USD M** | **ICER shift** |
|:---|:---|:---|:---|:---|:---|:---|:---|
| **Contrast** | **Vaccine** | **Scenario** | **Cases change** | **Deaths change** | **QALY change** | **Cost change, USD M** | **ICER shift** |
| ODE to Markov cohort check | PCV10-GSK | Dynamic rollout IRR | -23.1 | -0.1 | -1.9 | +0.01 | +1.8 |
| ODE to Markov cohort check | PCV10-GSK | Immediate mature IRR | -26.1 | +0.2 | +6.1 | -0.00 | +6.0 |
| ODE to Markov cohort check | PCV13 | Dynamic rollout IRR | -22.9 | -0.1 | -3.8 | +0.01 | +1.2 |
| ODE to Markov cohort check | PCV13 | Immediate mature IRR | -26.8 | 0.0 | +0.4 | +0.01 | +3.0 |
| ODE to Markov cohort check | PCV10-SII | Immediate mature IRR | -26.2 | -0.1 | -3.7 | +0.01 | +0.3 |
| Programme realism added to calendar PDE | PCV10-GSK | Rollout-to-stable programme | -2,550.6 | -9.2 | -262.7 | +0.80 | +43.1 |
| Programme realism added to calendar PDE | PCV13 | Rollout-to-stable programme | -2,635.5 | -10.2 | -293.3 | +0.81 | +33.6 |
| Programme realism added to calendar PDE | PCV10-SII | Rollout-to-stable programme | -3,428.3 | -13.6 | -391.2 | +1.22 | +47.8 |

Refinement-impact contrasts.

</div>

<div id="tab:supp_stable">

| **Vaccine** | **Endpoint** | **Baseline count** | **Vaccine count** | **Averted count** | **Reduction** |
|:---|:---|:---|:---|:---|:---|
| **Vaccine** | **Endpoint** | **Baseline count** | **Vaccine count** | **Averted count** | **Reduction** |
| PCV10-GSK | Meningitis | 818 | 237 | 581 | 71.0% |
| PCV10-GSK | Severe pneumonia | 46,914 | 12,518 | 34,396 | 73.3% |
| PCV10-GSK | Non-severe pneumonia | 73,687 | 19,663 | 54,024 | 73.3% |
| PCV10-GSK | Sepsis | 8,444 | 2,380 | 6,064 | 71.8% |
| PCV10-GSK | Acute otitis media | 195,947 | 55,331 | 140,616 | 71.8% |
| PCV10-GSK | Disease deaths | 1,369 | 376 | 993 | 72.5% |
| PCV13 | Meningitis | 818 | 208 | 610 | 74.6% |
| PCV13 | Severe pneumonia | 46,914 | 10,514 | 36,400 | 77.6% |
| PCV13 | Non-severe pneumonia | 73,687 | 16,514 | 57,173 | 77.6% |
| PCV13 | Sepsis | 8,444 | 2,067 | 6,377 | 75.5% |
| PCV13 | Acute otitis media | 195,947 | 45,896 | 150,051 | 76.6% |
| PCV13 | Disease deaths | 1,369 | 322 | 1,047 | 76.5% |
| PCV10-SII | Meningitis | 818 | 201 | 617 | 75.4% |
| PCV10-SII | Severe pneumonia | 46,914 | 10,552 | 36,362 | 77.5% |
| PCV10-SII | Non-severe pneumonia | 73,687 | 16,575 | 57,112 | 77.5% |
| PCV10-SII | Sepsis | 8,444 | 2,020 | 6,424 | 76.1% |
| PCV10-SII | Acute otitis media | 195,947 | 45,513 | 150,434 | 76.8% |
| PCV10-SII | Disease deaths | 1,369 | 318 | 1,051 | 76.8% |

Stable annual endpoint reductions in the programme-realism calendar PDE.

</div>

## S7. Complete factorial, interactions and herd counterfactuals

<div id="tab:allfactorial">

| Product / DWG   |     Cases | Deaths |      QALY | Cost diff. |      ICER |
|:----------------|----------:|-------:|----------:|-----------:|----------:|
| Product / DWG   |     Cases | Deaths |      QALY | Cost diff. |      ICER |
| PCV10-GSK / 000 | 14,609.53 | 62.709 | 1,588.599 |    -4.3970 | -2,767.82 |
| PCV10-GSK / 001 | 13,975.89 | 59.475 | 1,495.785 |    -4.0931 | -2,736.44 |
| PCV10-GSK / 010 | 14,237.41 | 62.205 | 1,576.151 |    -4.3223 | -2,742.30 |
| PCV10-GSK / 011 | 13,602.81 | 58.969 | 1,483.257 |    -4.0182 | -2,709.01 |
| PCV10-GSK / 100 | 14,267.83 | 61.362 | 1,553.801 |    -4.4155 | -2,841.76 |
| PCV10-GSK / 101 | 13,451.04 | 57.292 | 1,437.867 |    -4.0346 | -2,805.96 |
| PCV10-GSK / 110 | 13,950.77 | 60.931 | 1,543.141 |    -4.3518 | -2,820.09 |
| PCV10-GSK / 111 | 13,131.33 | 56.854 | 1,427.015 |    -3.9701 | -2,782.08 |
| PCV13 / 000     | 16,083.81 | 67.454 | 1,713.288 |    -4.3544 | -2,541.53 |
| PCV13 / 001     | 14,921.70 | 63.435 | 1,597.914 |    -3.9154 | -2,450.32 |
| PCV13 / 010     | 15,828.49 | 66.912 | 1,699.871 |    -4.2939 | -2,526.03 |
| PCV13 / 011     | 14,665.24 | 62.889 | 1,584.403 |    -3.8546 | -2,432.83 |
| PCV13 / 100     | 15,844.32 | 66.016 | 1,676.106 |    -4.4592 | -2,660.45 |
| PCV13 / 101     | 14,428.82 | 61.044 | 1,534.373 |    -3.9251 | -2,558.08 |
| PCV13 / 110     | 15,626.75 | 65.551 | 1,664.615 |    -4.4076 | -2,647.80 |
| PCV13 / 111     | 14,208.26 | 60.572 | 1,522.665 |    -3.8725 | -2,543.26 |
| PCV10-SII / 000 | 18,190.67 | 76.218 | 1,966.395 |    -5.7286 | -2,913.27 |
| PCV10-SII / 001 | 15,640.00 | 66.183 | 1,674.940 |    -4.7100 | -2,812.06 |
| PCV10-SII / 010 | 18,007.00 | 75.847 | 1,957.255 |    -5.6862 | -2,905.18 |
| PCV10-SII / 011 | 15,455.10 | 65.807 | 1,665.688 |    -4.6672 | -2,801.94 |
| PCV10-SII / 100 | 18,031.29 | 75.339 | 1,944.052 |    -5.7980 | -2,982.44 |
| PCV10-SII / 101 | 15,128.65 | 63.882 | 1,613.026 |    -4.6425 | -2,878.13 |
| PCV10-SII / 110 | 17,874.85 | 75.021 | 1,936.228 |    -5.7618 | -2,975.78 |
| PCV10-SII / 111 | 14,969.07 | 63.555 | 1,604.952 |    -4.6053 | -2,869.43 |

Complete monthly factorial outcomes versus NV per 100,000 births; cost differences in million USD.

</div>

All within-product mechanism changes below use the five-state 000 reference. An absolute change is $`Y(z)-Y(000)`$ and its relative change is $`100\{Y(z)-Y(000)\}/|Y(000)|`$. For incremental costs a positive change means less saving; for a negative ICER it is a descriptive ratio change. The baseline is product-specific.

<div id="tab:factorialchangecases_averted_per100k">

| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
|:----------------|----------------:|------------------:|-----------:|
| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
| PCV10-GSK / 000 |      14,609.530 |             0.000 |      0.000 |
| PCV10-GSK / 001 |      13,975.886 |          -633.644 |     -4.337 |
| PCV10-GSK / 010 |      14,237.411 |          -372.119 |     -2.547 |
| PCV10-GSK / 011 |      13,602.814 |        -1,006.716 |     -6.891 |
| PCV10-GSK / 100 |      14,267.830 |          -341.700 |     -2.339 |
| PCV10-GSK / 101 |      13,451.036 |        -1,158.494 |     -7.930 |
| PCV10-GSK / 110 |      13,950.772 |          -658.758 |     -4.509 |
| PCV10-GSK / 111 |      13,131.328 |        -1,478.202 |    -10.118 |
| PCV13 / 000     |      16,083.812 |             0.000 |      0.000 |
| PCV13 / 001     |      14,921.700 |        -1,162.112 |     -7.225 |
| PCV13 / 010     |      15,828.495 |          -255.317 |     -1.587 |
| PCV13 / 011     |      14,665.241 |        -1,418.571 |     -8.820 |
| PCV13 / 100     |      15,844.319 |          -239.493 |     -1.489 |
| PCV13 / 101     |      14,428.817 |        -1,654.995 |    -10.290 |
| PCV13 / 110     |      15,626.746 |          -457.066 |     -2.842 |
| PCV13 / 111     |      14,208.263 |        -1,875.549 |    -11.661 |
| PCV10-SII / 000 |      18,190.672 |             0.000 |      0.000 |
| PCV10-SII / 001 |      15,640.004 |        -2,550.668 |    -14.022 |
| PCV10-SII / 010 |      18,007.004 |          -183.669 |     -1.010 |
| PCV10-SII / 011 |      15,455.097 |        -2,735.576 |    -15.038 |
| PCV10-SII / 100 |      18,031.285 |          -159.387 |     -0.876 |
| PCV10-SII / 101 |      15,128.650 |        -3,062.022 |    -16.833 |
| PCV10-SII / 110 |      17,874.846 |          -315.827 |     -1.736 |
| PCV10-SII / 111 |      14,969.072 |        -3,221.600 |    -17.710 |

Within-product absolute and relative mechanism contrasts: Cases averted. Outcomes other than ICER are per 100,000 recorded-window births.

</div>

<div id="tab:factorialchangedeaths_averted_per100k">

| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
|:----------------|----------------:|------------------:|-----------:|
| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
| PCV10-GSK / 000 |          62.709 |             0.000 |      0.000 |
| PCV10-GSK / 001 |          59.475 |            -3.233 |     -5.156 |
| PCV10-GSK / 010 |          62.205 |            -0.504 |     -0.803 |
| PCV10-GSK / 011 |          58.969 |            -3.740 |     -5.964 |
| PCV10-GSK / 100 |          61.362 |            -1.346 |     -2.147 |
| PCV10-GSK / 101 |          57.292 |            -5.417 |     -8.638 |
| PCV10-GSK / 110 |          60.931 |            -1.777 |     -2.834 |
| PCV10-GSK / 111 |          56.854 |            -5.855 |     -9.337 |
| PCV13 / 000     |          67.454 |             0.000 |      0.000 |
| PCV13 / 001     |          63.435 |            -4.019 |     -5.958 |
| PCV13 / 010     |          66.912 |            -0.543 |     -0.804 |
| PCV13 / 011     |          62.889 |            -4.565 |     -6.768 |
| PCV13 / 100     |          66.016 |            -1.439 |     -2.133 |
| PCV13 / 101     |          61.044 |            -6.410 |     -9.503 |
| PCV13 / 110     |          65.551 |            -1.903 |     -2.822 |
| PCV13 / 111     |          60.572 |            -6.883 |    -10.203 |
| PCV10-SII / 000 |          76.218 |             0.000 |      0.000 |
| PCV10-SII / 001 |          66.183 |           -10.035 |    -13.167 |
| PCV10-SII / 010 |          75.847 |            -0.372 |     -0.487 |
| PCV10-SII / 011 |          65.807 |           -10.411 |    -13.659 |
| PCV10-SII / 100 |          75.339 |            -0.879 |     -1.154 |
| PCV10-SII / 101 |          63.882 |           -12.336 |    -16.185 |
| PCV10-SII / 110 |          75.021 |            -1.197 |     -1.571 |
| PCV10-SII / 111 |          63.555 |           -12.663 |    -16.614 |

Within-product absolute and relative mechanism contrasts: Pneumococcal deaths averted. Outcomes other than ICER are per 100,000 recorded-window births.

</div>

<div id="tab:factorialchangeqalys_gained_per100k">

| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
|:----------------|----------------:|------------------:|-----------:|
| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
| PCV10-GSK / 000 |       1,588.599 |             0.000 |      0.000 |
| PCV10-GSK / 001 |       1,495.785 |           -92.814 |     -5.842 |
| PCV10-GSK / 010 |       1,576.151 |           -12.448 |     -0.784 |
| PCV10-GSK / 011 |       1,483.257 |          -105.341 |     -6.631 |
| PCV10-GSK / 100 |       1,553.801 |           -34.798 |     -2.190 |
| PCV10-GSK / 101 |       1,437.867 |          -150.732 |     -9.488 |
| PCV10-GSK / 110 |       1,543.141 |           -45.458 |     -2.862 |
| PCV10-GSK / 111 |       1,427.015 |          -161.584 |    -10.171 |
| PCV13 / 000     |       1,713.288 |             0.000 |      0.000 |
| PCV13 / 001     |       1,597.914 |          -115.374 |     -6.734 |
| PCV13 / 010     |       1,699.871 |           -13.417 |     -0.783 |
| PCV13 / 011     |       1,584.403 |          -128.885 |     -7.523 |
| PCV13 / 100     |       1,676.106 |           -37.182 |     -2.170 |
| PCV13 / 101     |       1,534.373 |          -178.915 |    -10.443 |
| PCV13 / 110     |       1,664.615 |           -48.673 |     -2.841 |
| PCV13 / 111     |       1,522.665 |          -190.623 |    -11.126 |
| PCV10-SII / 000 |       1,966.395 |             0.000 |      0.000 |
| PCV10-SII / 001 |       1,674.940 |          -291.455 |    -14.822 |
| PCV10-SII / 010 |       1,957.255 |            -9.140 |     -0.465 |
| PCV10-SII / 011 |       1,665.688 |          -300.708 |    -15.292 |
| PCV10-SII / 100 |       1,944.052 |           -22.343 |     -1.136 |
| PCV10-SII / 101 |       1,613.026 |          -353.369 |    -17.970 |
| PCV10-SII / 110 |       1,936.228 |           -30.167 |     -1.534 |
| PCV10-SII / 111 |       1,604.952 |          -361.443 |    -18.381 |

Within-product absolute and relative mechanism contrasts: QALYs gained. Outcomes other than ICER are per 100,000 recorded-window births.

</div>

<div id="tab:factorialchangeincremental_costs_per100k">

| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
|:----------------|----------------:|------------------:|-----------:|
| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
| PCV10-GSK / 000 |          -4.397 |             0.000 |      0.000 |
| PCV10-GSK / 001 |          -4.093 |             0.304 |      6.910 |
| PCV10-GSK / 010 |          -4.322 |             0.075 |      1.698 |
| PCV10-GSK / 011 |          -4.018 |             0.379 |      8.615 |
| PCV10-GSK / 100 |          -4.416 |            -0.019 |     -0.423 |
| PCV10-GSK / 101 |          -4.035 |             0.362 |      8.241 |
| PCV10-GSK / 110 |          -4.352 |             0.045 |      1.027 |
| PCV10-GSK / 111 |          -3.970 |             0.427 |      9.709 |
| PCV13 / 000     |          -4.354 |             0.000 |      0.000 |
| PCV13 / 001     |          -3.915 |             0.439 |     10.081 |
| PCV13 / 010     |          -4.294 |             0.060 |      1.388 |
| PCV13 / 011     |          -3.855 |             0.500 |     11.478 |
| PCV13 / 100     |          -4.459 |            -0.105 |     -2.407 |
| PCV13 / 101     |          -3.925 |             0.429 |      9.860 |
| PCV13 / 110     |          -4.408 |            -0.053 |     -1.221 |
| PCV13 / 111     |          -3.873 |             0.482 |     11.066 |
| PCV10-SII / 000 |          -5.729 |             0.000 |      0.000 |
| PCV10-SII / 001 |          -4.710 |             1.019 |     17.781 |
| PCV10-SII / 010 |          -5.686 |             0.042 |      0.741 |
| PCV10-SII / 011 |          -4.667 |             1.061 |     18.529 |
| PCV10-SII / 100 |          -5.798 |            -0.069 |     -1.211 |
| PCV10-SII / 101 |          -4.643 |             1.086 |     18.960 |
| PCV10-SII / 110 |          -5.762 |            -0.033 |     -0.579 |
| PCV10-SII / 111 |          -4.605 |             1.123 |     19.609 |

Within-product absolute and relative mechanism contrasts: Incremental cost (million scenario USD). Outcomes other than ICER are per 100,000 recorded-window births.

</div>

<div id="tab:factorialchangeicer">

| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
|:----------------|----------------:|------------------:|-----------:|
| Product / DWG   | Value versus NV | Change versus 000 | Change (%) |
| PCV10-GSK / 000 |      -2,767.817 |             0.000 |      0.000 |
| PCV10-GSK / 001 |      -2,736.444 |            31.372 |      1.133 |
| PCV10-GSK / 010 |      -2,742.303 |            25.514 |      0.922 |
| PCV10-GSK / 011 |      -2,709.010 |            58.806 |      2.125 |
| PCV10-GSK / 100 |      -2,841.759 |           -73.943 |     -2.672 |
| PCV10-GSK / 101 |      -2,805.957 |           -38.140 |     -1.378 |
| PCV10-GSK / 110 |      -2,820.088 |           -52.271 |     -1.889 |
| PCV10-GSK / 111 |      -2,782.081 |           -14.264 |     -0.515 |
| PCV13 / 000     |      -2,541.535 |             0.000 |      0.000 |
| PCV13 / 001     |      -2,450.320 |            91.215 |      3.589 |
| PCV13 / 010     |      -2,526.031 |            15.504 |      0.610 |
| PCV13 / 011     |      -2,432.832 |           108.703 |      4.277 |
| PCV13 / 100     |      -2,660.449 |          -118.914 |     -4.679 |
| PCV13 / 101     |      -2,558.082 |           -16.547 |     -0.651 |
| PCV13 / 110     |      -2,647.795 |          -106.260 |     -4.181 |
| PCV13 / 111     |      -2,543.262 |            -1.727 |     -0.068 |
| PCV10-SII / 000 |      -2,913.270 |             0.000 |      0.000 |
| PCV10-SII / 001 |      -2,812.055 |           101.215 |      3.474 |
| PCV10-SII / 010 |      -2,905.180 |             8.090 |      0.278 |
| PCV10-SII / 011 |      -2,801.938 |           111.333 |      3.822 |
| PCV10-SII / 100 |      -2,982.445 |           -69.174 |     -2.374 |
| PCV10-SII / 101 |      -2,878.132 |            35.138 |      1.206 |
| PCV10-SII / 110 |      -2,975.776 |           -62.506 |     -2.146 |
| PCV10-SII / 111 |      -2,869.426 |            43.844 |      1.505 |

Within-product absolute and relative mechanism contrasts: ICER (scenario USD/QALY). Outcomes other than ICER are per 100,000 recorded-window births.

</div>

<div id="tab:interactioncases_averted_per100k">

| Product   | Term | Order | Coefficient |
|:----------|:-----|------:|------------:|
| Product   | Term | Order | Coefficient |
| PCV10-GSK | G    |     1 |    -633.644 |
| PCV10-GSK | W    |     1 |    -372.119 |
| PCV10-GSK | WG   |     2 |   -0.953339 |
| PCV10-GSK | D    |     1 |      -341.7 |
| PCV10-GSK | DG   |     2 |     -183.15 |
| PCV10-GSK | DW   |     2 |     55.0607 |
| PCV10-GSK | DWG  |     3 |    -1.69673 |
| PCV13     | G    |     1 |    -1162.11 |
| PCV13     | W    |     1 |    -255.317 |
| PCV13     | WG   |     2 |    -1.14189 |
| PCV13     | D    |     1 |    -239.493 |
| PCV13     | DG   |     2 |     -253.39 |
| PCV13     | DW   |     2 |     37.7441 |
| PCV13     | DWG  |     3 |    -1.83889 |
| PCV10-SII | G    |     1 |    -2550.67 |
| PCV10-SII | W    |     1 |    -183.669 |
| PCV10-SII | WG   |     2 |    -1.23887 |
| PCV10-SII | D    |     1 |    -159.387 |
| PCV10-SII | DG   |     2 |    -351.967 |
| PCV10-SII | DW   |     2 |      27.229 |
| PCV10-SII | DWG  |     3 |    -1.89922 |

Factorial interaction coefficients for cases_averted_per100k. Raw coefficients retain output units; ICER contrasts are descriptive ratio changes, without an additive component-benefit interpretation.

</div>

<div id="tab:interactiondeaths_averted_per100k">

| Product   | Term | Order | Coefficient |
|:----------|:-----|------:|------------:|
| Product   | Term | Order | Coefficient |
| PCV10-GSK | G    |     1 |    -3.23326 |
| PCV10-GSK | W    |     1 |   -0.503584 |
| PCV10-GSK | WG   |     2 | -0.00289071 |
| PCV10-GSK | D    |     1 |    -1.34631 |
| PCV10-GSK | DG   |     2 |   -0.837151 |
| PCV10-GSK | DW   |     2 |   0.0724152 |
| PCV10-GSK | DWG  |     3 | -0.00409383 |
| PCV13     | G    |     1 |    -4.01918 |
| PCV13     | W    |     1 |   -0.542608 |
| PCV13     | WG   |     2 | -0.00339339 |
| PCV13     | D    |     1 |    -1.43877 |
| PCV13     | DG   |     2 |   -0.952135 |
| PCV13     | DW   |     2 |   0.0779951 |
| PCV13     | DWG  |     3 | -0.00453946 |
| PCV10-SII | G    |     1 |    -10.0354 |
| PCV10-SII | W    |     1 |   -0.371561 |
| PCV10-SII | WG   |     2 | -0.00405807 |
| PCV10-SII | D    |     1 |    -0.87945 |
| PCV10-SII | DG   |     2 |     -1.4211 |
| PCV10-SII | DW   |     2 |    0.053572 |
| PCV10-SII | DWG  |     3 | -0.00508768 |

Factorial interaction coefficients for deaths_averted_per100k. Raw coefficients retain output units; ICER contrasts are descriptive ratio changes, without an additive component-benefit interpretation.

</div>

<div id="tab:interactionqalys_gained_per100k">

| Product   | Term | Order | Coefficient |
|:----------|:-----|------:|------------:|
| Product   | Term | Order | Coefficient |
| PCV10-GSK | G    |     1 |    -92.8138 |
| PCV10-GSK | W    |     1 |     -12.448 |
| PCV10-GSK | WG   |     2 |  -0.0796344 |
| PCV10-GSK | D    |     1 |    -34.7977 |
| PCV10-GSK | DG   |     2 |    -23.1209 |
| PCV10-GSK | DW   |     2 |     1.78751 |
| PCV10-GSK | DWG  |     3 |   -0.111064 |
| PCV13     | G    |     1 |    -115.374 |
| PCV13     | W    |     1 |    -13.4174 |
| PCV13     | WG   |     2 |  -0.0935618 |
| PCV13     | D    |     1 |     -37.182 |
| PCV13     | DG   |     2 |     -26.359 |
| PCV13     | DW   |     2 |     1.92585 |
| PCV13     | DWG  |     3 |   -0.123254 |
| PCV10-SII | G    |     1 |    -291.455 |
| PCV10-SII | W    |     1 |    -9.14028 |
| PCV10-SII | WG   |     2 |   -0.111971 |
| PCV10-SII | D    |     1 |    -22.3433 |
| PCV10-SII | DG   |     2 |    -39.5703 |
| PCV10-SII | DW   |     2 |     1.31626 |
| PCV10-SII | DWG  |     3 |   -0.138259 |

Factorial interaction coefficients for qalys_gained_per100k. Raw coefficients retain output units; ICER contrasts are descriptive ratio changes, without an additive component-benefit interpretation.

</div>

<div id="tab:interactionincremental_costs_per100k">

| Product   | Term | Order | Coefficient |
|:----------|:-----|------:|------------:|
| Product   | Term | Order | Coefficient |
| PCV10-GSK | G    |     1 |      303818 |
| PCV10-GSK | W    |     1 |     74667.5 |
| PCV10-GSK | WG   |     2 |     305.559 |
| PCV10-GSK | D    |     1 |    -18578.5 |
| PCV10-GSK | DG   |     2 |       77120 |
| PCV10-GSK | DW   |     2 |    -10930.6 |
| PCV10-GSK | DWG  |     3 |     476.622 |
| PCV13     | G    |     1 |      438980 |
| PCV13     | W    |     1 |     60454.7 |
| PCV13     | WG   |     2 |      360.24 |
| PCV13     | D    |     1 |     -104814 |
| PCV13     | DG   |     2 |     95163.3 |
| PCV13     | DW   |     2 |    -8818.32 |
| PCV13     | DWG  |     3 |     520.035 |
| PCV10-SII | G    |     1 | 1.01862e+06 |
| PCV10-SII | W    |     1 |     42461.8 |
| PCV10-SII | WG   |     2 |     409.039 |
| PCV10-SII | D    |     1 |    -69386.8 |
| PCV10-SII | DG   |     2 |      136907 |
| PCV10-SII | DW   |     2 |    -6215.27 |
| PCV10-SII | DWG  |     3 |     556.248 |

Factorial interaction coefficients for incremental_costs_per100k. Raw coefficients retain output units; ICER contrasts are descriptive ratio changes, without an additive component-benefit interpretation.

</div>

<div id="tab:interactionicer">

| Product   | Term | Order | Coefficient |
|:----------|:-----|------:|------------:|
| Product   | Term | Order | Coefficient |
| PCV10-GSK | G    |     1 |     31.3724 |
| PCV10-GSK | W    |     1 |     25.5138 |
| PCV10-GSK | WG   |     2 |     1.92026 |
| PCV10-GSK | D    |     1 |    -73.9427 |
| PCV10-GSK | DG   |     2 |     4.43024 |
| PCV10-GSK | DW   |     2 |    -3.84228 |
| PCV10-GSK | DWG  |     3 |     0.28391 |
| PCV13     | G    |     1 |     91.2146 |
| PCV13     | W    |     1 |     15.5035 |
| PCV13     | WG   |     2 |     1.98497 |
| PCV13     | D    |     1 |    -118.914 |
| PCV13     | DG   |     2 |     11.1529 |
| PCV13     | DW   |     2 |    -2.84967 |
| PCV13     | DWG  |     3 |    0.181106 |
| PCV10-SII | G    |     1 |     101.215 |
| PCV10-SII | W    |     1 |     8.08974 |
| PCV10-SII | WG   |     2 |       2.028 |
| PCV10-SII | D    |     1 |    -69.1744 |
| PCV10-SII | DG   |     2 |     3.09744 |
| PCV10-SII | DW   |     2 |    -1.42121 |
| PCV10-SII | DWG  |     3 |  0.00967258 |

Factorial interaction coefficients for icer. Raw coefficients retain output units; ICER contrasts are descriptive ratio changes, without an additive component-benefit interpretation.

</div>

<div id="tab:shapley">

| Product   | Factor |   Change | Joint share (%) |
|:----------|:-------|---------:|----------------:|
| Product   | Factor |   Change | Joint share (%) |
| PCV10-GSK | D      |  -45.501 |          28.160 |
| PCV10-GSK | W      |  -11.631 |           7.198 |
| PCV10-GSK | G      | -104.451 |          64.642 |
| PCV13     | D      |  -49.440 |          25.936 |
| PCV13     | W      |  -12.542 |           6.580 |
| PCV13     | G      | -128.641 |          67.485 |
| PCV10-SII | D      |  -41.516 |          11.486 |
| PCV10-SII | W      |   -8.584 |           2.375 |
| PCV10-SII | G      | -311.343 |          86.139 |

Shapley allocations for QALY gains, including interactions. Values are signed changes per 100,000 births.

</div>

<div id="tab:herd">

| Product / DWG | QALY off | QALY on | Indirect QALY | Indirect cases | Indirect deaths | Savings |
|:---|---:|---:|---:|---:|---:|---:|
| Product / DWG | QALY off | QALY on | Indirect QALY | Indirect cases | Indirect deaths | Savings |
| PCV10-GSK / 000 | 636.06 | 1588.60 | 952.54 | 7892.54 | 37.595 | 3.1387 |
| PCV10-GSK / 001 | 636.06 | 1495.79 | 859.73 | 7258.90 | 34.362 | 2.8349 |
| PCV10-GSK / 010 | 606.12 | 1576.15 | 970.03 | 8345.35 | 38.308 | 3.2346 |
| PCV10-GSK / 011 | 606.12 | 1483.26 | 877.14 | 7710.76 | 35.072 | 2.9305 |
| PCV10-GSK / 100 | 563.01 | 1553.80 | 990.79 | 8276.42 | 39.125 | 3.2683 |
| PCV10-GSK / 101 | 563.01 | 1437.87 | 874.85 | 7459.62 | 35.055 | 2.8874 |
| PCV10-GSK / 110 | 537.39 | 1543.14 | 1005.75 | 8662.03 | 39.736 | 3.3501 |
| PCV10-GSK / 111 | 537.39 | 1427.02 | 889.63 | 7842.58 | 35.658 | 2.9684 |
| PCV13 / 000 | 702.06 | 1713.29 | 1011.23 | 10181.47 | 39.741 | 3.6184 |
| PCV13 / 001 | 702.06 | 1597.91 | 895.85 | 9019.36 | 35.721 | 3.1794 |
| PCV13 / 010 | 669.05 | 1699.87 | 1030.83 | 10645.72 | 40.540 | 3.7196 |
| PCV13 / 011 | 669.05 | 1584.40 | 915.36 | 9482.46 | 36.518 | 3.2803 |
| PCV13 / 100 | 621.42 | 1676.11 | 1054.68 | 10580.06 | 41.477 | 3.7598 |
| PCV13 / 101 | 621.42 | 1534.37 | 912.95 | 9164.56 | 36.505 | 3.2256 |
| PCV13 / 110 | 593.16 | 1664.61 | 1071.45 | 10975.42 | 42.161 | 3.8462 |
| PCV13 / 111 | 593.16 | 1522.66 | 929.50 | 9556.94 | 37.181 | 3.3111 |
| PCV10-SII / 000 | 575.89 | 1966.40 | 1390.50 | 13348.65 | 53.485 | 4.8630 |
| PCV10-SII / 001 | 575.89 | 1674.94 | 1099.05 | 10797.99 | 43.450 | 3.8444 |
| PCV10-SII / 010 | 548.81 | 1957.25 | 1408.44 | 13755.32 | 54.214 | 4.9532 |
| PCV10-SII / 011 | 548.81 | 1665.69 | 1116.87 | 11203.42 | 44.175 | 3.9341 |
| PCV10-SII / 100 | 509.75 | 1944.05 | 1434.30 | 13712.73 | 55.210 | 5.0001 |
| PCV10-SII / 101 | 509.75 | 1613.03 | 1103.28 | 10810.10 | 43.753 | 3.8446 |
| PCV10-SII / 110 | 486.57 | 1936.23 | 1449.66 | 14059.14 | 55.834 | 5.0771 |
| PCV10-SII / 111 | 486.57 | 1604.95 | 1118.38 | 11153.37 | 44.368 | 3.9206 |

All matched herd-off/on contrasts. QALYs are gains versus NV; savings are million USD per 100,000 births.

</div>

## S8. All schedule outcomes

Each strategy uses W=G=1. Recovery combines original opportunities with 0.20 opportunities in the declared older age windows; it changes completion as well as timing. The model does not identify a clinical optimum.

<div id="tab:allschedules">

| Product / strategy             |    Cases | Deaths |    QALY | Net cost |    Doses |
|:-------------------------------|---------:|-------:|--------:|---------:|---------:|
| Product / strategy             |    Cases | Deaths |    QALY | Net cost |    Doses |
| PCV10-GSK / 2-4-12             | 13602.81 | 58.969 | 1483.26 |  -4.0182 | 253865.5 |
| PCV10-GSK / 3-5-12             | 13503.57 | 58.261 | 1464.45 |  -3.9731 | 252281.9 |
| PCV10-GSK / 2-5-12             | 13549.05 | 58.585 | 1473.08 |  -3.9937 | 253034.2 |
| PCV10-GSK / 3-6-12             | 13450.76 | 57.884 | 1454.48 |  -3.9492 | 251433.1 |
| PCV10-GSK / 2-4-15             | 13548.23 | 58.728 | 1476.71 |  -4.0051 | 251731.7 |
| PCV10-GSK / 3-5-15             | 13448.66 | 58.017 | 1457.84 |  -3.9599 | 250148.2 |
| PCV10-GSK / delayed            | 13131.33 | 56.854 | 1427.02 |  -3.9701 | 217271.3 |
| PCV10-GSK / dose1_only         |  7299.21 | 31.121 |  779.54 |  -2.2898 |  87683.1 |
| PCV10-GSK / completion_catchup | 13380.33 | 57.765 | 1451.12 |  -3.9585 | 243849.1 |
| PCV13 / 2-4-12                 | 14665.24 | 62.889 | 1584.40 |  -3.8546 | 253867.9 |
| PCV13 / 3-5-12                 | 14589.36 | 62.125 | 1564.09 |  -3.8142 | 252284.6 |
| PCV13 / 2-5-12                 | 14624.22 | 62.475 | 1573.42 |  -3.8326 | 253036.7 |
| PCV13 / 3-6-12                 | 14549.17 | 61.719 | 1553.33 |  -3.7929 | 251435.8 |
| PCV13 / 2-4-15                 | 14600.32 | 62.614 | 1576.93 |  -3.8428 | 251734.3 |
| PCV13 / 3-5-15                 | 14524.02 | 61.847 | 1556.54 |  -3.8022 | 250150.9 |
| PCV13 / delayed                | 14208.26 | 60.572 | 1522.66 |  -3.8725 | 217273.3 |
| PCV13 / dose1_only             |  7785.69 | 33.334 |  836.09 |  -2.2910 |  87683.2 |
| PCV13 / completion_catchup     | 14448.17 | 61.565 | 1548.98 |  -3.8104 | 243851.7 |
| PCV10-SII / 2-4-12             | 15455.10 | 65.807 | 1665.69 |  -4.6672 | 253870.4 |
| PCV10-SII / 3-5-12             | 15401.47 | 65.289 | 1651.94 |  -4.6390 | 252287.8 |
| PCV10-SII / 2-5-12             | 15426.19 | 65.527 | 1658.28 |  -4.6519 | 253039.4 |
| PCV10-SII / 3-6-12             | 15373.23 | 65.015 | 1644.70 |  -4.6243 | 251439.3 |
| PCV10-SII / 2-4-15             | 15373.95 | 65.475 | 1656.53 |  -4.6444 | 251736.7 |
| PCV10-SII / 3-5-15             | 15319.85 | 64.953 | 1642.67 |  -4.6161 | 250154.1 |
| PCV10-SII / delayed            | 14969.07 | 63.555 | 1604.95 |  -4.6053 | 217275.4 |
| PCV10-SII / dose1_only         |  8009.97 | 34.293 |  863.80 |  -2.5657 |  87683.4 |
| PCV10-SII / completion_catchup | 15224.21 | 64.587 | 1632.69 |  -4.6034 | 243854.5 |

Complete schedule outcomes per 100,000 births; net costs in million USD.

</div>

## S9. Complete deterministic sensitivity

Changes are relative to product-specific full111. Cost percentages use the absolute baseline denominator; a positive change in a negative net cost means less saving. These are screening values, not PSA or empirical uncertainty bounds.

<div id="tab:sensitivityPCV10-GSK">

| Setting                | Value |   Cases | Deaths |    QALY | Net cost |  ICER |
|:-----------------------|------:|--------:|-------:|--------:|---------:|------:|
| Setting                | Value |   Cases | Deaths |    QALY | Net cost |  ICER |
| price                  |  0.75 | 13131.3 |  56.85 | 1427.02 |   -4.066 | -2849 |
| price                  |  1.25 | 13131.3 |  56.85 | 1427.02 |   -3.874 | -2715 |
| direct_VE              |  0.75 | 12485.1 |  54.25 | 1359.67 |   -3.746 | -2755 |
| direct_VE              |  1.25 | 13777.5 |  59.45 | 1494.36 |   -4.194 | -2806 |
| CFR                    |  0.75 | 13133.8 |  42.65 | 1071.76 |   -3.441 | -3211 |
| CFR                    |  1.25 | 13128.8 |  71.06 | 1782.17 |   -4.499 | -2524 |
| sequela_cost           |  0.75 | 13131.3 |  56.85 | 1427.02 |   -3.968 | -2781 |
| sequela_cost           |  1.25 | 13131.3 |  56.85 | 1427.02 |   -3.972 | -2783 |
| utility_loss           |  0.75 | 13131.3 |  56.85 | 1425.71 |   -3.970 | -2785 |
| utility_loss           |  1.25 | 13131.3 |  56.85 | 1428.32 |   -3.970 | -2780 |
| burden                 |  0.75 |  9853.4 |  42.65 | 1070.60 |   -2.778 | -2595 |
| burden                 |  1.25 | 16405.9 |  71.04 | 1783.21 |   -5.161 | -2894 |
| waning                 |     0 | 13451.0 |  57.29 | 1437.87 |   -4.035 | -2806 |
| waning                 |  1.25 | 13051.4 |  56.74 | 1424.30 |   -3.954 | -2776 |
| herd_strength          |   0.5 |  9210.7 |  39.02 |  982.16 |   -2.486 | -2531 |
| herd_strength          |     0 |  5288.7 |  21.20 |  537.39 |   -1.002 | -1864 |
| delay_probability      |  0.75 | 12210.6 |  52.81 | 1321.70 |   -3.768 | -2851 |
| delay_probability      |  1.25 | 13464.3 |  58.33 | 1466.32 |   -4.012 | -2736 |
| birth_seasonality      |     0 | 13134.0 |  56.87 | 1427.41 |   -3.972 | -2783 |
| healthcare_perspective |     0 | 13131.3 |  56.85 | 1427.02 |    0.133 |    94 |

Deterministic settings for PCV10-GSK. Values per 100,000 births; costs in million USD.

</div>

<div id="tab:sensitivityPCV13">

| Setting                | Value |   Cases | Deaths |    QALY | Net cost |  ICER |
|:-----------------------|------:|--------:|-------:|--------:|---------:|------:|
| Setting                | Value |   Cases | Deaths |    QALY | Net cost |  ICER |
| price                  |  0.75 | 14208.3 |  60.57 | 1522.66 |   -4.080 | -2680 |
| price                  |  1.25 | 14208.3 |  60.57 | 1522.66 |   -3.665 | -2407 |
| direct_VE              |  0.75 | 13740.1 |  57.77 | 1450.00 |   -3.665 | -2528 |
| direct_VE              |  1.25 | 14676.4 |  63.38 | 1595.33 |   -4.080 | -2557 |
| CFR                    |  0.75 | 14210.9 |  45.44 | 1143.61 |   -3.308 | -2893 |
| CFR                    |  1.25 | 14205.7 |  75.70 | 1901.61 |   -4.437 | -2333 |
| sequela_cost           |  0.75 | 14208.3 |  60.57 | 1522.66 |   -3.871 | -2542 |
| sequela_cost           |  1.25 | 14208.3 |  60.57 | 1522.66 |   -3.874 | -2544 |
| utility_loss           |  0.75 | 14208.3 |  60.57 | 1521.26 |   -3.873 | -2546 |
| utility_loss           |  1.25 | 14208.3 |  60.57 | 1524.07 |   -3.873 | -2541 |
| burden                 |  0.75 | 10661.6 |  45.44 | 1142.35 |   -2.593 | -2270 |
| burden                 |  1.25 | 17751.3 |  75.69 | 1902.75 |   -5.151 | -2707 |
| waning                 |     0 | 14428.8 |  61.04 | 1534.37 |   -3.925 | -2558 |
| waning                 |  1.25 | 14153.1 |  60.45 | 1519.74 |   -3.859 | -2540 |
| herd_strength          |   0.5 |  9430.7 |  41.98 | 1057.87 |   -2.217 | -2096 |
| herd_strength          |     0 |  4651.3 |  23.39 |  593.16 |   -0.561 |  -946 |
| delay_probability      |  0.75 | 13251.3 |  56.24 | 1409.56 |   -3.731 | -2647 |
| delay_probability      |  1.25 | 14532.6 |  62.19 | 1565.76 |   -3.873 | -2473 |
| birth_seasonality      |     0 | 14211.3 |  60.59 | 1523.10 |   -3.875 | -2544 |
| healthcare_perspective |     0 | 14208.3 |  60.57 | 1522.66 |    0.535 |   351 |

Deterministic settings for PCV13. Values per 100,000 births; costs in million USD.

</div>

<div id="tab:sensitivityPCV10-SII">

| Setting                | Value |   Cases | Deaths |    QALY | Net cost |  ICER |
|:-----------------------|------:|--------:|-------:|--------:|---------:|------:|
| Setting                | Value |   Cases | Deaths |    QALY | Net cost |  ICER |
| price                  |  0.75 | 14969.1 |  63.56 | 1604.95 |   -4.701 | -2929 |
| price                  |  1.25 | 14969.1 |  63.56 | 1604.95 |   -4.510 | -2810 |
| direct_VE              |  0.75 | 14636.8 |  61.65 | 1555.73 |   -4.463 | -2868 |
| direct_VE              |  1.25 | 15301.3 |  65.46 | 1654.17 |   -4.748 | -2870 |
| CFR                    |  0.75 | 14971.8 |  47.67 | 1205.45 |   -4.011 | -3327 |
| CFR                    |  1.25 | 14966.4 |  79.43 | 2004.34 |   -5.200 | -2594 |
| sequela_cost           |  0.75 | 14969.1 |  63.56 | 1604.95 |   -4.603 | -2868 |
| sequela_cost           |  1.25 | 14969.1 |  63.56 | 1604.95 |   -4.607 | -2871 |
| utility_loss           |  0.75 | 14969.1 |  63.56 | 1603.42 |   -4.605 | -2872 |
| utility_loss           |  1.25 | 14969.1 |  63.56 | 1606.48 |   -4.605 | -2867 |
| burden                 |  0.75 | 11232.6 |  47.68 | 1204.08 |   -3.255 | -2703 |
| burden                 |  1.25 | 18701.7 |  79.42 | 2005.58 |   -5.955 | -2969 |
| waning                 |     0 | 15128.7 |  63.88 | 1613.03 |   -4.643 | -2878 |
| waning                 |  1.25 | 14929.2 |  63.47 | 1602.93 |   -4.596 | -2867 |
| herd_strength          |   0.5 |  9393.5 |  41.37 | 1045.70 |   -2.645 | -2530 |
| herd_strength          |     0 |  3815.7 |  19.19 |  486.57 |   -0.685 | -1407 |
| delay_probability      |  0.75 | 13924.9 |  59.02 | 1485.78 |   -4.354 | -2931 |
| delay_probability      |  1.25 | 15316.1 |  65.15 | 1647.84 |   -4.658 | -2827 |
| birth_seasonality      |     0 | 14974.3 |  63.58 | 1605.63 |   -4.608 | -2870 |
| healthcare_perspective |     0 | 14969.1 |  63.56 | 1604.95 |    0.048 |    30 |

Deterministic settings for PCV10-SII. Values per 100,000 births; costs in million USD.

</div>

## S10. AOM and mortality evidence-sensitive analyses

The primary AOM study recruited 482 children, with 466 spontaneous perforations and 16 complicated presentations; this is a selected clinical population. Culture and additional PCR yield distinct detection proportions, neither a national exclusive causal fraction. The frozen proxy is retained as a scenario constant for reproducible comparisons. The national incidence input also involves transfer from a Philippine assessment; age coverage and care-seeking should be independently calibrated. A case with $`\theta=1`$ applies PCV effects to the entire baseline as a mathematical bound, not an all-cause unaffected-component reconstruction.

<div id="tab:aom">

| Product        | Fraction |  AOM events | QALY gain | Net cost (million) |
|:---------------|---------:|------------:|----------:|-------------------:|
| Product        | Fraction |  AOM events | QALY gain | Net cost (million) |
| No vaccination |    0.303 | 1,971,313.6 |      0.00 |             0.0000 |
| PCV10-GSK      |    0.303 |   903,528.2 |   1427.02 |            -3.9701 |
| PCV13          |    0.303 |   803,725.8 |   1522.66 |            -3.8725 |
| PCV10-SII      |    0.303 |   732,338.2 |   1604.95 |            -4.6053 |
| No vaccination |    0.301 | 1,958,286.0 |      0.00 |             0.0000 |
| PCV10-GSK      |    0.301 |   897,559.8 |   1427.01 |            -3.9631 |
| PCV13          |    0.301 |   798,417.0 |   1522.66 |            -3.8649 |
| PCV10-SII      |    0.301 |   727,501.2 |   1604.95 |            -4.5971 |
| No vaccination |    0.411 | 2,675,112.7 |      0.00 |             0.0000 |
| PCV10-GSK      |    0.411 | 1,225,909.5 |   1427.24 |            -4.3480 |
| PCV13          |    0.411 | 1,090,475.6 |   1522.91 |            -4.2871 |
| PCV10-SII      |    0.411 |   993,595.4 |   1605.21 |            -5.0478 |
| No vaccination |    0.000 |         0.0 |      0.00 |             0.0000 |
| PCV10-GSK      |    0.000 |         0.0 |   1426.40 |            -2.9121 |
| PCV13          |    0.000 |         0.0 |   1521.99 |            -2.7120 |
| PCV10-SII      |    0.000 |         0.0 |   1604.23 |            -3.3664 |
| No vaccination |    1.000 | 6,524,141.2 |      0.00 |             0.0000 |
| PCV10-GSK      |    1.000 | 2,987,161.5 |   1428.45 |            -6.4166 |
| PCV13          |    1.000 | 2,656,862.6 |   1524.23 |            -6.5560 |
| PCV10-SII      |    1.000 | 2,420,510.7 |   1606.63 |            -7.4697 |

AOM attribution settings. Endpoint counts are total ten-year AOM events; gains/net costs are per 100,000 births.

</div>

The revised baseline subtracts the frozen no-vaccination disease hazard from the constant residual post-neonatal all-cause hazard. This prevents the specified mortality overlap while retaining explicit disease hazards. Alternative rows are historical inherited, neonatal-residual and zero-background accounting scenarios; their age allocation is not an empirical cause-deleted life table.

<div id="tab:mortality">

| Product        | Background                      | QALY gain | Net cost (million) |
|:---------------|:--------------------------------|----------:|-------------------:|
| Product        | Background                      | QALY gain | Net cost (million) |
| No vaccination | cause_deleted                   |      0.00 |             0.0000 |
| PCV10-GSK      | cause_deleted                   |   1427.02 |            -3.9701 |
| PCV13          | cause_deleted                   |   1522.66 |            -3.8725 |
| PCV10-SII      | cause_deleted                   |   1604.95 |            -4.6053 |
| No vaccination | historical_inherited            |      0.00 |             0.0000 |
| PCV10-GSK      | historical_inherited            |   1422.56 |            -3.9537 |
| PCV13          | historical_inherited            |   1517.94 |            -3.8555 |
| PCV10-SII      | historical_inherited            |   1600.14 |            -4.5871 |
| No vaccination | historical_neonatal_conditioned |      0.00 |             0.0000 |
| PCV10-GSK      | historical_neonatal_conditioned |   1426.32 |            -3.9677 |
| PCV13          | historical_neonatal_conditioned |   1521.92 |            -3.8701 |
| PCV10-SII      | historical_neonatal_conditioned |   1604.19 |            -4.6026 |
| No vaccination | zero_background                 |      0.00 |             0.0000 |
| PCV10-GSK      | zero_background                 |   1429.91 |            -3.9811 |
| PCV13          | zero_background                 |   1525.73 |            -3.8841 |
| PCV10-SII      | zero_background                 |   1608.06 |            -4.6175 |

Mortality accounting alternatives. Background is fixed by the designated convention; gains and costs use matched NV comparators.

</div>

## S11. Trajectories and late-period summaries

<div id="tab:stable">

| Product   | Endpoint          | NV events | Vaccine events | Reduction (%) |
|:----------|:------------------|----------:|---------------:|--------------:|
| Product   | Endpoint          | NV events | Vaccine events | Reduction (%) |
| PCV10-GSK | New_Meningitis    |     820.3 |          237.6 |         71.03 |
| PCV10-GSK | New_Pneumo_Sev    |  47,073.1 |       12,560.3 |         73.32 |
| PCV10-GSK | New_Pneumo_NonSev |  73,937.4 |       19,728.4 |         73.32 |
| PCV10-GSK | New_Sepsis        |   8,472.8 |        2,387.3 |         71.82 |
| PCV10-GSK | New_AOM           | 197,131.4 |       55,667.6 |         71.76 |
| PCV10-GSK | New_Deaths_Pneumo |   1,373.8 |          377.5 |         72.52 |
| PCV13     | New_Meningitis    |     820.3 |          208.4 |         74.59 |
| PCV13     | New_Pneumo_Sev    |  47,073.1 |       10,549.3 |         77.59 |
| PCV13     | New_Pneumo_NonSev |  73,937.4 |       16,569.7 |         77.59 |
| PCV13     | New_Sepsis        |   8,472.8 |        2,074.0 |         75.52 |
| PCV13     | New_AOM           | 197,131.4 |       46,175.3 |         76.58 |
| PCV13     | New_Deaths_Pneumo |   1,373.8 |          322.8 |         76.50 |
| PCV10-SII | New_Meningitis    |     820.3 |          201.7 |         75.42 |
| PCV10-SII | New_Pneumo_Sev    |  47,073.1 |       10,587.9 |         77.51 |
| PCV10-SII | New_Pneumo_NonSev |  73,937.4 |       16,630.4 |         77.51 |
| PCV10-SII | New_Sepsis        |   8,472.8 |        2,027.0 |         76.08 |
| PCV10-SII | New_AOM           | 197,131.4 |       45,789.4 |         76.77 |
| PCV10-SII | New_Deaths_Pneumo |   1,373.8 |          319.4 |         76.75 |

Year-10 endpoint counts and reductions; no fitted transmission equilibrium is claimed.

</div>

<div id="tab:trajectory">

| Scenario                |       C |       r | H (IPD) |
|:------------------------|--------:|--------:|--------:|
| Scenario                |       C |       r | H (IPD) |
| pcv10_gsk_D0W0G0        | 0.34912 | 1.00000 | 0.43046 |
| pcv10_gsk_D0W0G1        | 0.34912 | 1.00000 | 0.43046 |
| pcv10_gsk_D0W1G0        | 0.34912 | 1.00000 | 0.43046 |
| pcv10_gsk_D0W1G1        | 0.34912 | 1.00000 | 0.43046 |
| pcv10_gsk_D1W0G0        | 0.31052 | 1.00000 | 0.43046 |
| pcv10_gsk_D1W0G1        | 0.31052 | 1.00000 | 0.43046 |
| pcv10_gsk_D1W1G0        | 0.31052 | 1.00000 | 0.43046 |
| pcv10_gsk_D1W1G1        | 0.31052 | 1.00000 | 0.43046 |
| pcv10_gsk_D0W0_herd_off | 0.34912 | 1.00000 | 1.00000 |
| pcv10_gsk_D0W1_herd_off | 0.34912 | 1.00000 | 1.00000 |
| pcv10_gsk_D1W0_herd_off | 0.31053 | 1.00000 | 1.00000 |
| pcv10_gsk_D1W1_herd_off | 0.31053 | 1.00000 | 1.00000 |
| pcv13_D0W0G0            | 0.34912 | 1.00000 | 0.41914 |
| pcv13_D0W0G1            | 0.34912 | 1.00000 | 0.41914 |
| pcv13_D0W1G0            | 0.34912 | 1.00000 | 0.41914 |
| pcv13_D0W1G1            | 0.34912 | 1.00000 | 0.41914 |
| pcv13_D1W0G0            | 0.31052 | 1.00000 | 0.41914 |
| pcv13_D1W0G1            | 0.31052 | 1.00000 | 0.41914 |
| pcv13_D1W1G0            | 0.31052 | 1.00000 | 0.41914 |
| pcv13_D1W1G1            | 0.31052 | 1.00000 | 0.41914 |
| pcv13_D0W0_herd_off     | 0.34912 | 1.00000 | 1.00000 |
| pcv13_D0W1_herd_off     | 0.34912 | 1.00000 | 1.00000 |
| pcv13_D1W0_herd_off     | 0.31053 | 1.00000 | 1.00000 |
| pcv13_D1W1_herd_off     | 0.31053 | 1.00000 | 1.00000 |
| pcv10sii_D0W0G0         | 0.34912 | 1.00000 | 0.36312 |
| pcv10sii_D0W0G1         | 0.34912 | 1.00000 | 0.36312 |
| pcv10sii_D0W1G0         | 0.34912 | 1.00000 | 0.36312 |
| pcv10sii_D0W1G1         | 0.34912 | 1.00000 | 0.36312 |
| pcv10sii_D1W0G0         | 0.31052 | 1.00000 | 0.36312 |
| pcv10sii_D1W0G1         | 0.31052 | 1.00000 | 0.36312 |
| pcv10sii_D1W1G0         | 0.31052 | 1.00000 | 0.36312 |
| pcv10sii_D1W1G1         | 0.31052 | 1.00000 | 0.36312 |
| pcv10sii_D0W0_herd_off  | 0.34912 | 1.00000 | 1.00000 |
| pcv10sii_D0W1_herd_off  | 0.34912 | 1.00000 | 1.00000 |
| pcv10sii_D1W0_herd_off  | 0.31053 | 1.00000 | 1.00000 |
| pcv10sii_D1W1_herd_off  | 0.31053 | 1.00000 | 1.00000 |

Final-month proxy and feedback for core vaccinated scenarios.

</div>

**Extended trajectory figure S1.** [Open vector figure](../figures/trajectories/pcv10_gsk_population_trajectory.pdf). Product-specific no-vaccination, on-time and delayed programme trajectories.

**Extended trajectory figure S2.** [Open vector figure](../figures/trajectories/pcv13_population_trajectory.pdf). Product-specific no-vaccination, on-time and delayed programme trajectories.

**Extended trajectory figure S3.** [Open vector figure](../figures/trajectories/pcv10_sii_population_trajectory.pdf). Product-specific no-vaccination, on-time and delayed programme trajectories.

**Extended trajectory figure S4.** [Open vector figure](../figures/trajectories/pcv10_gsk_dose_integrals.pdf). Product-specific dose-state integrals for on-time and delayed delivery.

**Extended trajectory figure S5.** [Open vector figure](../figures/trajectories/pcv13_dose_integrals.pdf). Product-specific dose-state integrals for on-time and delayed delivery.

**Extended trajectory figure S6.** [Open vector figure](../figures/trajectories/pcv10_sii_dose_integrals.pdf). Product-specific dose-state integrals for on-time and delayed delivery.

**Extended trajectory figure S7.** [Open vector figure](../figures/trajectories/pcv10_gsk_age_calendar_dose_density.pdf). State-resolved age–calendar densities for no vaccination, on-time and delayed delivery.

**Extended trajectory figure S8.** [Open vector figure](../figures/trajectories/pcv13_age_calendar_dose_density.pdf). State-resolved age–calendar densities for no vaccination, on-time and delayed delivery.

**Extended trajectory figure S9.** [Open vector figure](../figures/trajectories/pcv10_sii_age_calendar_dose_density.pdf). State-resolved age–calendar densities for no vaccination, on-time and delayed delivery.

## S12. Full numerical ledger and reproducibility

<div id="tab:allrefinement">

| Product   | h (month) |    QALY gain | Gain minus reference | Difference (%) |
|:----------|----------:|-------------:|---------------------:|---------------:|
| Product   | h (month) |    QALY gain | Gain minus reference | Difference (%) |
| PCV10-GSK |         1 | 1,427.015296 |             5.482497 |       0.385675 |
| PCV10-GSK |       1/2 | 1,424.093420 |             2.560621 |       0.180131 |
| PCV10-GSK |       1/4 | 1,422.631457 |             1.098658 |       0.077287 |
| PCV10-GSK |       1/8 | 1,421.899061 |             0.366262 |       0.025765 |
| PCV10-GSK |      1/16 | 1,421.532799 |             0.000000 |       0.000000 |
| PCV13     |         1 | 1,522.664903 |             6.319425 |       0.416754 |
| PCV13     |       1/2 | 1,519.296993 |             2.951516 |       0.194647 |
| PCV13     |       1/4 | 1,517.611727 |             1.266250 |       0.083507 |
| PCV13     |       1/8 | 1,516.767616 |             0.422138 |       0.027839 |
| PCV13     |      1/16 | 1,516.345478 |             0.000000 |       0.000000 |
| PCV10-SII |         1 | 1,604.952058 |             8.456476 |       0.529690 |
| PCV10-SII |       1/2 | 1,600.446399 |             3.950817 |       0.247468 |
| PCV10-SII |       1/4 | 1,598.190696 |             1.695114 |       0.106177 |
| PCV10-SII |       1/8 | 1,597.060740 |             0.565158 |       0.035400 |
| PCV10-SII |      1/16 | 1,596.495582 |             0.000000 |       0.000000 |

QALY-gain refinement ledger for full111, per 100,000 recorded-window births. Differences are relative to the finest computed grid, $`h=1/16`$ month; zero at that grid is definitional and does not assert zero numerical error.

</div>

<div id="tab:manufactured">

|   k |       h | L1 mass error | Population error | Observed order |
|----:|--------:|--------------:|-----------------:|---------------:|
|   k |       h | L1 mass error | Population error | Observed order |
|   1 | 1.00000 |    234.829319 |        22.036801 |              – |
|   2 | 0.50000 |    121.523356 |        10.645246 |         1.0497 |
|   4 | 0.25000 |     61.230118 |         5.231310 |         1.0250 |
|   8 | 0.12500 |     30.685289 |         2.592947 |         1.0126 |
|  16 | 0.06250 |     15.357069 |         1.290804 |         1.0063 |
|  32 | 0.03125 |      7.681285 |         0.643985 |         1.0032 |

Manufactured transport test. Errors are evaluated against exact cell integrals and exact total population.

</div>

<div id="tab:ctmc">

| Time | State | Exact probability | ODE probability |    MC mean |      MC SE |
|-----:|------:|------------------:|----------------:|-----------:|-----------:|
| Time | State | Exact probability | ODE probability |    MC mean |      MC SE |
| 0.25 |     0 |        0.94885432 |      0.94885432 | 0.94882796 | 0.00002293 |
| 0.25 |     1 |        0.04772851 |      0.04772851 | 0.04775583 | 0.00002229 |
| 0.25 |     2 |        0.00090151 |      0.00090151 | 0.00090228 | 0.00000299 |
| 0.25 |     3 |        0.00000606 |      0.00000606 | 0.00000643 | 0.00000025 |
| 0.25 |     4 |        0.00250960 |      0.00250960 | 0.00250750 | 0.00000511 |
|  1.0 |     0 |        0.81058425 |      0.81058425 | 0.81056713 | 0.00003946 |
|  1.0 |     1 |        0.16607066 |      0.16607066 | 0.16611760 | 0.00003707 |
|  1.0 |     2 |        0.01282821 |      0.01282821 | 0.01281158 | 0.00001121 |
|  1.0 |     3 |        0.00035388 |      0.00035388 | 0.00035366 | 0.00000193 |
|  1.0 |     4 |        0.01016300 |      0.01016300 | 0.01015003 | 0.00000995 |
|  3.0 |     0 |        0.53259180 |      0.53259180 | 0.53254137 | 0.00004900 |
|  3.0 |     1 |        0.34370836 |      0.34370836 | 0.34373853 | 0.00004717 |
|  3.0 |     2 |        0.08457593 |      0.08457593 | 0.08461266 | 0.00002796 |
|  3.0 |     3 |        0.00750111 |      0.00750111 | 0.00748953 | 0.00000879 |
|  3.0 |     4 |        0.03162280 |      0.03162280 | 0.03161791 | 0.00001698 |
| 12.0 |     0 |        0.08045961 |      0.08045961 | 0.08038851 | 0.00002767 |
| 12.0 |     1 |        0.26112781 |      0.26112781 | 0.26114936 | 0.00004262 |
| 12.0 |     2 |        0.34243499 |      0.34243499 | 0.34246058 | 0.00004714 |
| 12.0 |     3 |        0.16992181 |      0.16992181 | 0.16992462 | 0.00003791 |
| 12.0 |     4 |        0.14605579 |      0.14605579 | 0.14607693 | 0.00003607 |

Continuous-time chain probability verification.

</div>

The complete application rerun passed 750 runtime and 7,544 independent checks. The method experiment script passed 19,283 checks, predominantly per-step population balances; check counts are not independent scientific evidence. The retained cohort audit uses 4,000 replicates in total; this revision rechecks its stored replicate-level summaries. The syntax audit parses retained Python sources without executing legacy scripts that can write original outputs. Source hashes, parameter records, plot tables and reference metadata responses are retained in the separate code/data archive. The historical cohort PSA summary reported 200 draws; its original draws and seed are unavailable, so it is not evidence for current calendar-model PSA. At the time of the original 2026-10-27 record, the computational package was local. That computational audit did not establish permission to redistribute third-party unpublished aggregates or substitute for public deposition. For the current public release and depositor permissions statement, see [Data and code availability](availability.md) and [Licensing and source rights](../LICENSES.md). Original-source citations support identified published claims; exact constants without primary numerical confirmation remain flagged. Author and institutional declarations must be confirmed before submission.

## S13. Additional verification and sensitivity evidence from model development

The development record contains several distinct estimands. We preserve them here as verification or exploratory evidence, rather than pooling their outputs with the present calendar application. New reconstructions use frozen copies of the original code and inputs. The accompanying source manifest records original paths, SHA-256 hashes and retained files. The reconstructions do not replace the present 000, 011 or 111 application ledger.

#### Exact propagation, accounting checks and local sensitivity.

For the early two-subcohort model with coefficients fixed inside a month, alive mass obeys $`A_g(m+\tau)=A_g(m)\exp(-\mu_{gm}\tau)`$, and a cumulative reward with source rate $`\alpha_{gm}`$ receives $`\alpha_{gm}A_g(m)\phi(\mu_{gm},\tau)`$, where $`\phi(\mu,\tau)=(1-\exp(-\mu\tau))/\mu`$ and $`\phi(0,\tau)=\tau`$. This is the origin of the frozen-coefficient reaction step used in the present transport algorithm. The earlier sequela-stock rule also has an integrating-factor representation conditional on its specified survival-average approximation.

For each of meningitis, pneumonia, sepsis and AOM, the 60 monthly age weights sum to one to within $`10^{-12}`$. Recomputing the PCV13 dynamic cohort over 12, 24, 36, 48 and 60 months leaves every output in the common prefix unchanged, with zero recorded difference. This tests horizon bookkeeping; it does not validate the numerical values of the weights. In six scenarios and 120 reported scenario–metric comparisons, exact propagation versus the original adaptive RK45 implementation gives a largest absolute difference of 1.57344 events in a cumulative nonsevere-pneumonia count and a largest scaled difference of $`9.98598\times10^{-4}`$. The scaling uses $`\max(|Y_{\rm RK45}|,1)`$ in the denominator. Exact propagation and the original adaptive solver therefore agree closely at the reported scale, but not at machine precision. The characteristic rewrite below uses the same adaptive numerical convention as the original ODE and addresses a different equivalence check.

Early local elasticity is calculated as
``` math
E_{Y,\theta}=\frac{\log|Y(\theta_{+})|-\log|Y(\theta_{-})|}{\log\theta_{+}-\log\theta_{-}},
```
with a $`10^{-12}`$ numerical floor for zero endpoint magnitudes. Perturbations are 90% and 110% of each retained parameter, except programme coverage (95% and 105% of its retained value). These are deterministic screening choices. Table <a href="#tab:legacy_elasticity" data-reference-type="ref" data-reference="tab:legacy_elasticity">35</a> reports the reconstructed PCV13 dynamic-cohort elasticities. They concern the early IRR-to-VE mapping, not the present three-product nonlocal feedback. In particular, coverage is almost offset by that mapping, and the resulting near-zero elasticity must not be interpreted as a biological lack of coverage effect.

<div id="tab:legacy_elasticity">

| Parameter         | Cases averted | Deaths averted |  QALY gain | Cost magnitude |
|:------------------|--------------:|---------------:|-----------:|---------------:|
| inc_pneumonia     |       0.33688 |        0.48575 |    0.48059 |        0.68593 |
| inc_meningitis    |     0.0021499 |       0.018222 |   0.028195 |       0.027791 |
| inc_sepsis        |      0.023176 |         0.4922 |    0.49429 |        0.40641 |
| inc_aom_all       |       0.63319 |              0 |  0.0005099 |        0.38294 |
| cfr_meningitis    |   -1.7195e-05 |       0.018222 |   0.025249 |       0.016445 |
| cfr_pneumo_severe |   -0.00040501 |         0.4863 |     0.4803 |        0.30437 |
| cfr_sepsis        |   -0.00041007 |        0.49224 |    0.49424 |        0.31261 |
| vaccine_coverage  |    0.00039285 |    -0.00030855 |  2.463e-07 |     0.00019165 |
| prob_hear         |             0 |              0 | 0.00032192 |      0.0014947 |
| prob_hemi         |             0 |              0 |  0.0014032 |     0.00027528 |

Reconstructed early PCV13 cohort local elasticities. Cost elasticity uses the magnitude of the incremental cost; a positive value does not denote increasing expenditure. The early societal ledger is retained only for this separate benchmark.

</div>

Structural screening preserves the analytic propagation while changing the vaccine-effect map. Direct-effect-only mapping leaves nonvaccinated children at baseline risk; uniform-population-IRR applies the same annual IRR to both subcohorts; partial-herd screens additionally reduce nonvaccinated risk by 25% or 50% of the annual $`1-\mathrm{IRR}`$ reduction. The four reconstructed outcomes are retained in Table <a href="#tab:legacy_structural" data-reference-type="ref" data-reference="tab:legacy_structural">36</a>. The population is one cohort of 1,384,610 births followed for five years, rather than ten calendar years of successive cohorts. The relatively small uniform-IRR change and larger partial-herd changes motivated explicit programme feedback in the later model; they are not an estimate of the current indirect contribution.

<div id="tab:legacy_structural">

| Mapping                | Cases averted | Deaths averted | QALY gain | Cost (million) |
|:-----------------------|--------------:|---------------:|----------:|---------------:|
| direct_effect_only     |    199,647.41 |         765.28 | 21,902.09 |       -49.9173 |
| uniform_population_irr |    199,655.24 |         765.26 | 21,902.10 |       -49.9183 |
| partial_herd_25pct     |    204,640.12 |         784.41 | 22,446.44 |       -51.7878 |
| partial_herd_50pct     |    209,631.94 |         803.54 | 22,962.20 |       -53.6137 |

Early cohort vaccine-effect structural screens reconstructed from frozen code. Cost is incremental societal cost in millions of 2023-aligned USD under the retained price-conversion proxy. These entries are historical-model illustrations, not present calendar estimates.

</div>

#### PCV-specific ODE–characteristic and independent transport checks.

The original PCV characteristic comparison initializes one age-zero cohort, applies the same neonatal adjustment, introduces no subsequent births and uses no warm-up. Age follows $`a=t`$ for 60 months. The ODE and characteristic rewrite share inputs, endpoint rules, reporting boundaries and adaptive solver settings. New reconstruction of baseline, dynamic GSK, dynamic PCV13 and three mature-effect immediate scenarios gives zero recorded difference in all 120 scenario–metric monthly comparisons, and zero difference in all 120 corresponding rounded annual comparisons. This is a controlled rewrite check under matched conditions. It does not establish external validation or convergence of the full nonlocal calendar model.

A separately written simplified monthly scalar-versus-age-bin transport check compares baseline and mature PCV13 using six quantities: alive population, meningitis, pneumonia, sepsis, AOM and deaths. All 12 maximum-difference checks are below $`10^{-8}`$; the largest is $`1.3969838619\times10^{-9}`$ in alive population. This simplified check applies explicit monthly event accounting without dose transitions, continuing births or gradual herd feedback. Its error cannot be assigned to the original six-scenario characteristic comparison or to the current nonlinear solver.

#### Retained closed-cohort sensitivity benchmark.

The independently specified mature-effect PCV13 closed cohort was previously screened in 17 settings. Reconstructing the saved table reproduces all five outcome columns exactly at the stored precision. Table <a href="#tab:legacy_cohort_dsa" data-reference-type="ref" data-reference="tab:legacy_cohort_dsa">37</a> retains these results as a single-cohort, 2024-dollar benchmark. The total-population IRR is applied uniformly in this benchmark; thus coverage alters programme cost but not the disease-rate multiplier. Lowering the IRR multiplier corresponds to stronger protection. The 75%/125% ranges are declared perturbations, not verified confidence intervals. The current delayed-dose calendar application has a different target history, mechanism, price-year convention and normalization; these outcomes must not be substituted for its 000, 011 or 111 results.

<div id="tab:legacy_cohort_dsa">

| Parameter / level | Cases averted | Deaths averted | QALY gain | Cost (million) |
|:---|---:|---:|---:|---:|
| base_case / base | 254,454.2 | 1,069.3 | 30,777.1 | -79.574 |
| vaccine_price / low_75pct | 254,454.2 | 1,069.3 | 30,777.1 | -83.685 |
| mature_population_effect / low_75pct | 272,332.8 | 1,144.4 | 32,939.3 | -86.948 |
| CFR / low_75pct | 254,567.4 | 802.1 | 23,111.3 | -68.718 |
| sequela_cost / low_75pct | 254,454.2 | 1,069.3 | 30,777.1 | -79.540 |
| utility_loss / low_75pct | 254,454.2 | 1,069.3 | 30,752.6 | -79.574 |
| disease_burden / low_75pct | 190,925.5 | 802.1 | 23,086.2 | -53.357 |
| vaccine_price / high_125pct | 254,454.2 | 1,069.3 | 30,777.1 | -75.464 |
| mature_population_effect / high_125pct | 236,576.8 | 994.2 | 28,615.1 | -72.201 |
| CFR / high_125pct | 254,340.9 | 1,336.4 | 38,440.7 | -90.427 |
| sequela_cost / high_125pct | 254,454.2 | 1,069.3 | 30,777.1 | -79.608 |
| utility_loss / high_125pct | 254,454.2 | 1,069.3 | 30,801.7 | -79.574 |
| disease_burden / high_125pct | 317,926.2 | 1,336.4 | 38,465.8 | -105.778 |
| coverage / 80pct | 254,454.2 | 1,069.3 | 30,777.1 | -82.394 |
| coverage / 95pct | 254,454.2 | 1,069.3 | 30,777.1 | -78.165 |
| discount_rate / 0pct | 254,454.2 | 1,069.3 | 77,747.3 | -84.615 |
| discount_rate / 5pct | 254,454.2 | 1,069.3 | 20,010.7 | -76.511 |

Historical mature-effect PCV13 closed-cohort one-way sensitivity reconstructed from frozen code. Full-cohort totals use 1,384,610 births; cost is incremental societal cost in millions of 2024 USD. These are separate model-development evidence.

</div>

#### Exploratory probabilistic screening and evidence limits.

The development presentation reports an exploratory 200-draw PCV13 dynamic-cohort probabilistic analysis, with all sampled points in the dominant region and a cost-effectiveness acceptability curve (CEAC) of one on the tested threshold grid. The original draw-level output and run seed could not be recovered, and this exploratory analysis was not rerun or independently audited here. We preserve the record of that work without treating it as probabilistic evidence for the present age–calendar application. The current application remains a deterministic analysis; its terminal-occupancy samples and count-microsimulation replicates quantify simulation variation rather than joint parameter uncertainty.

#### SII bridge calculation and source separation.

The retained relative multiplier $`m_{\rm SII}=0.8203164086`$ is applied as $`\mathrm{VE}_{\rm SII,e}=m_{\rm SII}\mathrm{VE}_{\rm PCV13,e}`$, producing 0.49383048 for IPD, 0.15586012 for pneumonia and 0.37734555 for AOM. These three products are computationally reproducible from the frozen ledger. The associated workbook retains summary weights of 109% across visible categories, with 76 percentage points in shared vaccine serotypes and 33 in other categories; consequently $`76/109=0.6972477064`$. It also retains a weighted mean ratio of 0.8414473684. The original serotype-by-serotype geometric mean concentration (GMC) extraction needed to independently reconstruct the 0.8203164086 weighted-GMC summary is unavailable. Therefore the multiplication is verified, whereas its empirical weighting remains an unverified retained immunogenicity proxy. It is not a newly estimated clinical endpoint efficacy. Published immunogenicity studies motivate the bridge, but source identification alone does not verify the retained summary. The earlier GSK-based SII profile and the later PCV13-scaled profile are alternative historical assumptions, not interchangeable estimates.

<div id="tab:added_coverage">

| Research strand | Evidence status | Role in the present article |
|:---|:---|:---|
| Exact piecewise cohort propagation | Reconstructed from frozen source | Analytic/numerical verification, separate from current calendar application |
| Age-weight normalization and horizon-prefix invariance | Reconstructed from frozen source | Accounting and causal-horizon checks |
| Local cohort elasticities and structural screens | Reconstructed from frozen source | Early mechanism motivation; not calendar-PDE uncertainty estimates |
| Original PCV cohort characteristic validation | Reconstructed from frozen source | Same-input, same-estimand rewrite; no independent epidemiological validation |
| Independent simplified closed-cohort transport limit | Reconstructed and matched to saved historical file | Independent scalar versus age-bin bookkeeping; separate from original Week4 validation |
| Closed-cohort deterministic sensitivity | Reconstructed and matched to saved historical file | 2024-valued single-cohort mature-effect benchmark; not current rollout sensitivity |
| Exploratory cohort PSA | Historical report only; raw draws/seed not recovered; not rerun | Recorded historical work; excluded as uncertainty evidence for present calendar PDE |
| SII immunogenicity bridge | Computational multiplication verified; row-level original evidence unavailable | Retained proxy; full serotype/GMC extraction cannot be independently reconstructed |

Coverage and evidential roles of previously omitted model-development work. Frozen source files, reconstruction outputs and source hashes accompany the submission code.

</div>

## S14. Parameter provenance and transformations

The machine-readable parameter ledger separates three questions: whether a reference has the correct identity; whether its population, endpoint and observation period support a transfer; and whether the exact computational value can be reconstructed. A resolved DOI establishes the first question only. All retained input values, source roles and unresolved transfers are recorded in [`parameter_sources.csv`](../research/sources/parameter_sources.csv). The numerical comparisons and conversion identities in [`input_comparison.csv`](../research/sources/input_comparison.csv) are consistency checks, rather than independent clinical validation. Unavailable administrative costs and caregiver-survey aggregates are declared numerical scenario constants. The application consequently tests a modelling method with a transparent valuation ledger.

#### Disease-age allocation and exposure units.

For completed age month $`a\in\{0,\ldots,59\}`$, let $`w_e(a)`$ be a monthly allocation weight with $`\sum_{a=0}^{59}w_e(a)=1`$. The band masses and the uniform allocation within each band are given in Table <a href="#tab:age-allocation" data-reference-type="ref" data-reference="tab:age-allocation">39</a>. With a retained incidence input $`i_e`$ per 100,000 child-years, the computational monthly probability is
``` math
p_e^0(a)=1-\exp\{-5i_e w_e(a)/10^5\}.
\label{eq:age-incidence-allocation}
```
Thus $`5i_e/10^5`$ is the integrated five-year exposure for an ideal survivor under this allocation; $`i_e`$ itself is not a monthly probability. Probability modifiers are applied before conversion to the event hazard $`-\log(1-p_e)`$. The pneumonia and sepsis allocation is $`(0.50,0.25,0.15,0.10)`$ across the four bands; the meningitis allocation is $`(0.65,0.18,0.07,0.10)`$. These are frozen age-allocation scenarios. No independently verified national age-specific extract is available for these weights or the AOM allocation.

<div id="tab:age-allocation">

| Endpoint           | 0–11 | 12–23 | 24–35 | 36–59 |
|:-------------------|-----:|------:|------:|------:|
| Pneumonia / sepsis | 0.50 |  0.25 |  0.15 |  0.10 |
| Meningitis         | 0.65 |  0.18 |  0.07 |  0.10 |

Retained disease-age allocation. Each entry is the total mass in the specified age band; divide by its width in months for $`w_e(a)`$.

</div>

AOM uses 0.1869 at 0–11 months, 0.1889 at 12–23 months and 0.6242 at 24–59 months.

#### Acute duration and valuation.

For endpoint $`e`$, the undiscounted acute QALY loss is $`E_e d_e(1-u_e)/365`$, where $`E_e`$ is the event count, $`d_e`$ is duration in days and $`u_e`$ is the episode utility. The retained durations are 8.74 days for inpatient pneumonia, 7 for outpatient pneumonia, 18.41 for meningitis, 12.44 for sepsis and 7 for AOM. These durations were embedded in the economic configuration and are now exposed in the ledger; their exact primary estimates remain unverified. Utilities are transferred valuation scenarios contextualised by the original literature [Bennett2000,Oh1996,Oostenbrink2002,vanHoek2012](reference_key_map.csv). In Bennett and colleagues’ original elicitation, 0.9941 describes local infection, 0.9921 hospitalization, 0.9768 meningitis with recovery and 0.3903 severe brain damage [Bennett2000](reference_key_map.csv). Their present mappings to inpatient pneumonia, sepsis, meningitis and hemiplegia, respectively, are declared state-valuation transfers; the hemiplegia value specifically uses a severe-brain-damage proxy. The original chained-standard-gamble elicitation described whole disease outcomes with specified courses; its conversion to per-time health-state utility has not been independently established. Duration and transferred utility are consequently treated as separate scenario inputs. A duration and a utility are reported separately because the same average utility can imply a different episode loss when interpreted over a different observation window.

#### Cost-year transformation.

The inherited computational convention maps an amount labelled source year $`s>2023`$ to
``` math
c_{2023}^{\mathrm{proxy}}=c_s(1.03)^{-(s-2023)}.
\label{eq:proxy-cost-year}
```
Values with unspecified source year are retained unchanged. This 3% alignment factor is distinct from the 3% discount rate applied to future outcomes. It is a declared proxy and is not a reconstruction using a Vietnamese consumer-price index or an exchange-rate series. Medical and most sequela constants labelled 2024 are divided by 1.03; nonmedical, caregiver and special-education constants labelled 2025 are divided by $`1.03^2`$. The two PCV10 procurement scenarios use $`2/1.03=1.94174757`$ USD per dose. PCV13 uses the unchanged 4.22 USD constant whose source-year conversion is unresolved. Monetary results are consequently reported in scenario USD with the declared transformations, and a common real price year remains unresolved.

The delivery constant is $`2.10/1.03=2.03883495`$ USD per administered dose. Mai and colleagues examined tetanus-toxoid and tetanus–diphtheria delivery, providing context for programme costs [Mai2022](reference_key_map.csv); they do not independently identify this infant-PCV unit cost. The primary study reported fiscal estimates of 1.76, 3.86 and 1.49 USD per dose for facility, outreach and school delivery in its 2017 setting. Procurement cost is multiplied by $`(1+w)`$ with $`w=0.05`$, a retained mark-up convention. Delivery and adverse-event costs are applied to administered doses. Severe adverse-event cost 4,800.355 USD is the arithmetic midpoint of 4,717.59 and 4,883.12; the retained event rate is 4.01 per 100,000 doses. Original records and rate denominators require verification before policy use.

#### Lifetime weights and the calendar ledger.

Let $`d=0.03`$, assumed life expectancy $`L=74`$ years and assumed mean fatal-event age $`a_d=2.5`$ years. Write $`R=L-a_d`$, $`n=\lfloor R\rfloor`$ and $`\tau=R-n`$. The event-relative future QALY weight used by the code is
``` math
W_Q(a_d)=\sum_{k=0}^{n-1}(1+d)^{-(k+1/2)}
       +\tau(1+d)^{-(n+\tau/2)}=29.74234703.
\label{eq:remaining-life-weight}
```
Death losses occurring in calendar year $`y`$ receive the additional factor $`(1+d)^{-(y-1/2)}`$. This is a fixed life-expectancy approximation without a cause-specific survival table. Acute and sequela quantities are recorded within the ten-year programme window, whereas death losses receive this residual-life weight. The economic functional therefore has a mixed horizon; a fully lifetime sequela evaluation requires a post-window continuation model.

The 46 stored annual productivity contributions, for working ages $`k=15,\ldots,60`$, can be independently reconstructed as
``` math
q_k=\frac{12\,m\,(1-u_k)}{(1.03)^k},\qquad
 m=224.08583333333333,
\label{eq:productivity-weights}
```
where $`u_k=0.0792`$ at 15–24 years, $`0.0160`$ at 25–49 years and $`0.0141`$ at 50–60 years. Here $`m`$ is the unrounded monthly income inferred algebraically from the stored contributions. It should not be read as an independently verified national income observation. The raw sum is 42,392.92005858 USD; using the rounded displayed income 224.09 instead gives 42,393.70831531 USD. The code preserves the full stored array and uses $`Q_P=\sum_kq_k/1.03=41,158.17481416`$ scenario USD after the declared proxy adjustment. Because $`q_k`$ was discounted from birth, its fatal-event valuation is
``` math
W_P(a_d)=(1+d)^{a_d}\sum_{k>\lfloor a_d\rfloor}q_k/1.03
          =44,314.83832662.
\label{eq:event-productivity}
```
It then receives the same calendar-event discount as Eq. <a href="#eq:remaining-life-weight" data-reference-type="ref" data-reference="eq:remaining-life-weight">[eq:remaining-life-weight]</a>. The original statistical series, monetary conversion and employment-weight transfer have not been independently verified. exposes all contributions and the reconstruction error.

#### Pneumococcal AOM and source denominator.

The AOM endpoint is pneumococcal-attributable acute otitis media. Its baseline computational probability is $`p_{\mathrm{AOM}}^0(a)=\theta[1-\exp\{-5i_{\mathrm{AOM}}w_{\mathrm{AOM}}(a)/10^5\}]`$, with frozen proxy $`\theta=0.303`$. In the cited hospital study, 482 children aged at most five years had acute suppurative otitis media; pneumococcus represented 30.1% in the culture-related summary and 41.1% in the combined detection results [Do2025](reference_key_map.csv). The selected setting and laboratory denominators leave national causal attribution unidentified. Accordingly, 0.301 and 0.411 are operational sensitivity settings, not bounds of a population confidence interval. The exact 0.303 value was not independently recovered from the primary article. This distinction is retained once in the endpoint definition and parameter ledger.

#### Product bridging and evidence roles.

The Gambian trial supplies original evidence on serotype-specific post-booster immunogenicity [Adigweme2023](reference_key_map.csv). The retained SII multiplier 0.8203164085517026 is a summary proxy: the original row-level serotype weights, GMC extraction and aggregation rule needed to recover it are unavailable. Multiplying the PCV13 endpoint profile by that frozen constant is fully reproducible; identifying the constant from its claimed source is a separate unresolved step. The source data [`sii_bridge.csv`](../research/sources/sii_bridge.csv) record the reproducible endpoint multiplication. The visible weight totals, 76 and 109, give $`76/109=0.6972477064`$ but do not provide a validated compositional serotype distribution. Product-specific target histories, endpoint profiles and procurement proxies are therefore evaluated as three scenarios, without a head-to-head clinical ranking.

#### Input definitions and separate valuation states.

The computational ledger preserves endpoint-specific values at full precision. Outpatient-pneumonia nonmedical and caregiver constants are 18.98 and 124.49 USD before the disclosed conversion; the corresponding AOM constants are 21.21 and 139.12 USD. These distinct scenario inputs require original questionnaires and administrative records for empirical verification. Historical all-cause AOM and mature-effect output ledgers use different outcome definitions, denominators and accounting assumptions; their role is a separate benchmark.

<div id="tab:sii-bridge-reconstruction">

| Endpoint  | PCV13 profile | Summary multiplier |  SII profile |
|:----------|--------------:|-------------------:|-------------:|
| IPD       |         0.602 |       0.8203164086 | 0.4938304779 |
| Pneumonia |         0.190 |       0.8203164086 | 0.1558601176 |
| AOM       |         0.460 |       0.8203164086 | 0.3773455480 |

Reproducible multiplication of the retained SII summary multiplier. Original row-level identification of the multiplier remains unresolved.

</div>

## S15. Migration from the preceding manuscript

The preceding version is preserved unchanged. The revised calendar baseline reconciles mortality inputs and costs incident sequela flows. Numerical differences below result from those accounting changes; product inputs and D/W/G definitions are retained.

<div id="tab:ledger_map">

| Ledger | Protection convention | QALY gain |
|:---|:---|---:|
| Earlier age-profile on schedule | Archived earlier protection model | 1,811.1 (rounded) |
| Previous methods | DWG=000 | 1,707.471 |
| Previous methods | DWG=011 | 1,579.455 |
| Previous methods | DWG=111 | 1,517.943 |
| Revised methods | DWG=000 | 1,713.288 |
| Revised methods | DWG=011 | 1,584.403 |
| Revised methods | DWG=111 | 1,522.665 |

PCV13 QALY ledger map per 100,000 births. The earlier age-profile model uses a distinct protection implementation. Five-state 000, timing-reference 011 and full111 answer separate contrasts. The revised rows add the disclosed accounting corrections.

</div>

<div id="tab:migrationqalys_gained_per100k">

| Product / scenario         |  Previous |   Revised | Change | Change (%) |
|:---------------------------|----------:|----------:|-------:|-----------:|
| Product / scenario         |  Previous |   Revised | Change | Change (%) |
| No vaccination             |     0.000 |     0.000 |  0.000 |          – |
| PCV10-GSK / 000            | 1,583.209 | 1,588.599 |  5.390 |     0.3404 |
| PCV10-GSK / 001            | 1,491.024 | 1,495.785 |  4.761 |     0.3193 |
| PCV10-GSK / 010            | 1,570.857 | 1,576.151 |  5.294 |     0.3370 |
| PCV10-GSK / 011            | 1,478.594 | 1,483.257 |  4.664 |     0.3154 |
| PCV10-GSK / 100            | 1,548.516 | 1,553.801 |  5.285 |     0.3413 |
| PCV10-GSK / 101            | 1,433.330 | 1,437.867 |  4.537 |     0.3165 |
| PCV10-GSK / 110            | 1,537.938 | 1,543.141 |  5.203 |     0.3383 |
| PCV10-GSK / 111            | 1,422.563 | 1,427.015 |  4.452 |     0.3130 |
| PCV10-GSK / D0W0, herd off |   633.803 |   636.059 |  2.256 |     0.3560 |
| PCV10-GSK / D0W1, herd off |   604.096 |   606.120 |  2.024 |     0.3350 |
| PCV10-GSK / D1W0, herd off |   560.990 |   563.014 |  2.024 |     0.3608 |
| PCV10-GSK / D1W1, herd off |   535.560 |   537.386 |  1.826 |     0.3409 |
| PCV13 / 000                | 1,707.471 | 1,713.288 |  5.817 |     0.3407 |
| PCV13 / 001                | 1,592.861 | 1,597.914 |  5.053 |     0.3172 |
| PCV13 / 010                | 1,694.157 | 1,699.871 |  5.714 |     0.3373 |
| PCV13 / 011                | 1,579.455 | 1,584.403 |  4.948 |     0.3133 |
| PCV13 / 100                | 1,670.402 | 1,676.106 |  5.705 |     0.3415 |
| PCV13 / 101                | 1,529.560 | 1,534.373 |  4.813 |     0.3147 |
| PCV13 / 110                | 1,658.998 | 1,664.615 |  5.617 |     0.3385 |
| PCV13 / 111                | 1,517.943 | 1,522.665 |  4.722 |     0.3111 |
| PCV13 / D0W0, herd off     |   699.569 |   702.060 |  2.491 |     0.3561 |
| PCV13 / D0W1, herd off     |   666.811 |   669.045 |  2.235 |     0.3351 |
| PCV13 / D1W0, herd off     |   619.190 |   621.425 |  2.235 |     0.3609 |
| PCV13 / D1W1, herd off     |   591.149 |   593.165 |  2.016 |     0.3410 |
| PCV10-SII / 000            | 1,959.678 | 1,966.395 |  6.717 |     0.3428 |
| PCV10-SII / 001            | 1,669.818 | 1,674.940 |  5.122 |     0.3067 |
| PCV10-SII / 010            | 1,950.609 | 1,957.255 |  6.646 |     0.3407 |
| PCV10-SII / 011            | 1,660.639 | 1,665.688 |  5.048 |     0.3040 |
| PCV10-SII / 100            | 1,937.406 | 1,944.052 |  6.646 |     0.3430 |
| PCV10-SII / 101            | 1,608.147 | 1,613.026 |  4.879 |     0.3034 |
| PCV10-SII / 110            | 1,929.642 | 1,936.228 |  6.586 |     0.3413 |
| PCV10-SII / 111            | 1,600.137 | 1,604.952 |  4.815 |     0.3009 |
| PCV10-SII / D0W0, herd off |   573.850 |   575.894 |  2.043 |     0.3561 |
| PCV10-SII / D0W1, herd off |   546.980 |   548.813 |  1.833 |     0.3351 |
| PCV10-SII / D1W0, herd off |   507.918 |   509.751 |  1.833 |     0.3609 |
| PCV10-SII / D1W1, herd off |   484.916 |   486.570 |  1.654 |     0.3410 |

Previous and revised QALY gains per 100,000 births. Changes use the absolute previous denominator. A dash denotes an undefined ratio or percentage at a zero denominator.

</div>

<div id="tab:migrationincremental_costs_per100k">

| Product / scenario         | Previous | Revised | Change | Change (%) |
|:---------------------------|---------:|--------:|-------:|-----------:|
| Product / scenario         | Previous | Revised | Change | Change (%) |
| No vaccination             |    0.000 |   0.000 |  0.000 |          – |
| PCV10-GSK / 000            |   -4.377 |  -4.397 | -0.020 |    -0.4530 |
| PCV10-GSK / 001            |   -4.075 |  -4.093 | -0.018 |    -0.4337 |
| PCV10-GSK / 010            |   -4.303 |  -4.322 | -0.019 |    -0.4466 |
| PCV10-GSK / 011            |   -4.001 |  -4.018 | -0.017 |    -0.4263 |
| PCV10-GSK / 100            |   -4.396 |  -4.416 | -0.020 |    -0.4450 |
| PCV10-GSK / 101            |   -4.018 |  -4.035 | -0.017 |    -0.4225 |
| PCV10-GSK / 110            |   -4.333 |  -4.352 | -0.019 |    -0.4396 |
| PCV10-GSK / 111            |   -3.954 |  -3.970 | -0.016 |    -0.4159 |
| PCV10-GSK / D0W0, herd off |   -1.250 |  -1.258 | -0.008 |    -0.6538 |
| PCV10-GSK / D0W1, herd off |   -1.081 |  -1.088 | -0.007 |    -0.6266 |
| PCV10-GSK / D1W0, herd off |   -1.140 |  -1.147 | -0.007 |    -0.6481 |
| PCV10-GSK / D1W1, herd off |   -0.995 |  -1.002 | -0.006 |    -0.6226 |
| PCV13 / 000                |   -4.333 |  -4.354 | -0.021 |    -0.4907 |
| PCV13 / 001                |   -3.897 |  -3.915 | -0.018 |    -0.4660 |
| PCV13 / 010                |   -4.273 |  -4.294 | -0.021 |    -0.4862 |
| PCV13 / 011                |   -3.837 |  -3.855 | -0.018 |    -0.4604 |
| PCV13 / 100                |   -4.438 |  -4.459 | -0.021 |    -0.4766 |
| PCV13 / 101                |   -3.908 |  -3.925 | -0.018 |    -0.4489 |
| PCV13 / 110                |   -4.387 |  -4.408 | -0.021 |    -0.4727 |
| PCV13 / 111                |   -3.855 |  -3.873 | -0.017 |    -0.4438 |
| PCV13 / D0W0, herd off     |   -0.729 |  -0.736 | -0.007 |    -1.0053 |
| PCV13 / D0W1, herd off     |   -0.568 |  -0.574 | -0.006 |    -1.0566 |
| PCV13 / D1W0, herd off     |   -0.693 |  -0.699 | -0.007 |    -0.9604 |
| PCV13 / D1W1, herd off     |   -0.556 |  -0.561 | -0.006 |    -0.9948 |
| PCV10-SII / 000            |   -5.703 |  -5.729 | -0.026 |    -0.4478 |
| PCV10-SII / 001            |   -4.691 |  -4.710 | -0.019 |    -0.4120 |
| PCV10-SII / 010            |   -5.661 |  -5.686 | -0.025 |    -0.4450 |
| PCV10-SII / 011            |   -4.648 |  -4.667 | -0.019 |    -0.4082 |
| PCV10-SII / 100            |   -5.773 |  -5.798 | -0.025 |    -0.4411 |
| PCV10-SII / 101            |   -4.624 |  -4.643 | -0.019 |    -0.4018 |
| PCV10-SII / 110            |   -5.737 |  -5.762 | -0.025 |    -0.4387 |
| PCV10-SII / 111            |   -4.587 |  -4.605 | -0.018 |    -0.3983 |
| PCV10-SII / D0W0, herd off |   -0.859 |  -0.866 | -0.006 |    -0.7372 |
| PCV10-SII / D0W1, herd off |   -0.728 |  -0.733 | -0.005 |    -0.7216 |
| PCV10-SII / D1W0, herd off |   -0.792 |  -0.798 | -0.006 |    -0.7242 |
| PCV10-SII / D1W1, herd off |   -0.680 |  -0.685 | -0.005 |    -0.7083 |

Previous and revised Incremental costs (million scenario USD) per 100,000 births. Changes use the absolute previous denominator. A dash denotes an undefined ratio or percentage at a zero denominator.

</div>

<div id="tab:migrationcases_averted_per100k">

| Product / scenario         |   Previous |    Revised | Change | Change (%) |
|:---------------------------|-----------:|-----------:|-------:|-----------:|
| Product / scenario         |   Previous |    Revised | Change | Change (%) |
| No vaccination             |      0.000 |      0.000 |  0.000 |          – |
| PCV10-GSK / 000            | 14,538.735 | 14,609.530 | 70.795 |     0.4869 |
| PCV10-GSK / 001            | 13,910.182 | 13,975.886 | 65.704 |     0.4723 |
| PCV10-GSK / 010            | 14,169.753 | 14,237.411 | 67.658 |     0.4775 |
| PCV10-GSK / 011            | 13,540.269 | 13,602.814 | 62.545 |     0.4619 |
| PCV10-GSK / 100            | 14,198.684 | 14,267.830 | 69.145 |     0.4870 |
| PCV10-GSK / 101            | 13,388.189 | 13,451.036 | 62.847 |     0.4694 |
| PCV10-GSK / 110            | 13,884.290 | 13,950.772 | 66.482 |     0.4788 |
| PCV10-GSK / 111            | 13,071.188 | 13,131.328 | 60.140 |     0.4601 |
| PCV10-GSK / D0W0, herd off |  6,681.498 |  6,716.991 | 35.493 |     0.5312 |
| PCV10-GSK / D0W1, herd off |  5,863.542 |  5,892.057 | 28.514 |     0.4863 |
| PCV10-GSK / D1W0, herd off |  5,959.514 |  5,991.414 | 31.899 |     0.5353 |
| PCV10-GSK / D1W1, herd off |  5,262.768 |  5,288.744 | 25.977 |     0.4936 |
| PCV13 / 000                | 16,005.104 | 16,083.812 | 78.708 |     0.4918 |
| PCV13 / 001                | 14,852.154 | 14,921.700 | 69.545 |     0.4683 |
| PCV13 / 010                | 15,751.934 | 15,828.495 | 76.560 |     0.4860 |
| PCV13 / 011                | 14,597.868 | 14,665.241 | 67.372 |     0.4615 |
| PCV13 / 100                | 15,766.753 | 15,844.319 | 77.567 |     0.4920 |
| PCV13 / 101                | 14,361.986 | 14,428.817 | 66.830 |     0.4653 |
| PCV13 / 110                | 15,551.003 | 15,626.746 | 75.743 |     0.4871 |
| PCV13 / 111                | 14,143.304 | 14,208.263 | 64.958 |     0.4593 |
| PCV13 / D0W0, herd off     |  5,871.283 |  5,902.343 | 31.060 |     0.5290 |
| PCV13 / D0W1, herd off     |  5,157.801 |  5,182.776 | 24.975 |     0.4842 |
| PCV13 / D1W0, herd off     |  5,236.341 |  5,264.256 | 27.915 |     0.5331 |
| PCV13 / D1W1, herd off     |  4,628.575 |  4,651.326 | 22.751 |     0.4915 |
| PCV10-SII / 000            | 18,100.680 | 18,190.672 | 89.993 |     0.4972 |
| PCV10-SII / 001            | 15,568.124 | 15,640.004 | 71.880 |     0.4617 |
| PCV10-SII / 010            | 17,918.566 | 18,007.004 | 88.438 |     0.4936 |
| PCV10-SII / 011            | 15,384.798 | 15,455.097 | 70.299 |     0.4569 |
| PCV10-SII / 100            | 17,942.088 | 18,031.285 | 89.197 |     0.4971 |
| PCV10-SII / 101            | 15,059.628 | 15,128.650 | 69.022 |     0.4583 |
| PCV10-SII / 110            | 17,786.968 | 17,874.846 | 87.878 |     0.4941 |
| PCV10-SII / 111            | 14,901.419 | 14,969.072 | 67.653 |     0.4540 |
| PCV10-SII / D0W0, herd off |  4,816.538 |  4,842.018 | 25.480 |     0.5290 |
| PCV10-SII / D0W1, herd off |  4,231.191 |  4,251.679 | 20.488 |     0.4842 |
| PCV10-SII / D1W0, herd off |  4,295.655 |  4,318.555 | 22.900 |     0.5331 |
| PCV10-SII / D1W1, herd off |  3,797.038 |  3,815.702 | 18.663 |     0.4915 |

Previous and revised Cases averted per 100,000 births. Changes use the absolute previous denominator. A dash denotes an undefined ratio or percentage at a zero denominator.

</div>

<div id="tab:migrationdeaths_averted_per100k">

| Product / scenario         | Previous | Revised | Change | Change (%) |
|:---------------------------|---------:|--------:|-------:|-----------:|
| Product / scenario         | Previous | Revised | Change | Change (%) |
| No vaccination             |    0.000 |   0.000 |  0.000 |          – |
| PCV10-GSK / 000            |   62.496 |  62.709 |  0.212 |     0.3395 |
| PCV10-GSK / 001            |   59.285 |  59.475 |  0.190 |     0.3206 |
| PCV10-GSK / 010            |   61.997 |  62.205 |  0.208 |     0.3360 |
| PCV10-GSK / 011            |   58.783 |  58.969 |  0.186 |     0.3166 |
| PCV10-GSK / 100            |   61.154 |  61.362 |  0.208 |     0.3403 |
| PCV10-GSK / 101            |   57.110 |  57.292 |  0.182 |     0.3180 |
| PCV10-GSK / 110            |   60.726 |  60.931 |  0.205 |     0.3372 |
| PCV10-GSK / 111            |   56.676 |  56.854 |  0.178 |     0.3144 |
| PCV10-GSK / D0W0, herd off |   25.024 |  25.114 |  0.089 |     0.3573 |
| PCV10-GSK / D0W1, herd off |   23.817 |  23.897 |  0.080 |     0.3357 |
| PCV10-GSK / D1W0, herd off |   22.157 |  22.237 |  0.080 |     0.3621 |
| PCV10-GSK / D1W1, herd off |   21.123 |  21.196 |  0.072 |     0.3416 |
| PCV13 / 000                |   67.226 |  67.454 |  0.228 |     0.3397 |
| PCV13 / 001                |   63.234 |  63.435 |  0.201 |     0.3187 |
| PCV13 / 010                |   66.688 |  66.912 |  0.224 |     0.3362 |
| PCV13 / 011                |   62.692 |  62.889 |  0.197 |     0.3146 |
| PCV13 / 100                |   65.792 |  66.016 |  0.224 |     0.3405 |
| PCV13 / 101                |   60.852 |  61.044 |  0.192 |     0.3163 |
| PCV13 / 110                |   65.331 |  65.551 |  0.220 |     0.3374 |
| PCV13 / 111                |   60.383 |  60.572 |  0.189 |     0.3126 |
| PCV13 / D0W0, herd off     |   27.615 |  27.714 |  0.099 |     0.3571 |
| PCV13 / D0W1, herd off     |   26.283 |  26.372 |  0.088 |     0.3355 |
| PCV13 / D1W0, herd off     |   24.450 |  24.539 |  0.088 |     0.3619 |
| PCV13 / D1W1, herd off     |   23.311 |  23.390 |  0.080 |     0.3415 |
| PCV10-SII / 000            |   75.959 |  76.218 |  0.259 |     0.3408 |
| PCV10-SII / 001            |   65.979 |  66.183 |  0.204 |     0.3088 |
| PCV10-SII / 010            |   75.591 |  75.847 |  0.256 |     0.3386 |
| PCV10-SII / 011            |   65.607 |  65.807 |  0.201 |     0.3060 |
| PCV10-SII / 100            |   75.083 |  75.339 |  0.256 |     0.3410 |
| PCV10-SII / 101            |   63.688 |  63.882 |  0.195 |     0.3057 |
| PCV10-SII / 110            |   74.767 |  75.021 |  0.254 |     0.3392 |
| PCV10-SII / 111            |   63.363 |  63.555 |  0.192 |     0.3031 |
| PCV10-SII / D0W0, herd off |   22.652 |  22.733 |  0.081 |     0.3571 |
| PCV10-SII / D0W1, herd off |   21.560 |  21.632 |  0.072 |     0.3355 |
| PCV10-SII / D1W0, herd off |   20.057 |  20.129 |  0.073 |     0.3619 |
| PCV10-SII / D1W1, herd off |   19.122 |  19.187 |  0.065 |     0.3415 |

Previous and revised Deaths averted per 100,000 births. Changes use the absolute previous denominator. A dash denotes an undefined ratio or percentage at a zero denominator.

</div>

<div id="tab:migrationicer">

| Product / scenario         |   Previous |    Revised | Change | Change (%) |
|:---------------------------|-----------:|-----------:|-------:|-----------:|
| Product / scenario         |   Previous |    Revised | Change | Change (%) |
| No vaccination             |          – |          – |      – |          – |
| PCV10-GSK / 000            | -2,764.714 | -2,767.817 | -3.102 |    -0.1122 |
| PCV10-GSK / 001            | -2,733.329 | -2,736.444 | -3.115 |    -0.1140 |
| PCV10-GSK / 010            | -2,739.310 | -2,742.303 | -2.993 |    -0.1093 |
| PCV10-GSK / 011            | -2,706.019 | -2,709.010 | -2.991 |    -0.1105 |
| PCV10-GSK / 100            | -2,838.824 | -2,841.759 | -2.936 |    -0.1034 |
| PCV10-GSK / 101            | -2,802.995 | -2,805.957 | -2.962 |    -0.1057 |
| PCV10-GSK / 110            | -2,817.245 | -2,820.088 | -2.843 |    -0.1009 |
| PCV10-GSK / 111            | -2,779.229 | -2,782.081 | -2.852 |    -0.1026 |
| PCV10-GSK / D0W0, herd off | -1,972.390 | -1,978.241 | -5.852 |    -0.2967 |
| PCV10-GSK / D0W1, herd off | -1,789.326 | -1,794.526 | -5.200 |    -0.2906 |
| PCV10-GSK / D1W0, herd off | -2,031.822 | -2,037.638 | -5.817 |    -0.2863 |
| PCV10-GSK / D1W1, herd off | -1,858.710 | -1,863.928 | -5.218 |    -0.2807 |
| PCV13 / 000                | -2,537.740 | -2,541.535 | -3.795 |    -0.1495 |
| PCV13 / 001                | -2,446.691 | -2,450.320 | -3.629 |    -0.1483 |
| PCV13 / 010                | -2,522.288 | -2,526.031 | -3.743 |    -0.1484 |
| PCV13 / 011                | -2,429.269 | -2,432.832 | -3.563 |    -0.1467 |
| PCV13 / 100                | -2,656.872 | -2,660.449 | -3.577 |    -0.1346 |
| PCV13 / 101                | -2,554.663 | -2,558.082 | -3.419 |    -0.1338 |
| PCV13 / 110                | -2,644.261 | -2,647.795 | -3.534 |    -0.1337 |
| PCV13 / 111                | -2,539.900 | -2,543.262 | -3.362 |    -0.1324 |
| PCV13 / D0W0, herd off     | -1,041.606 | -1,048.345 | -6.738 |    -0.6469 |
| PCV13 / D0W1, herd off     |   -852.243 |   -858.371 | -6.128 |    -0.7190 |
| PCV13 / D1W0, herd off     | -1,118.853 | -1,125.537 | -6.684 |    -0.5974 |
| PCV13 / D1W1, herd off     |   -940.323 |   -946.450 | -6.127 |    -0.6516 |
| PCV10-SII / 000            | -2,910.224 | -2,913.270 | -3.046 |    -0.1047 |
| PCV10-SII / 001            | -2,809.106 | -2,812.055 | -2.949 |    -0.1050 |
| PCV10-SII / 010            | -2,902.164 | -2,905.180 | -3.016 |    -0.1039 |
| PCV10-SII / 011            | -2,799.031 | -2,801.938 | -2.907 |    -0.1038 |
| PCV10-SII / 100            | -2,979.534 | -2,982.445 | -2.911 |    -0.0977 |
| PCV10-SII / 101            | -2,875.312 | -2,878.132 | -2.821 |    -0.0981 |
| PCV10-SII / 110            | -2,972.890 | -2,975.776 | -2.886 |    -0.0971 |
| PCV10-SII / 111            | -2,866.643 | -2,869.426 | -2.783 |    -0.0971 |
| PCV10-SII / D0W0, herd off | -1,497.486 | -1,503.173 | -5.688 |    -0.3798 |
| PCV10-SII / D0W1, herd off | -1,330.499 | -1,335.624 | -5.125 |    -0.3852 |
| PCV10-SII / D1W0, herd off | -1,559.682 | -1,565.328 | -5.646 |    -0.3620 |
| PCV10-SII / D1W1, herd off | -1,402.044 | -1,407.175 | -5.132 |    -0.3660 |

Previous and revised ICER (scenario USD/QALY) in scenario USD per QALY. Changes use the absolute previous denominator. A dash denotes an undefined ratio or percentage at a zero denominator.

</div>

## S16. Probability-to-hazard conversion, nominal CFR and sequela-source semantics

This added audit describes the frozen calculation precisely without changing its coefficients, trajectories or economic results. For a state-specific monthly disease probability $`p`$ and nominal case-fatality input $`f`$, the implementation separately constructs the episode intensity $`\lambda=-\log(1-p)`$ and disease-death intensity $`\mu=-\log(1-fp)`$. Its episode and death rewards are integrated against the same alive-state exposure. Consequently their ratio is $`f_{\mathrm{eff}}=\mu/\lambda`$, which is generally different from $`f`$. For $`0<p<1`$ and $`0<f<1`$, convexity of $`-\log(1-p)`$ gives $`\mu<f\lambda`$. The frozen CFR is therefore a nominal parameter in a probability conversion; it is not exactly the conditional episode-to-death proportion of a joint disease-course model. In the small-probability limit, $`f_{\mathrm{eff}}/f=1-(1-f)p/2+O(p^2)`$.

Pneumonia episodes are split into severe and non-severe rewards using the fixed severe fraction $`s`$. The death-probability conversion uses $`p\,s\,f_{\mathrm{sev}}`$, whereas the severe-episode intensity is $`s\lambda`$. Its effective severe CFR is therefore $`\mu/(s\lambda)`$. Meningitis and sepsis use the unsplit expression above; AOM has no disease-death reward in this audit. No new clinical CFR evidence is inferred from these mathematical relationships.

#### Frozen-exposure diagnostic and reconstruction checks.

The executable [CFR translation audit](../checks/revision_20261103/cfr_translation_audit.py) reads all 144 saved configurations: 120 monthly states per configuration, each with five alive dose states and 60 attained-age bins. It freezes the stored endpoint-specific $`H`$ and target $`H^*`$ values, including the fixed-target direct-VE perturbations, and rebuilds the original monthly coefficients in memory. For saved post-birth, post-dose state $`A`$, the one-month constant-coefficient exposure is

``` math
E=A\frac{1-e^{-\mu_{\rm all}}}{\mu_{\rm all}},
\qquad E=A\quad\text{when }\mu_{\rm all}=0.
```

The audit uses `expm1` for this expression and to recover $`p=-\operatorname{expm1}(-\lambda)`$, then `log1p` to reconstruct the component death intensities. These stable functions are confined to the diagnostic; they do not replace the frozen solver's arithmetic. The original aggregate death coefficient is used for the direct saved-flow reconstruction. Across 207,360 monthly comparisons for six disease/death rewards and six incident-sequela rewards, the largest absolute discrepancy was $`1.054\times10^{-10}`$ events and the largest relative discrepancy at a nonzero reference was $`5.4985\times10^{-13}`$. All 1,728 configuration-by-metric checks met the unchanged $`\mathrm{rtol}=10^{-12}`$, $`\mathrm{atol}=10^{-7}`$ criteria; none failed. The separate stable reconstruction of component death hazards differs from the original summed death coefficient by at most $`1.627\times10^{-16}`$ per state-month; cumulative reconstructed component deaths differ from the saved total by at most $`3.810\times10^{-8}`$ events. Input and model-source SHA-256 values are verified unchanged before and after the audit. The audit's 2,595 checks are identities and reconstruction checks, not independent observations or a new full trajectory rerun.

#### Magnitude under the saved trajectories.

For each state and month with positive exposure and episode probability, the relative difference $`f_{\mathrm{eff}}/f-1`$ has the following extrema over all 144 configurations. For severe pneumonia, $`f_{\mathrm{eff}}`$ is defined relative to severe episodes and $`f`$ is the nominal severe CFR. These are dimensionless ratios, not percentage points.

| Endpoint | Minimum $`f_{\mathrm{eff}}/f-1`$ | Maximum $`f_{\mathrm{eff}}/f-1`$ |
|:---|---:|---:|
| Meningitis | $`-1.9657010\times10^{-5}`$ | $`-1.7127199\times10^{-7}`$ |
| Severe pneumonia | $`-2.2897325\times10^{-3}`$ | $`-3.0116278\times10^{-5}`$ |
| Sepsis | $`-1.4852824\times10^{-4}`$ | $`-1.6824814\times10^{-6}`$ |

Replacing only the death reward by $`f\lambda E`$ on the same frozen exposures, while retaining the original trajectories, gives the diagnostic totals below over the 120-month recorded window. The denominator is 13,846,100 gross recorded-window births. The frozen death column uses the component-hazard reconstruction; it agrees with the saved totals to numerical precision. The percentage column is the hypothetical reward increase relative to the reconstructed frozen death reward.

| Configuration | Frozen disease-death reward | CFR-times-episode reward on the same exposure | Difference | Relative increase (%) | Difference per 100,000 recorded-window births |
|:---|---:|---:|---:|---:|---:|
| No vaccination | 13,737.60998 | 13,746.52425 | 8.91427 | 0.064890 | 0.064381 |
| PCV10-GSK, D/W/G=111 | 5,865.58364 | 5,867.65475 | 2.07111 | 0.035310 | 0.014958 |
| PCV13, D/W/G=111 | 5,350.78557 | 5,352.61724 | 1.83167 | 0.034232 | 0.013229 |
| PCV10-SII, D/W/G=111 | 4,937.68391 | 4,939.29781 | 1.61390 | 0.032685 | 0.011656 |

Across all 144 configurations the hypothetical increase ranges from 0.014384% to 0.081108%, or 0.003353 to 0.100562 deaths per 100,000 recorded-window births. The configurations include repeated trajectories with valuation-only perturbations; these ranges do not describe independent empirical samples. Replacing a death hazard in a complete model would also change survival, future exposure and indirect-protection feedback. These frozen-exposure differences therefore do not estimate the effect of a fully revised model, its QALYs or its costs.

#### Nominal-CFR sequela sources.

The incident meningitis sequela-source proxy uses $`(1-f)\lambda E`$ before multiplying by the outcome-specific sequela probabilities. It is not algebraically identical to $`(\lambda-\mu)E`$. The latter minus the former is $`E(f\lambda-\mu)`$. Summed over the recorded window, this source difference is 0.002979 individuals for no vaccination, 0.000712 for PCV10-GSK full111, 0.000638 for PCV13 full111 and 0.000523 for PCV10-SII full111; over all 144 configurations it ranges from 0.000202 to 0.004654 individuals. These are differences between two source conventions, before outcome-specific sequela probabilities and downstream survival adjustments. They are not observed survivor counts or a recomputation of sequela stocks, QALYs or costs. The source and death-reward conventions remain disclosed approximations rather than a jointly calibrated clinical progression model.

Full precision, extreme-cell locations, endpoint totals and strict reconstruction residuals are available in [the audit JSON](../checks/revision_20261103/cfr_translation_audit.json), [all configurations](../checks/revision_20261103/cfr_translation_runs.csv), [the four core configurations](../checks/revision_20261103/cfr_translation_core.csv), [endpoint totals](../checks/revision_20261103/cfr_translation_endpoints.csv), [extreme cells](../checks/revision_20261103/cfr_translation_extrema.csv) and [reconstruction checks](../checks/revision_20261103/cfr_translation_reconstruction.csv). Reproduce the audit from the repository root with `python checks/revision_20261103/cfr_translation_audit.py`.
