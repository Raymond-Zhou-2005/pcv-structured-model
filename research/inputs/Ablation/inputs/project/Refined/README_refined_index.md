# Refined PCV Model Index

Updated after folder organization on 2026-06-23.

## Naming Convention

- Week 1-4 files are preserved as historical development files.
- Week 5 and later files use shorter feature-based names while keeping the `Week N_...` prefix.
- Generated outputs are consolidated under `Refined/Outputs` instead of being stored as separate development subfolders.
- Product naming uses `PCV10-SII` in current model/reporting assets.

## Main Refined Files

- `Week 1.pdf`
- `Week 1_assumption_tiers.py`
- `Week 1_exact_sensitivity.py`
- `Week 2.pdf`
- `Week 2_markov_refined.py`
- `Week 2_ode_refined.py`
- `Week 2_parameter_table.xlsx`
- `Week 2_qaly_refined.py`
- `Week 3.pdf`
- `Week 3_age_time_original.py`
- `Week 3_age_time_refined.py`
- `Week 4.pdf`
- `Week 4_calendar_original.py`
- `Week 4_pde_validation.py`
- `Week 5_calendar_pde_original.py`
- `Week 5_calendar_pde_parameter_table.xlsx`
- `Week 5_calendar_pde_refined.py`
- `Week 5_calendar_pde_rollout_to_stable.pdf`
- `Week 6_calendar_pde_program_realism.py`
- `Week 6_calendar_pde_program_realism_update.pdf`
- `Week 6_calendar_pde_program_realism_update_English.pptx`
- `Week 6_calendar_pde_program_realism_update_original.pptx`
- `Week 7_mechanism_model_refinement_qaly_impact.pdf`
- `Week 7_mechanism_qaly_impact_brief_original.pptx`
- `Week 7_mechanism_qaly_impact_brief_revised.pptx`
- `Week 7_qaly_refinement_impact.py`
- `Week 7_qaly_refinement_methods_and_results.md`

## Consolidated Outputs

- `Week 6_calendar_pde_program_realism`: 27 files
- `Week 7_qaly_refinement`: 23 files

Detailed output inventory: `Outputs/refined_outputs_index.csv`.

## Refined Reference Files

- `Artical/refined_article_index.csv` indexes the newly added mechanism and Vietnam-specific PCV references used to update the refined manuscript references.
- The three Vietnam-specific PDF sources are high relevance to the refined model because they inform carriage calibration, schedule/dose-state assumptions, and endpoint-specific direct VE refinement.
- The US Pharmacist and StatPearls/PMID 29939571 files are background references for vaccine products and mechanism, not Vietnam-specific parameter sources.

## Notes

- `Week 6_calendar_pde_program_realism.py` is the current most advanced epidemiological model.
- `Week 7_qaly_refinement_impact.py` builds QALY and cost-effectiveness comparisons from original QALY outputs and consolidated calendar PDE outputs.
- `Paper_Release` remains the publication-oriented bundle and has been patched to read the organized output paths.
