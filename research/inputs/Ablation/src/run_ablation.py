"""Matched 2^3 programme-realism experiment with explicit accounting corrections.

No original model file is modified. G=0 preserves annual IRR targets and
removes only their additional carriage ramp; G=1 is original Week 6.
Frozen Refined inputs are unchanged. The adapters below cause-delete baseline
disease mortality from a constant postneonatal all-cause rate and charge one-off
sequela costs to incident flows. These are disclosed accounting conventions;
the background profile is not an empirically identified cause-deleted life table.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "inputs" / "project"
OUTPUT = ROOT / "outputs"
PRODUCTS = {"PCV10-GSK": ("pcv10_gsk", "pcv10"),
            "PCV13": ("pcv13", "pcv13"), "PCV10-SII": ("pcv10sii", "pcv10sii")}
BITS = list(itertools.product((0, 1), repeat=3))

def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

econ = load_module(FROZEN / "Refined" / "Week 7_qaly_refinement_impact.py", "release_economics")
w6 = econ.delayed_dose_source
w5 = econ.on_schedule_source
BASE_P = copy.deepcopy(w6.params_vietnam_micro)
ECON_P = econ.build_calendar_params()

# Calibrate once at the frozen base incidence/CFR values, then hold this
# background fixed when those disease inputs are perturbed.
_INHERITED_COEFFICIENTS = w6.build_week6_month_coefficients
_Q_NEONATAL = 1.0 - w6.base.neonatal_survival_multiplier(BASE_P)
_Q_U5 = float(BASE_P["u5_mortality_rate"]) / 1000.0
_POSTNEONATAL_HAZARD = -math.log((1.0 - _Q_U5) / (1.0 - _Q_NEONATAL)) / 60.0
_NV_DISEASE_HAZARD = _INHERITED_COEFFICIENTS(
    BASE_P, 0, np.ones((w6.N_GROUPS, 60)), False, None)["mu_dis"][0].copy()
_BACKGROUND_HAZARD = _POSTNEONATAL_HAZARD - _NV_DISEASE_HAZARD
if np.any(_BACKGROUND_HAZARD < 0):
    raise ValueError("Cause-deleted background would be negative; revise the allocation convention.")
BASE_P["mortality_accounting"] = "cause_deleted"
MORTALITY_ACCOUNTING = {
    "convention": "constant_postneonatal_all_cause_minus_frozen_NV_disease_hazard",
    "neonatal_probability": _Q_NEONATAL, "under_five_probability": _Q_U5,
    "postneonatal_monthly_all_cause_hazard": _POSTNEONATAL_HAZARD,
    "baseline_NV_disease_hazard_by_month_of_age": _NV_DISEASE_HAZARD.tolist(),
    "cause_deleted_background_hazard_by_month_of_age": _BACKGROUND_HAZARD.tolist(),
    "background_fixed_during_disease_perturbations": True,
    "status": "accounting convention, not an identified empirical life table",
}

def revised_coefficients(p, calendar_month, A_start, intervention=False, vaccine_type=None):
    cc = _INHERITED_COEFFICIENTS(p, calendar_month, A_start, intervention, vaccine_type)
    mode = p.get("mortality_accounting", "cause_deleted")
    if mode == "cause_deleted":
        cc["mu_all"] = cc["mu_dis"] + _BACKGROUND_HAZARD[None, :]
    elif mode == "historical_neonatal_conditioned":
        cc["mu_all"] = cc["mu_dis"] + _POSTNEONATAL_HAZARD
    elif mode == "zero_background":
        cc["mu_all"] = cc["mu_dis"].copy()
    elif mode != "historical_inherited":
        raise ValueError(f"Unknown mortality accounting mode: {mode}")
    return cc

w6.build_week6_month_coefficients = revised_coefficients

def incident_sequela_flows(p, flows):
    """Incident source counts; distinct from surviving stock increments."""
    m = (1.0 - float(p["cfr_meningitis"])) * flows["New_Meningitis"]
    hearing = float(p["prob_hear"]) * m
    return {
        "New_Seq_Cognitive": float(p["prob_cog"]) * m,
        "New_Seq_Hearing": hearing,
        "New_Seq_Epilepsy": float(p["prob_epi"]) * m,
        "New_Seq_Hemiplegia": float(p["prob_hemi"]) * m,
        "New_Seq_Hearing_Implant": float(p.get("prob_cochlear_implant", 0.0)) * hearing,
        "New_Seq_Hearing_SpEd": float(p.get("prob_special_education", 0.0)) * hearing,
    }

_INHERITED_ECONOMICS = econ.compute_calendar_economics
def corrected_calendar_economics(model_label, scenario, vaccine_label, vaccine_key,
                                 epi, vaccine_costs, params):
    details = _INHERITED_ECONOMICS(model_label, scenario, vaccine_label, vaccine_key,
                                  epi, vaccine_costs, params)
    mapping = {"hearing": ("New_Seq_Hearing", "hearing_case"),
               "epilepsy": ("New_Seq_Epilepsy", "epilepsy_case"),
               "hemiplegia": ("New_Seq_Hemiplegia", "hemiplegia_case"),
               "implant": ("New_Seq_Hearing_Implant", "cochlear_implant_case")}
    for ix, row in epi.reset_index(drop=True).iterrows():
        correction = 0.0
        for label, (column, cost) in mapping.items():
            if column not in row:
                raise ValueError("Revised economic ledger requires actual incident sequela flows")
            old_column = f"new_seq_{label}_proxy"
            old = float(details.loc[ix, old_column])
            incident = float(row[column])
            correction += (incident - old) * float(params["costs_sequelae"][cost])
            details.loc[ix, f"legacy_net_increment_seq_{label}"] = old
            details.loc[ix, f"new_seq_{label}_incident"] = incident
        event_df = econ.discount_factor(int(row["Year"]) - 0.5, float(params["discount_rate"]))
        details.loc[ix, "raw_sequela_cost"] += correction
        details.loc[ix, "raw_total_cost"] += correction
        details.loc[ix, "pv_sequela_cost"] += correction * event_df
        details.loc[ix, "pv_total_cost"] += correction * event_df
    details = details.drop(columns=[f"new_seq_{x}_proxy" for x in mapping])
    details["sequela_oneoff_counting"] = "actual_incident_source_flow"
    return details

econ.compute_calendar_economics = corrected_calendar_economics

def scenario_id(vaccine, d, w, g):
    return f"{PRODUCTS[vaccine][0]}_D{d}W{w}G{g}"

def parameters(d=0, w=0):
    p = copy.deepcopy(BASE_P)
    p["record_months"] = 120
    if not d:
        p["routine_dose1_transition_probs"] = {2: 1.0}
        p["routine_dose2_transition_probs"] = {4: 1.0}
        p["routine_booster_transition_probs"] = {12: 1.0}
    if not w:
        p["primary_waning_by_age"] = {}
        p["booster_waning_by_age"] = {}
    return p

def yearly_unrounded(monthly):
    """Same flow sum and year-end stock convention as Refined, no round(0)."""
    stocks = [c for c in w6.OUTPUT_COLUMNS if c.startswith("Stock_")]
    flows = [c for c in monthly if c.startswith("New_")]
    annual = monthly.groupby("Year")[flows].sum()
    annual[stocks] = monthly.groupby("Year")[stocks].last()
    annual["Year"] = annual.index.astype(int)
    extra = [c for c in flows if c not in w6.OUTPUT_COLUMNS]
    return annual.reset_index(drop=True)[w6.OUTPUT_COLUMNS + extra]

def simulate(vaccine, d, w, g, name, *, herd_off=False, rtol=1e-8, atol=1e-7):
    p = parameters(d, w)
    shares = w6.base.normalize_birth_shares(p["birth_month_shares"])
    A, S, warm = w6.build_week6_initial_population(p, shares, name, False)
    original_herd = w6.herd_multiplier_for_endpoint
    def selected_herd(pp, vv, endpoint, year, state):
        actual, target, relative = original_herd(pp, vv, endpoint, year, state)
        return (1.0 if herd_off else actual if g else target), target, relative
    w6.herd_multiplier_for_endpoint = selected_herd
    history, dose_history, herd_history, states, balance = [], [], [], [], []
    cumulative = {k: 0.0 for k in w6.FLOW_NAMES}
    cumulative_all_deaths = 0.0
    try:
        for month in range(p["record_months"]):
            initial_alive = float(A.sum())
            gross_births = p["birth_cohort"] * shares[month % 12]
            live_births, neonatal_deaths = w6.add_record_births_week6(A, gross_births, p, vaccine is not None, vaccine)
            before_events = float(A.sum())
            counts = econ.apply_delayed_vaccination_events_counted(A, p, month, vaccine is not None, vaccine)
            after_events = float(A.sum())
            diag = w6._dose_diagnostics(A, p)
            states.append(A.copy())
            herd_row = {"Month": month, "Year": month // 12 + 1,
                        "carriage_index": diag["Carriage_Index"],
                        "relative_carriage_index": diag["Relative_Carriage_Index"],
                        "applied_build_fraction": 0.0 if herd_off or vaccine is None else diag["Relative_Carriage_Index"] if g else 1.0}
            for endpoint in ("ipd", "pneumonia", "aom"):
                actual, target, _ = selected_herd(p, vaccine, endpoint, month // 12, A)
                herd_row[f"H_{endpoint}"] = actual
                herd_row[f"H_target_{endpoint}"] = target
            herd_history.append(herd_row)
            A_after, S, flows, all_deaths = w6.solve_one_month_week6(
                A, S, month, p, vaccine is not None, vaccine, rtol=rtol, atol=atol)
            cumulative_all_deaths += neonatal_deaths + all_deaths
            for key in cumulative:
                cumulative[key] += flows[key]
            row = {"Month": month, "Year": month // 12 + 1,
                   **flows, **incident_sequela_flows(p, flows), "Stock_Alive": float(A_after.sum()),
                   "Stock_Dead": cumulative_all_deaths,
                   **{f"Stock_{k}": v for k, v in cumulative.items()}, **S,
                   "Monthly_Total_Cases": sum(flows[k] for k in w6.FLOW_NAMES if k != "New_Deaths_Pneumo"),
                   **diag, "Gross_Births": gross_births}
            history.append(row)
            dose_history.append({"Month": month, "Year": month // 12 + 1, **counts})
            A, aged_out = w6.age_shift_week6(A_after)
            balance.append({"Month": month,
                            "population_balance_error": float(A.sum()) - (initial_alive + live_births - all_deaths - aged_out),
                            "vaccination_balance_error": after_events - before_events,
                            "minimum_state": float(A.min())})
    finally:
        w6.herd_multiplier_for_endpoint = original_herd
    return {"monthly": pd.DataFrame(history), "doses": pd.DataFrame(dose_history),
            "herd": pd.DataFrame(herd_history), "balance": pd.DataFrame(balance),
            "states": np.asarray(states), "parameters": p, "warmup": warm}

def economics(result, vaccine, name, legacy_rounding=False):
    yearly = yearly_unrounded(result["monthly"])
    if legacy_rounding:
        yearly = yearly.round(0)
    counts = {int(year): grp[["dose1", "dose2", "booster"]].sum().to_dict()
              for year, grp in result["doses"].groupby("Year")}
    key = PRODUCTS[vaccine][1] if vaccine else None
    label = econ.MODEL_CALENDAR_DELAYED
    costs = pd.DataFrame(econ.cost_rows_from_dose_counts(label, name, vaccine or "Baseline", key, counts, ECON_P))
    details = econ.compute_calendar_economics(label, name, vaccine or "Baseline", key, yearly, costs, ECON_P)
    return yearly, costs, details

def summarize(result, details):
    m = result["monthly"]
    return {"cases_total": float(m["Monthly_Total_Cases"].sum()),
            "deaths_total": float(m["New_Deaths_Pneumo"].sum()),
            "qaly_loss": float(details["pv_total_qaly_loss"].sum()),
            "costs_total": float(details["pv_total_cost"].sum()),
            "vaccine_costs": float(details["pv_vaccine_program_cost"].sum()),
            "disease_costs": float(details["pv_disease_cost"].sum()),
            "sequela_costs": float(details["pv_sequela_cost"].sum()),
            "productivity_costs": float(details["pv_death_prod_cost"].sum()),
            "doses_total": float(details["administered_doses"].sum()),
            "births_denominator": float(m["Gross_Births"].sum())}

def save_result(name, result):
    dest = OUTPUT / "runs" / name
    dest.mkdir(parents=True, exist_ok=True)
    for key in ("monthly", "doses", "herd", "balance"):
        result[key].to_csv(dest / f"{key}.csv", index=False)
    np.savez_compressed(dest / "states.npz", alive_after_births_and_doses=result["states"],
                        age_month=np.arange(60), group_names=np.asarray(w6.GROUP_NAMES))
    (dest / "parameters.json").write_text(json.dumps(result["parameters"], indent=2), encoding="utf-8")
    (dest / "warmup.json").write_text(json.dumps(result["warmup"], indent=2), encoding="utf-8")
    (dest / "mortality_accounting.json").write_text(json.dumps(MORTALITY_ACCOUNTING, indent=2), encoding="utf-8")

def add_contrasts(frame):
    frame = frame.copy()
    nv = frame.loc[frame["scenario"] == "no_vaccination"].iloc[0]
    for metric in ("cases_total", "deaths_total", "qaly_loss", "costs_total", "vaccine_costs", "disease_costs", "sequela_costs", "productivity_costs", "doses_total"):
        frame[f"{metric}_per100k"] = frame[metric] / frame["births_denominator"] * 100000.0
    frame["cases_averted"] = nv.cases_total - frame.cases_total
    frame["deaths_averted"] = nv.deaths_total - frame.deaths_total
    frame["qalys_gained"] = nv.qaly_loss - frame.qaly_loss
    frame["incremental_costs"] = frame.costs_total - nv.costs_total
    for metric in ("cases_averted", "deaths_averted", "qalys_gained", "incremental_costs"):
        frame[f"{metric}_per100k"] = frame[metric] / frame["births_denominator"] * 100000.0
    frame["icer"] = np.where(abs(frame.qalys_gained) > 1e-10, frame.incremental_costs / frame.qalys_gained, np.nan)
    frame["economic_status"] = np.select(
        [(frame.qalys_gained > 0) & (frame.incremental_costs < 0),
         (frame.qalys_gained < 0) & (frame.incremental_costs > 0)],
        ["dominant_vs_no_vaccination", "dominated_vs_no_vaccination"], default="tradeoff_or_reference")
    metrics = ["cases_total_per100k", "deaths_total_per100k", "qaly_loss_per100k", "costs_total_per100k", "vaccine_costs_per100k",
               "cases_averted_per100k", "deaths_averted_per100k", "qalys_gained_per100k", "incremental_costs_per100k", "icer"]
    for vaccine in PRODUCTS:
        b = frame.loc[frame.scenario == scenario_id(vaccine, 0, 0, 0)].iloc[0]
        mask = frame.vaccine == vaccine
        for metric in metrics:
            frame.loc[mask, f"{metric}_delta_vs_on_schedule"] = frame.loc[mask, metric] - b[metric]
            frame.loc[mask, f"{metric}_pct_vs_on_schedule"] = 100 * (frame.loc[mask, metric] - b[metric]) / abs(b[metric]) if abs(b[metric]) > 1e-12 else np.nan
    return frame

def decompose(frame):
    interactions, shapley, indirect = [], [], []
    metrics = ["cases_total_per100k", "deaths_total_per100k", "qaly_loss_per100k", "costs_total_per100k",
               "cases_averted_per100k", "deaths_averted_per100k", "qalys_gained_per100k", "incremental_costs_per100k", "icer"]
    for vaccine in PRODUCTS:
        rows = {(d, w, g): frame.loc[frame.scenario == scenario_id(vaccine, d, w, g)].iloc[0] for d, w, g in BITS}
        base = rows[(0, 0, 0)]
        for metric in metrics:
            betas = {}
            for subset in BITS[1:]:
                val = 0.0
                for lower in BITS:
                    if all(l <= s for l, s in zip(lower, subset)):
                        val += (-1) ** (sum(subset) - sum(lower)) * rows[lower][metric]
                betas[subset] = val
                interactions.append({"vaccine": vaccine, "metric": metric,
                    "term": "".join(c for c, active in zip("DWG", subset) if active),
                    "order": sum(subset), "interaction_delta": val,
                    "pct_of_abs_on_schedule": 100 * val / abs(base[metric])})
            for i, factor in enumerate("DWG"):
                value = sum(v / sum(s) for s, v in betas.items() if s[i])
                joint = rows[(1, 1, 1)][metric] - base[metric]
                shapley.append({"vaccine": vaccine, "metric": metric, "factor": factor,
                                "shapley_delta": value, "joint_delta": joint,
                                "share_of_joint_pct": 100 * value / joint if abs(joint) > 1e-12 else np.nan})
        for d, w, g in BITS:
            on = rows[(d, w, g)]
            off = frame.loc[frame.scenario == f"{PRODUCTS[vaccine][0]}_D{d}W{w}_herd_off"].iloc[0]
            row = {"vaccine": vaccine, "D": d, "W": w, "G": g, "scenario": on.scenario, "paired_herd_off": off.scenario}
            for raw, benefit in (("cases_total_per100k", "cases_averted_per100k"),
                                 ("deaths_total_per100k", "deaths_averted_per100k"),
                                 ("qaly_loss_per100k", "qalys_gained_per100k")):
                contribution = off[raw] - on[raw]
                row[f"indirect_{benefit}"] = contribution
                row[f"direct_only_{benefit}"] = off[benefit]
                row[f"total_{benefit}"] = on[benefit]
                row[f"indirect_share_{benefit}_pct"] = 100 * contribution / on[benefit] if abs(on[benefit]) > 1e-12 else np.nan
            row["indirect_net_cost_savings_per100k"] = off.costs_total_per100k - on.costs_total_per100k
            row["incremental_vaccine_cost_due_to_herd_per100k"] = on.vaccine_costs_per100k - off.vaccine_costs_per100k
            indirect.append(row)
    return pd.DataFrame(interactions), pd.DataFrame(shapley), pd.DataFrame(indirect)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-existing", action="store_true", help="Use persisted summaries; rebuild contrasts only")
    args = ap.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    start = time.time()
    specs = [("no_vaccination", None, 0, 0, 0, True)]
    for vaccine in PRODUCTS:
        specs.extend((scenario_id(vaccine, d, w, g), vaccine, d, w, g, False) for d, w, g in BITS)
        specs.extend((f"{PRODUCTS[vaccine][0]}_D{d}W{w}_herd_off", vaccine, d, w, 0, True) for d, w in itertools.product((0, 1), repeat=2))
    rows, validation, legacy = [], [], []
    if args.skip_existing:
        rows = pd.read_csv(OUTPUT / "raw_summary.csv").to_dict("records")
    else:
        for index, (name, vaccine, d, w, g, off) in enumerate(specs):
            result = simulate(vaccine, d, w, g, name, herd_off=off)
            save_result(name, result)
            annual, costs, details = economics(result, vaccine, name)
            dest = OUTPUT / "runs" / name
            annual.to_csv(dest / "yearly_unrounded.csv", index=False)
            costs.to_csv(dest / "vaccine_costs.csv", index=False)
            details.to_csv(dest / "economic_details.csv", index=False)
            rows.append({"scenario": name, "vaccine": vaccine or "No vaccination", "D": d, "W": w, "G": g,
                         "herd_off": off, "design": "no_vaccination" if vaccine is None else "herd_off" if off else "factorial",
                         **summarize(result, details)})
            for measure in ("population_balance_error", "vaccination_balance_error"):
                error = result["balance"][measure].abs().max()
                validation.append({"check": measure, "scenario": name, "max_abs_error": error, "tolerance": 1e-7, "pass": bool(error < 1e-7)})
            validation.append({"check": "nonnegative_states", "scenario": name, "max_abs_error": result["balance"].minimum_state.min(), "tolerance": 0, "pass": bool(result["balance"].minimum_state.min() >= 0)})
            # Frozen equality is valid only under the historical accounting mode.
            if (d, w, g, off) == (1, 1, 1, False) or vaccine is None:
                if result['parameters'].get('mortality_accounting') == 'historical_inherited':
                    suffix = "baseline" if vaccine is None else f"{PRODUCTS[vaccine][0]}_rollout_to_stable"
                    path = FROZEN / "Refined" / "Outputs" / "Week 6_calendar_pde_program_realism" / "csv" / f"calendar_pde_program_realism_{suffix}_monthly.csv"
                    frozen = pd.read_csv(path)
                    columns = [c for c in frozen if c in result["monthly"]]
                    diff = (frozen[columns] - result["monthly"][columns]).abs().to_numpy().max()
                    validation.append({"check": "historical_week6_monthly_reproduction", "scenario": name, "max_abs_error": diff, "tolerance": 1e-6, "pass": bool(diff < 1e-6)})
                elif vaccine is None:
                    q=1-(1-_Q_NEONATAL)*math.exp(-60*_POSTNEONATAL_HAZARD)
                    validation.append({"check": "NV_under_five_probability_calibration", "scenario": name, "max_abs_error": abs(q-_Q_U5), "tolerance": 1e-14, "pass": bool(abs(q-_Q_U5)<1e-14)})
                _, _, ld = economics(result, vaccine, name, legacy_rounding=True)
                legacy.append({"scenario": name, "vaccine": vaccine or "No vaccination", **summarize(result, ld)})
            print(f"[{index + 1}/{len(specs)}] {name}: {rows[-1]['qaly_loss']:.6f} QALY loss; elapsed={time.time()-start:.1f}s", flush=True)
        pd.DataFrame(rows).to_csv(OUTPUT / "raw_summary.csv", index=False)
        pd.DataFrame(validation).to_csv(OUTPUT / "run_validation.csv", index=False)
        pd.DataFrame(legacy).to_csv(OUTPUT / "legacy_rounded_economics.csv", index=False)
    frame = add_contrasts(pd.DataFrame(rows))
    frame.to_csv(OUTPUT / "all_scenarios_summary.csv", index=False)
    frame.loc[frame.design == "factorial"].to_csv(OUTPUT / "ablation_summary.csv", index=False)
    interactions, shapley, indirect = decompose(frame)
    interactions.to_csv(OUTPUT / "interactions.csv", index=False)
    shapley.to_csv(OUTPUT / "shapley.csv", index=False)
    indirect.to_csv(OUTPUT / "herd_contributions.csv", index=False)
    (OUTPUT / "economic_parameters.json").write_text(json.dumps(ECON_P, indent=2, default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x)), encoding="utf-8")
    (OUTPUT / "run_metadata.json").write_text(json.dumps({"python": sys.version, "platform": platform.platform(),
        "elapsed_seconds": time.time() - start, "timestamp_local": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "randomness": "none; deterministic", "specs": specs, "primary_precision": "unrounded annual aggregation",
        "currency": "2023-aligned USD, inherited 3% price-year proxy", "discount_rate": .03,
        "git": "source is not a git repository; SHA-256 input manifest is authoritative"}, indent=2), encoding="utf-8")
    print("Completed 24 factorial + 12 paired herd-off + 1 no-vaccination runs.")

if __name__ == "__main__":
    main()
