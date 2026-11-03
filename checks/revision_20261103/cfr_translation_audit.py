"""Audit probability-to-hazard and CFR reward semantics on frozen PCV states.

Run from any directory: python checks/revision_20261103/cfr_translation_audit.py
The original 144 trajectories, monthly coefficients and solver are not changed.
This diagnostic reconstructs monthly constant-coefficient exposure from saved
states and saved herd multipliers. CFR-times-episode counts are compared on that
same frozen exposure; this is not a rerun of an alternative mortality model.
All writes are confined to this script's directory and occur only in main().
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd


CORE = ["no_vaccination", "pcv10_gsk_D1W1G1", "pcv13_D1W1G1",
        "pcv10sii_D1W1G1"]
RTOL, ATOL = 1e-12, 1e-7
FLOW_COLUMNS = ["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev",
                "New_Sepsis", "New_AOM", "New_Deaths_Pneumo",
                "New_Seq_Cognitive", "New_Seq_Hearing", "New_Seq_Epilepsy",
                "New_Seq_Hemiplegia", "New_Seq_Hearing_Implant",
                "New_Seq_Hearing_SpEd"]
SOURCES = ["research/inputs/Ablation/src/run_ablation.py"] + [
    "research/inputs/Ablation/inputs/project/Refined/" + name for name in [
        "Week 4_calendar_original.py", "Week 5_calendar_pde_original.py",
        "Week 5_calendar_pde_refined.py", "Week 6_calendar_pde_program_realism.py",
        "Week 7_qaly_refinement_impact.py"]]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def locate_extreme(array, mask, scenario, month, probability, effective, target,
                   kind):
    """Keep the exact nonempty cell that attains a ratio extreme."""
    eligible = np.where(mask, array, np.inf if kind == "min" else -np.inf)
    index = np.unravel_index(
        np.argmin(eligible) if kind == "min" else np.argmax(eligible), array.shape)
    return dict(value=float(array[index]), scenario=scenario, month=int(month),
                dose_state=int(index[0]), age_month=int(index[1]),
                monthly_episode_probability=float(probability[index]),
                effective_CFR=float(effective[index]), nominal_CFR=float(target))


def calculate(repo):
    """Return audit tables and checks without writing files or solving an ODE."""
    sys.path.insert(0, str(repo / "research/inputs/Ablation/src"))
    import run_ablation as r
    w6 = r.w6
    run_root = repo / "research/inputs/Final/outputs/runs"
    directories = sorted(p for p in run_root.iterdir() if p.is_dir())
    source_hashes = {name: sha(repo / name) for name in SOURCES}
    inputs = {p.relative_to(repo).as_posix(): sha(p) for directory in directories
              for name in ["states.npz", "parameters.json", "monthly.csv", "herd.csv"]
              for p in [directory / name]}
    checks = []

    def check(name, value, detail=""):
        checks.append(dict(check=name, passed=bool(value), detail=detail))

    check("expected_144_configurations", len(directories) == 144, len(directories))
    runs, endpoint_rows, reconstruction_rows = [], [], []
    extrema = {scope: {disease: {"min": None, "max": None}
                      for disease in ["meningitis", "severe_pneumonia", "sepsis"]}
               for scope in ["all_144", "four_core"]}
    for run_index, directory in enumerate(directories):
        scenario = directory.name
        p = json.loads((directory / "parameters.json").read_text(encoding="utf-8"))
        monthly = pd.read_csv(directory / "monthly.csv")
        herd = pd.read_csv(directory / "herd.csv")
        with np.load(directory / "states.npz") as archive:
            states = archive["alive_after_births_and_doses"].copy()
        check("saved_state_shape", states.shape == (120, 5, 60), scenario)
        check("saved_month_alignment", len(monthly) == len(herd) == 120 and
              np.array_equal(monthly.Month.to_numpy(), np.arange(120)) and
              np.array_equal(herd.Month.to_numpy(), np.arange(120)), scenario)
        check("finite_nonnegative_saved_states", np.isfinite(states).all() and
              (states >= 0).all(), scenario)
        vaccine = next((v for v, (prefix, _) in r.PRODUCTS.items()
                        if scenario.startswith(prefix + "_")), None)
        reconstructed = {name: [] for name in FLOW_COLUMNS}
        endpoint_totals = {d: dict(cases=0., hazard_deaths=0., cfr_episode_deaths=0.,
                                  positive_exposure_cells=0)
                           for d in ["meningitis", "severe_pneumonia", "sepsis"]}
        local_extrema = {d: {"min": None, "max": None} for d in endpoint_totals}
        men_original, men_alternative, men_shortfall = 0., 0., 0.
        max_coefficient_error = 0.
        for month, A in enumerate(states):
            original_herd = w6.herd_multiplier_for_endpoint
            # Recorded H and H* preserve fixed-herd VE perturbations exactly.
            def saved_herd(pp, vv, endpoint, year, aa):
                return (float(herd.loc[month, "H_" + endpoint]),
                        float(herd.loc[month, "H_target_" + endpoint]),
                        float(herd.loc[month, "relative_carriage_index"]))
            w6.herd_multiplier_for_endpoint = saved_herd
            try:
                c = w6.build_week6_month_coefficients(
                    p, month, A, vaccine is not None, vaccine)
            finally:
                w6.herd_multiplier_for_endpoint = original_herd
            factor = np.ones_like(c["mu_all"])
            np.divide(-np.expm1(-c["mu_all"]), c["mu_all"],
                      out=factor, where=c["mu_all"] != 0)
            E = A * factor
            prop_severe = float(p["prop_pneumo_severe"])
            deaths = np.zeros_like(E)
            month_endpoint = {}
            for disease, hazard_key, f, episode_fraction in [
                    ("meningitis", "lam_men", float(p["cfr_meningitis"]), 1.),
                    ("severe_pneumonia", "lam_pneu", float(p["cfr_pneumo_severe"]), prop_severe),
                    ("sepsis", "lam_sep", float(p["cfr_sepsis"]), 1.)]:
                lam = c[hazard_key]
                probability = -np.expm1(-lam)
                monthly_death_probability = probability * episode_fraction * f
                mu = -np.log1p(-monthly_death_probability)
                deaths += mu
                event_hazard = episode_fraction * lam
                mask = (E > 0) & (probability > 0) & (event_hazard > 0) & (f > 0)
                effective = np.zeros_like(lam)
                np.divide(mu, event_hazard, out=effective, where=event_hazard > 0)
                discrepancy = np.zeros_like(lam)
                if f:
                    discrepancy = effective / f - 1
                cases = float(np.sum(E * event_hazard))
                hazard_deaths = float(np.sum(E * mu))
                cfr_deaths = float(np.sum(E * (f * event_hazard)))
                totals = endpoint_totals[disease]
                totals["cases"] += cases
                totals["hazard_deaths"] += hazard_deaths
                totals["cfr_episode_deaths"] += cfr_deaths
                totals["positive_exposure_cells"] += int(mask.sum())
                month_endpoint[disease] = dict(cases=cases, hazard_deaths=hazard_deaths,
                                               cfr_deaths=cfr_deaths)
                if mask.any():
                    for kind in ["min", "max"]:
                        candidate = locate_extreme(discrepancy, mask, scenario, month,
                                                   probability, effective, f, kind)
                        old = local_extrema[disease][kind]
                        if old is None or (candidate["value"] < old["value"] if kind == "min"
                                           else candidate["value"] > old["value"]):
                            local_extrema[disease][kind] = candidate
            max_coefficient_error = max(max_coefficient_error,
                                        float(np.max(np.abs(deaths - c["mu_dis"]))))
            men_cases = month_endpoint["meningitis"]["cases"]
            original = (1 - float(p["cfr_meningitis"])) * men_cases
            alternative = men_cases - month_endpoint["meningitis"]["hazard_deaths"]
            # Positive quantity: (episodes - hazard deaths) minus legacy source.
            shortfall = (month_endpoint["meningitis"]["cfr_deaths"] -
                         month_endpoint["meningitis"]["hazard_deaths"])
            men_original += original
            men_alternative += alternative
            men_shortfall += shortfall
            flows = dict(New_Meningitis=men_cases,
                         New_Pneumo_Sev=month_endpoint["severe_pneumonia"]["cases"],
                         New_Pneumo_NonSev=float(np.sum(E * c["lam_pneu"])) * (1 - prop_severe),
                         New_Sepsis=month_endpoint["sepsis"]["cases"],
                         New_AOM=float(np.sum(E * c["lam_aom"])),
                         # The original mu_dis is used for this saved-flow check.
                         New_Deaths_Pneumo=float(np.sum(E * c["mu_dis"])))
            flows.update(r.incident_sequela_flows(p, flows))
            for name in FLOW_COLUMNS:
                reconstructed[name].append(flows[name])
        for name in FLOW_COLUMNS:
            actual = np.asarray(reconstructed[name])
            saved = monthly[name].to_numpy(dtype=float)
            difference = actual - saved
            relative = np.zeros_like(difference)
            np.divide(np.abs(difference), np.abs(saved), out=relative, where=saved != 0)
            passed = np.isclose(actual, saved, rtol=RTOL, atol=ATOL)
            worst = int(np.argmax(np.abs(difference)))
            reconstruction_rows.append(dict(
                scenario=scenario, metric=name, months=len(saved),
                max_absolute_error=float(np.max(np.abs(difference))),
                max_relative_error_nonzero_reference=float(relative.max()),
                worst_absolute_month=worst, reconstructed_at_worst=float(actual[worst]),
                saved_at_worst=float(saved[worst]), strict_pass=bool(passed.all()),
                failed_month_count=int((~passed).sum()), rtol=RTOL, atol=ATOL))
            check("frozen_monthly_flow_reconstruction", passed.all(),
                  dict(scenario=scenario, metric=name,
                       max_absolute_error=float(np.max(np.abs(difference))),
                       failed_months=np.flatnonzero(~passed).tolist()))
        check("component_death_hazards_reconstruct_mu_dis", max_coefficient_error < 1e-14,
              dict(scenario=scenario, max_absolute_error=max_coefficient_error,
                   atol=1e-14, note="binary64 conversion consistency, not a clinical tolerance"))
        total_hazard = sum(t["hazard_deaths"] for t in endpoint_totals.values())
        total_cfr = sum(t["cfr_episode_deaths"] for t in endpoint_totals.values())
        saved_total = float(monthly.New_Deaths_Pneumo.sum())
        births = float(monthly.Gross_Births.sum())
        saved_cfr = float((monthly.New_Meningitis * float(p["cfr_meningitis"]) +
                           monthly.New_Pneumo_Sev * float(p["cfr_pneumo_severe"]) +
                           monthly.New_Sepsis * float(p["cfr_sepsis"])).sum())
        runs.append(dict(
            scenario=scenario, core=scenario in CORE, vaccine=vaccine or "No vaccination",
            months=120, recorded_window_gross_births=births,
            saved_frozen_model_deaths=saved_total,
            reconstructed_frozen_hazard_deaths=total_hazard,
            hypothetical_CFR_times_episodes_same_exposure=total_cfr,
            hypothetical_minus_reconstructed_deaths=total_cfr - total_hazard,
            hypothetical_increase_pct_vs_reconstructed=100 * (total_cfr / total_hazard - 1),
            saved_CFR_times_saved_episodes=saved_cfr,
            saved_CFR_times_episodes_minus_saved_deaths=saved_cfr - saved_total,
            saved_reconstruction_residual=total_hazard - saved_total,
            frozen_deaths_per100k_recorded_births=total_hazard * 1e5 / births,
            hypothetical_deaths_per100k_recorded_births=total_cfr * 1e5 / births,
            difference_per100k_recorded_births=(total_cfr - total_hazard) * 1e5 / births,
            meningitis_original_survivor_source=men_original,
            meningitis_episodes_minus_hazard_deaths=men_alternative,
            meningitis_survivor_source_shortfall=men_shortfall,
            meningitis_survivor_source_shortfall_pct=100 * men_shortfall / men_original,
            maximum_component_hazard_conversion_error=max_coefficient_error))
        check("hazard_reward_not_above_CFR_episode_reward", total_cfr >= total_hazard - 1e-10,
              scenario)
        check("meningitis_source_difference_identity",
              abs((men_alternative - men_original) - men_shortfall) < 1e-9, scenario)
        for disease, totals in endpoint_totals.items():
            endpoint_rows.append(dict(
                scenario=scenario, core=scenario in CORE, endpoint=disease, **totals,
                hypothetical_minus_hazard_deaths=totals["cfr_episode_deaths"] - totals["hazard_deaths"],
                hypothetical_increase_pct_vs_hazard=100 *
                (totals["cfr_episode_deaths"] / totals["hazard_deaths"] - 1),
                effective_CFR_ratio_minus_one_min=local_extrema[disease]["min"]["value"],
                effective_CFR_ratio_minus_one_max=local_extrema[disease]["max"]["value"]))
            for scope in ["all_144"] + (["four_core"] if scenario in CORE else []):
                for kind in ["min", "max"]:
                    candidate = local_extrema[disease][kind]
                    old = extrema[scope][disease][kind]
                    if old is None or (candidate["value"] < old["value"] if kind == "min"
                                       else candidate["value"] > old["value"]):
                        extrema[scope][disease][kind] = candidate
        if (run_index + 1) % 12 == 0:
            print(f"CFR translation audit: {run_index + 1}/{len(directories)} configurations", flush=True)
    check("input_files_unchanged", all(sha(repo / name) == digest for name, digest in inputs.items()))
    check("model_sources_unchanged", all(sha(repo / name) == digest for name, digest in source_hashes.items()))
    frames = dict(runs=pd.DataFrame(runs), endpoints=pd.DataFrame(endpoint_rows),
                  reconstruction=pd.DataFrame(reconstruction_rows))
    frames["core"] = frames["runs"].set_index("scenario").loc[CORE].reset_index()
    extreme_rows = [dict(scope=scope, endpoint=disease, extreme=kind, **record)
                    for scope, diseases in extrema.items() for disease, kinds in diseases.items()
                    for kind, record in kinds.items()]
    frames["extrema"] = pd.DataFrame(extreme_rows)
    recon = frames["reconstruction"]
    groups = {}
    for scope, frame in [("all_144", frames["runs"]), ("four_core", frames["core"])]:
        groups[scope] = dict(
            configurations=len(frame),
            hypothetical_increase_pct_range=[float(frame.hypothetical_increase_pct_vs_reconstructed.min()),
                                             float(frame.hypothetical_increase_pct_vs_reconstructed.max())],
            difference_per100k_recorded_births_range=[float(frame.difference_per100k_recorded_births.min()),
                                                      float(frame.difference_per100k_recorded_births.max())],
            meningitis_survivor_source_shortfall_range=[float(frame.meningitis_survivor_source_shortfall.min()),
                                                       float(frame.meningitis_survivor_source_shortfall.max())])
    summary = dict(
        audit="probability_to_hazard_CFR_translation", date="2026-11-03",
        scope="144 saved configurations x 120 months x five dose states x 60 attained-age bins; no ODE rerun",
        interpretation="Fixed-trajectory exposure/reward-definition diagnostic. Hypothetical CFR-times-episode rewards use the original frozen exposures; they are not the outcomes of a fully rerun alternative disease process and do not alter published results.",
        formulas=dict(episode_hazard="lambda=-log(1-p)",
                      death_hazard="mu=-log(1-p*f)",
                      exposure="E=A*(1-exp(-mu_all))/mu_all, with continuous limit A at zero",
                      ordinary_effective_CFR="mu/lambda",
                      severe_pneumonia_effective_CFR="mu/(prop_severe*lambda_pneumonia), where f=prop_severe*CFR_severe in the death-probability conversion",
                      nominal_meningitis_source="(1-CFR_meningitis)*lambda_meningitis*E",
                      source_shortfall="sum(E*(CFR_meningitis*lambda_meningitis-mu_meningitis)); alternative (episodes-hazard deaths) minus original source"),
        audit_numerics="expm1/log1p are used only inside this diagnostic; original coefficient and solver implementations are unchanged",
        probability_recovery="p=-expm1(-lambda) inverts the saved model coefficient to binary64 precision; component mu is reconstructed with -log1p, and its sum is explicitly checked against original mu_dis",
        reconstruction=dict(rtol=RTOL, atol=ATOL,
                            max_absolute_monthly_flow_error=float(recon.max_absolute_error.max()),
                            max_relative_monthly_flow_error_nonzero_reference=float(recon.max_relative_error_nonzero_reference.max()),
                            maximum_component_hazard_conversion_error=float(frames["runs"].maximum_component_hazard_conversion_error.max()),
                            max_absolute_cumulative_component_deaths_minus_saved=float(frames["runs"].saved_reconstruction_residual.abs().max()),
                            failed_scenario_metric_checks=int((~recon.strict_pass).sum()),
                            failed_monthly_flow_comparisons=int(recon.failed_month_count.sum())),
        groups=groups, extrema=extrema,
        core=frames["core"].to_dict(orient="records"),
        checks=len(checks), passed=sum(x["passed"] for x in checks),
        failed=[x for x in checks if not x["passed"]],
        source_sha256=source_hashes, input_sha256=inputs,
        historical_811_checks_rerun=False,
        cautions=["The 144 configurations include repeated valuation-only trajectories and are not independent empirical samples.",
                  "No joint disease-course state model, calibrated clinical CFR estimate, stochastic PCV risk model or parameter PSA is introduced.",
                  "Changes to death hazards would change survival and feedback in a complete rerun; this fixed-exposure diagnostic does not estimate that change.",
                  "Sequela incidence sources use nominal CFR; downstream stocks and economics are left unchanged."])
    return frames, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path)
    args = parser.parse_args()
    repo = args.repo.resolve() if args.repo else next(
        p for p in Path(__file__).resolve().parents if (p / "reproduce.py").is_file())
    frames, summary = calculate(repo)
    out = Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in frames.items():
        frame.to_csv(out / ("cfr_translation_" + name + ".csv"), index=False)
    summary["audit_script_sha256"] = sha(Path(__file__))
    (out / "cfr_translation_audit.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ["reconstruction", "groups", "checks", "passed", "failed"]}, indent=2), flush=True)
    if summary["failed"]:
        raise SystemExit("CFR audit contains failed checks; inspect the JSON and CSV evidence.")


if __name__ == "__main__":
    main()
