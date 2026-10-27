import copy
import math
import os
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================================================
# Flexible imports: support either *_new.py or original names
# ============================================================
try:
    from markov_new import (
        VietnamRandSimulation,
        summarize_micro,
        params_vietnam_micro as markov_params,
        _make_immediate_irrs_from_year5 as markov_make_immediate_irrs_from_year5,
    )
except ImportError:
    from markov import (
        VietnamRandSimulation,
        summarize_micro,
        params_vietnam_micro as markov_params,
        _make_immediate_irrs_from_year5 as markov_make_immediate_irrs_from_year5,
    )

try:
    from ode_new import (
        run_ode_and_build_monthly,
        aggregate_yearly_like_micro,
        params_vietnam_micro as ode_params,
        _make_immediate_irrs_from_year5 as ode_make_immediate_irrs_from_year5,
    )
except ImportError:
    from ode import (
        run_ode_and_build_monthly,
        aggregate_yearly_like_micro,
        params_vietnam_micro as ode_params,
        _make_immediate_irrs_from_year5 as ode_make_immediate_irrs_from_year5,
    )

# ===================== Configuration =====================
N_MARKOV_RUNS = 20
SEED_BASE: Optional[int] = None  # e.g. 20260222 for reproducibility

OUTPUT_CSV_DIR = "combined_csv_outputs"
OUTPUT_FIG_DIR = "combined_figure_outputs"
OUTPUT_RAW_MARKOV_DIR = "markov_raw_yearly_outputs"
EXPORT_ALL_MARKOV_RUNS_XLSX = True

# Beeswarm plotting controls
SCATTER_MAX_HALF_WIDTH = 0.34
SCATTER_LANE_STEP = 0.015
SCATTER_ALPHA = 0.18
SCATTER_SIZE = 8
DRAW_VIOLIN_BACKGROUND = False  # optional visual layer

SCENARIOS = [
    ("baseline", False, None),
    ("pcv10_dynamic", True, "PCV10-GSK"),
    ("pcv13_dynamic", True, "PCV13"),
    ("pcv10_year5_immediate", True, "PCV10-GSK"),
    ("pcv13_year5_immediate", True, "PCV13"),
    ("pcv10sii_year5_immediate", True, "PCV10-SII"),
]

METRICS = [
    "New_Meningitis",
    "New_Pneumo_Sev",
    "New_Pneumo_NonSev",
    "New_Sepsis",
    "New_AOM",
    "New_Deaths_Pneumo",
    "Stock_Alive",
    "Stock_Dead",
    "Stock_New_Meningitis",
    "Stock_New_Pneumo_Sev",
    "Stock_New_Pneumo_NonSev",
    "Stock_New_Sepsis",
    "Stock_New_AOM",
    "Stock_New_Deaths_Pneumo",
    "Stock_Seq_Cognitive",
    "Stock_Seq_Hearing",
    "Stock_Seq_Epilepsy",
    "Stock_Seq_Hemiplegia",
    "Stock_Seq_Hearing_Implant",
    "Stock_Seq_Hearing_SpEd",
]
COL_ORDER = ["Year"] + METRICS

for _d in [OUTPUT_CSV_DIR, OUTPUT_FIG_DIR, OUTPUT_RAW_MARKOV_DIR]:
    os.makedirs(_d, exist_ok=True)

plt.rcParams["font.size"] = 10
plt.rcParams["figure.figsize"] = (12, 7)


# ===================== Parameter helpers =====================
def apply_year5_immediate_irr_ode(params: dict, scenario_name: str, vaccine_type: Optional[str]) -> dict:
    if ("year5_immediate" not in scenario_name) or (vaccine_type is None):
        return params
    irr_patch = ode_make_immediate_irrs_from_year5(vaccine_type)
    base = dict(params.get("vaccine_irrs", {}))
    base.update(irr_patch)
    params["vaccine_irrs"] = base
    return params


def apply_year5_immediate_irr_markov(params: dict, scenario_name: str, vaccine_type: Optional[str]) -> dict:
    if ("year5_immediate" not in scenario_name) or (vaccine_type is None):
        return params
    irr_patch = markov_make_immediate_irrs_from_year5(vaccine_type)
    base = dict(params.get("vaccine_irrs", {}))
    base.update(irr_patch)
    params["vaccine_irrs"] = base
    return params


# ===================== Validation =====================
def _validate_yearly_df(df: pd.DataFrame, label: str) -> pd.DataFrame:
    out = df.copy()
    if "Year" not in out.columns:
        raise ValueError(f"{label}: missing 'Year' column")
    for c in COL_ORDER:
        if c not in out.columns:
            out[c] = np.nan
    out["Year"] = out["Year"].astype(int)
    out = out.reindex(columns=COL_ORDER)
    out = out.sort_values("Year").reset_index(drop=True)
    return out


# ===================== ODE =====================
def run_ode_yearly(scenario_name: str, intervention: bool, vaccine_type: Optional[str]) -> pd.DataFrame:
    params = copy.deepcopy(ode_params)
    params = apply_year5_immediate_irr_ode(params, scenario_name, vaccine_type)

    df_monthly = run_ode_and_build_monthly(params, intervention=intervention, vaccine_type=vaccine_type)
    df_yearly = aggregate_yearly_like_micro(df_monthly)
    return _validate_yearly_df(df_yearly, f"ODE-{scenario_name}")


# ===================== Markov multi-run =====================
def run_markov_yearly_multi(
    scenario_name: str,
    intervention: bool,
    vaccine_type: Optional[str],
) -> Tuple[pd.DataFrame, List[pd.DataFrame]]:
    params = copy.deepcopy(markov_params)
    params = apply_year5_immediate_irr_markov(params, scenario_name, vaccine_type)

    raw_runs: List[pd.DataFrame] = []
    failures = 0

    for run_idx in range(N_MARKOV_RUNS):
        if SEED_BASE is not None:
            np.random.seed(SEED_BASE + run_idx)

        if run_idx == 0 or (run_idx + 1) % 25 == 0 or (run_idx + 1) == N_MARKOV_RUNS:
            print(f"Markov {scenario_name}: run {run_idx + 1}/{N_MARKOV_RUNS}")

        try:
            model = VietnamRandSimulation(params)
            df_monthly = model.run(intervention=intervention, vaccine_type=vaccine_type)
            df_yearly = summarize_micro(df_monthly)
            df_yearly = _validate_yearly_df(df_yearly, f"Markov-{scenario_name}-run{run_idx+1}")
            raw_runs.append(df_yearly)
        except Exception as e:
            failures += 1
            print(f"[WARN] {scenario_name} run {run_idx + 1} failed: {e}")

    if len(raw_runs) == 0:
        raise RuntimeError(f"All Markov runs failed for scenario: {scenario_name}")
    if failures > 0:
        print(f"[WARN] {scenario_name}: {failures} runs failed, stats use successful runs only")

    all_runs = pd.concat(
        raw_runs,
        keys=[f"run_{i+1}" for i in range(len(raw_runs))],
        names=["Run", "Row"],
    )

    years = sorted(all_runs["Year"].unique().tolist())
    stats_df = pd.DataFrame({"Year": years})
    stats_df["Year"] = stats_df["Year"].astype(int)

    grouped = all_runs.groupby("Year", as_index=False)
    for m in METRICS:
        gm = grouped[m]
        stats_df[f"{m}_median"] = gm.median()[m].values
        stats_df[f"{m}_mean"] = gm.mean()[m].values
        stats_df[f"{m}_max"] = gm.max()[m].values
        stats_df[f"{m}_min"] = gm.min()[m].values
        stats_df[f"{m}_std"] = np.round(gm.std(ddof=0)[m].astype(float).values, 2)

    return stats_df, raw_runs


# ===================== Raw Markov exports =====================
def raw_runs_to_long_df(raw_runs: List[pd.DataFrame], scenario_name: str) -> pd.DataFrame:
    """
    Convert list of yearly run DataFrames into one long DataFrame.
    Columns: Scenario, Run, Year, <metrics...>
    """
    pieces = []
    for i, df in enumerate(raw_runs, start=1):
        rr = _validate_yearly_df(df, f"raw-export-{scenario_name}-run{i}")
        rr = rr.copy()
        rr.insert(0, "Run", i)
        rr.insert(0, "Scenario", scenario_name)
        pieces.append(rr)

    out = pd.concat(pieces, ignore_index=True)
    out = out[["Scenario", "Run", "Year"] + METRICS]
    return out


# ===================== Combined numeric table =====================
def combine_ode_markov_numeric(ode_df: pd.DataFrame, markov_stats_df: pd.DataFrame) -> pd.DataFrame:
    ode_df = _validate_yearly_df(ode_df, "ODE yearly")
    mk_stats = markov_stats_df.copy()
    if "Year" not in mk_stats.columns:
        raise ValueError("Markov stats missing Year")
    mk_stats["Year"] = mk_stats["Year"].astype(int)

    ode_renamed = ode_df[["Year"] + METRICS].rename(columns={m: f"ODE_{m}" for m in METRICS})

    required = ["Year"]
    for m in METRICS:
        for stat in ["median", "mean", "max", "min", "std"]:
            required.append(f"{m}_{stat}")
    for c in required:
        if c not in mk_stats.columns:
            mk_stats[c] = np.nan

    mk = mk_stats[required].copy().rename(
        columns={f"{m}_{stat}": f"Markov_{stat}_{m}" for m in METRICS for stat in ["median", "mean", "max", "min", "std"]}
    )

    combined = pd.merge(ode_renamed, mk, on="Year", how="left")

    ordered = ["Year"]
    for m in METRICS:
        ordered.append(f"ODE_{m}")
        for stat in ["median", "mean", "max", "min", "std"]:
            ordered.append(f"Markov_{stat}_{m}")
    for c in ordered:
        if c not in combined.columns:
            combined[c] = np.nan

    return combined[ordered].sort_values("Year").reset_index(drop=True)


# ===================== Compact cell table =====================
def _fmt_1(x: float) -> str:
    if pd.isna(x):
        return "NA"
    return f"{float(x):.1f}"


def _fmt_2(x: float) -> str:
    if pd.isna(x):
        return "NA"
    return f"{float(x):.2f}"


def _round_half_up_to_int_float(x: float) -> float:
    xv = float(x)
    if xv >= 0:
        return float(math.floor(xv + 0.5))
    return float(math.ceil(xv - 0.5))


def _fmt_markov_mean_display(x: float) -> str:
    if pd.isna(x):
        return "NA"
    return f"{_round_half_up_to_int_float(float(x)):.1f}"


def combine_ode_markov_compact_cells(ode_df: pd.DataFrame, markov_stats_df: pd.DataFrame) -> pd.DataFrame:
    num_df = combine_ode_markov_numeric(ode_df, markov_stats_df)
    rows = []
    for _, r in num_df.iterrows():
        row = {"Year": int(r["Year"])}
        for m in METRICS:
            row[m] = " | ".join(
                [
                    _fmt_1(r.get(f"ODE_{m}", np.nan)),
                    _fmt_1(r.get(f"Markov_median_{m}", np.nan)),
                    _fmt_markov_mean_display(r.get(f"Markov_mean_{m}", np.nan)),
                    _fmt_1(r.get(f"Markov_max_{m}", np.nan)),
                    _fmt_1(r.get(f"Markov_min_{m}", np.nan)),
                    _fmt_2(r.get(f"Markov_std_{m}", np.nan)),
                ]
            )
        rows.append(row)

    return pd.DataFrame(rows, columns=["Year"] + METRICS)


# ===================== Beeswarm spreading =====================
def _deterministic_unit(seed_value: float) -> float:
    """Deterministic pseudo-random in [0,1] based on a float input."""
    return 0.5 * (1.0 + math.sin(seed_value * 12.9898 + 78.233))


def _center_out_offsets(n: int, step: float) -> np.ndarray:
    """
    Generate offsets like [0, +1, -1, +2, -2, ...] * step.
    This keeps low multiplicity points near center.
    """
    if n <= 0:
        return np.array([], dtype=float)
    if n == 1:
        return np.array([0.0], dtype=float)

    offsets = [0.0]
    k = 1
    while len(offsets) < n:
        offsets.append(k * step)
        if len(offsets) < n:
            offsets.append(-k * step)
        k += 1
    return np.array(offsets[:n], dtype=float)


def _spread_duplicate_y_points(
    x_center: float,
    y_values: List[float],
    max_half_width: float = SCATTER_MAX_HALF_WIDTH,
    step: float = SCATTER_LANE_STEP,
) -> Tuple[List[float], List[float]]:
    """
    Beeswarm-style spreading for duplicated y-values at one x position.
    Desired behavior:
    - one point -> center
    - repeated equal values -> gradually spread left/right
    - avoid hard side walls that look like two vertical lines
    """
    if len(y_values) == 0:
        return [], []

    counts: Dict[float, int] = {}
    for y in y_values:
        key = float(round(float(y), 10))
        counts[key] = counts.get(key, 0) + 1

    xs: List[float] = []
    ys: List[float] = []

    # Process by y so the visual stack is stable
    for y_key in sorted(counts.keys()):
        n = counts[y_key]
        offsets = _center_out_offsets(n, step=step)

        if len(offsets) > 0:
            half_span = float(np.max(np.abs(offsets)))
            if half_span > max_half_width and half_span > 0:
                # Softly vary target width across y-levels to avoid repeated edge walls
                u = _deterministic_unit(y_key + 0.001 * n)
                target_half = max_half_width * (0.82 + 0.16 * u)
                offsets = offsets * (target_half / half_span)

            # Tiny deterministic micro-jitter in x so dense stacks look more cloud-like
            # It is small enough to preserve center-out structure.
            micro_amp = min(0.004, step * 0.22)
            for j in range(n):
                u2 = _deterministic_unit(y_key * 1.7 + (j + 1) * 0.37)
                offsets[j] += (u2 - 0.5) * 2.0 * micro_amp

        for off in offsets:
            xs.append(float(x_center) + float(off))
            ys.append(float(y_key))

    return xs, ys


# ===================== Plotting =====================
def plot_metric_trend_with_scatter(
    scenario_name: str,
    ode_df: pd.DataFrame,
    markov_stats_df: pd.DataFrame,
    markov_raw_runs: List[pd.DataFrame],
    metric: str,
) -> None:
    ode_df = _validate_yearly_df(ode_df, "ODE yearly (plot)")
    years = ode_df["Year"].astype(int).values
    ode_vals = ode_df[metric].astype(float).values

    # stats by year
    stats_map: Dict[int, Dict[str, float]] = {}
    for _, r in markov_stats_df.iterrows():
        y = int(r["Year"])
        stats_map[y] = {
            "median": float(r.get(f"{metric}_median", np.nan)),
            "mean": float(r.get(f"{metric}_mean", np.nan)),
            "std": float(r.get(f"{metric}_std", np.nan)),
        }

    markov_medians = np.array([stats_map.get(int(y), {}).get("median", np.nan) for y in years], dtype=float)
    markov_means = np.array([stats_map.get(int(y), {}).get("mean", np.nan) for y in years], dtype=float)
    markov_stds = np.array([stats_map.get(int(y), {}).get("std", np.nan) for y in years], dtype=float)

    # raw values by year
    year_to_vals: Dict[int, List[float]] = {int(y): [] for y in years}
    for run_df in markov_raw_runs:
        rr = _validate_yearly_df(run_df, "Markov run (plot)")
        for y in years:
            row = rr.loc[rr["Year"] == int(y), metric]
            if len(row) == 0:
                continue
            year_to_vals[int(y)].append(float(row.iloc[0]))

    scatter_x: List[float] = []
    scatter_y: List[float] = []
    for y in years:
        xs, ys = _spread_duplicate_y_points(
            x_center=float(y),
            y_values=year_to_vals.get(int(y), []),
            max_half_width=SCATTER_MAX_HALF_WIDTH,
            step=SCATTER_LANE_STEP,
        )
        scatter_x.extend(xs)
        scatter_y.extend(ys)

    fig, ax = plt.subplots()

    if DRAW_VIOLIN_BACKGROUND:
        violin_data = [year_to_vals.get(int(y), []) for y in years]
        parts = ax.violinplot(
            violin_data,
            positions=years.astype(float),
            widths=0.55,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for pc in parts["bodies"]:
            pc.set_alpha(0.10)

    ax.scatter(scatter_x, scatter_y, alpha=SCATTER_ALPHA, s=SCATTER_SIZE, label="Markov Raw Runs", zorder=1)
    ax.plot(years, ode_vals, "o-", linewidth=2.5, markersize=8, label="ODE Value", zorder=4)
    ax.plot(years, markov_means, "^--", linewidth=2.0, markersize=6, label="Markov Mean", zorder=2)
    ax.errorbar(
        years,
        markov_medians,
        yerr=markov_stds,
        fmt="s-",
        capsize=6,
        linewidth=2.0,
        markersize=7,
        label="Markov Median (±SD)",
        zorder=3,
    )

    ax.set_xlabel("Year", fontsize=14, fontweight="bold")
    ax.set_ylabel(metric, fontsize=14, fontweight="bold")
    ax.set_title(f"{metric} - Scenario: {scenario_name}", fontsize=16, fontweight="bold")
    ax.set_xticks(years)
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="best", fontsize=11)

    fig_path = os.path.join(OUTPUT_FIG_DIR, f"{scenario_name}_{metric}.png")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ===================== Raw Markov collection export =====================
def save_all_markov_runs_collection(scenario_to_raw_df: Dict[str, pd.DataFrame]) -> None:
    """
    Save all Markov yearly runs into a single CSV and an optional XLSX workbook.
    - CSV long format: all scenarios stacked
    - XLSX: one sheet per scenario + all_scenarios sheet
    """
    if len(scenario_to_raw_df) == 0:
        return

    all_long = pd.concat([scenario_to_raw_df[s] for s in scenario_to_raw_df], ignore_index=True)
    all_long = all_long[["Scenario", "Run", "Year"] + METRICS]

    csv_path = os.path.join(OUTPUT_RAW_MARKOV_DIR, f"all_markov_yearly_runs_{N_MARKOV_RUNS}.csv")
    all_long.to_csv(csv_path, index=False, encoding="utf-8")
    print(f"Saved raw Markov collection CSV: {csv_path}")

    if EXPORT_ALL_MARKOV_RUNS_XLSX:
        xlsx_path = os.path.join(OUTPUT_RAW_MARKOV_DIR, f"all_markov_yearly_runs_{N_MARKOV_RUNS}.xlsx")
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            all_long.to_excel(writer, sheet_name="all_scenarios", index=False)
            for scenario_name, df in scenario_to_raw_df.items():
                # Excel sheet names must be <= 31 chars
                sheet = scenario_name[:31]
                df.to_excel(writer, sheet_name=sheet, index=False)
        print(f"Saved raw Markov collection XLSX: {xlsx_path}")


# ===================== Main =====================
if __name__ == "__main__":
    print("===== Running Combined Workflow (ODE + Markov) =====")
    print(f"N_MARKOV_RUNS = {N_MARKOV_RUNS}")
    print(f"SEED_BASE = {SEED_BASE}")
    print(f"DRAW_VIOLIN_BACKGROUND = {DRAW_VIOLIN_BACKGROUND}")

    scenario_to_raw_runs_long: Dict[str, pd.DataFrame] = {}

    for scenario_name, intervention, vaccine_type in SCENARIOS:
        print("\n" + "=" * 60)
        print(f"Scenario: {scenario_name}")
        print("=" * 60)

        print("Running ODE...")
        ode_df = run_ode_yearly(scenario_name, intervention, vaccine_type)

        print("Running Markov Monte Carlo...")
        markov_stats_df, markov_raw_runs = run_markov_yearly_multi(scenario_name, intervention, vaccine_type)

        # Save per-scenario raw 1000-run yearly table
        raw_long_df = raw_runs_to_long_df(markov_raw_runs, scenario_name)
        scenario_to_raw_runs_long[scenario_name] = raw_long_df
        raw_csv_path = os.path.join(OUTPUT_RAW_MARKOV_DIR, f"markov_raw_yearly_{scenario_name}.csv")
        raw_long_df.to_csv(raw_csv_path, index=False, encoding="utf-8")
        print(f"Saved scenario raw Markov yearly CSV: {raw_csv_path}")

        # Save combined tables
        print("Saving combined tables...")
        combined_numeric_df = combine_ode_markov_numeric(ode_df, markov_stats_df)
        combined_compact_df = combine_ode_markov_compact_cells(ode_df, markov_stats_df)

        numeric_csv = os.path.join(OUTPUT_CSV_DIR, f"combined_{scenario_name}_numeric.csv")
        compact_csv = os.path.join(OUTPUT_CSV_DIR, f"combined_{scenario_name}.csv")

        combined_numeric_df.to_csv(numeric_csv, index=False, encoding="utf-8")
        combined_compact_df.to_csv(compact_csv, index=False, encoding="utf-8")

        print(f"Saved combined numeric CSV: {numeric_csv}")
        print(f"Saved combined compact CSV: {compact_csv}")

        # Plots
        print("Generating plots...")
        for metric in METRICS:
            print(f"  Plotting {metric}")
            plot_metric_trend_with_scatter(
                scenario_name=scenario_name,
                ode_df=ode_df,
                markov_stats_df=markov_stats_df,
                markov_raw_runs=markov_raw_runs,
                metric=metric,
            )

    # Save one master document collecting all 1000-run Markov results
    print("\nSaving master raw Markov collection document...")
    save_all_markov_runs_collection(scenario_to_raw_runs_long)

    print("\n===== All tasks completed =====")
    print(f"Combined table CSV dir: {os.path.abspath(OUTPUT_CSV_DIR)}")
    print(f"Plot dir: {os.path.abspath(OUTPUT_FIG_DIR)}")
    print(f"Raw Markov collection dir: {os.path.abspath(OUTPUT_RAW_MARKOV_DIR)}")
