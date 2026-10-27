import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =========================
# User settings
# =========================
CONFIG = {
    # Cohort / horizon
    "BASE_COHORT_YEAR": 2023,
    "HORIZON_YEARS": 5,
    "LIFE_EXPECTANCY_YEARS": 74.0,
    "DISCOUNT_RATE": 0.03,

    # Price-year alignment proxy
    # We align known 2024/2025 monetary inputs back to the 2023 birth cohort year
    # using 3% because no CPI/deflator table was provided.
    "PRICE_YEAR_ALIGNMENT_RATE": 0.03,

    # Vaccination program assumptions not explicitly provided in the parameter table
    "DOSES_PER_CHILD": 3,
    "VACCINE_EVENT_AGE_YEARS": 0.25,

    # Baseline strategy vaccine cost
    # None means baseline has no vaccine-program cost.
    # If your baseline represents status quo private market vaccination, set this to "pcv10" or "pcv13".
    "BASELINE_VACCINE_PROGRAM": None,

    # Input files (script also tries *_UPDATED.csv automatically)
    "FILES": [
        "combined_baseline.csv",
        "combined_pcv10_dynamic.csv",
        "combined_pcv10_year5_immediate.csv",
        "combined_pcv10sii_year5_immediate.csv",
        "combined_pcv13_dynamic.csv",
        "combined_pcv13_year5_immediate.csv",
    ],

    # Current combine format:
    # ODE | Markov_median | Markov_mean | Markov_max | Markov_min | Markov_sd
    "MODEL_INDEX": {"ODE": 0, "Markov": 1},

    # LOS values kept from your previous code
    "LOS_DAYS": {
        "inpatient_pneumonia": 8.74,
        "outpatient_pneumonia": 7.0,
        "meningitis": 18.41,
        "sepsis": 12.44,
        "aom": 7.0,
    },

    # Average age at death by model year (your requested setting)
    "AVG_DEATH_AGE_BY_YEAR": {
        1: 0.5,
        2: 1.5,
        3: 2.5,
        4: 3.5,
        5: 4.5,
    },
}


EXPECTED_COLS = [
    "Year",
    "New_Meningitis",
    "New_Pneumo_Sev",
    "New_Pneumo_NonSev",
    "New_Sepsis",
    "New_AOM",
    "New_Deaths_Pneumo",
    "Stock_Seq_Cognitive",
    "Stock_Seq_Hearing",
    "Stock_Seq_Epilepsy",
    "Stock_Seq_Hemiplegia",
    "Stock_Seq_Hearing_Implant",
    "Stock_Seq_Hearing_SpEd",
]


# =========================
# Small utilities
# =========================
def _norm_col(s: str) -> str:
    return "".join(ch for ch in str(s).lower() if ch.isalnum())


def parse_numeric_any(x) -> float:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        raise ValueError("Empty numeric cell")
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip().replace(",", "")
    m = re.search(r"[-+]?\d*\.?\d+", s)
    if not m:
        raise ValueError(f"Cannot parse numeric value from: {x}")
    return float(m.group(0))


def parse_range_midpoint(x) -> float:
    if x is None:
        raise ValueError("Empty range cell")
    s = str(x).strip().replace(",", "")
    # Cost ranges are written as positive hyphenated intervals, e.g. "680-1,360".
    # After comma removal the hyphen is a separator, not a negative sign.
    nums = re.findall(r"\d+(?:\.\d+)?", s)
    if len(nums) < 2:
        return parse_numeric_any(x)
    a, b = float(nums[0]), float(nums[1])
    return (a + b) / 2.0


def parse_cell_value(x, idx: int) -> float:
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    try:
        if "|" in s:
            parts = [p.strip() for p in s.split("|") if p.strip() != ""]
            if not parts:
                return np.nan
            if idx < len(parts):
                return float(parts[idx].replace(",", ""))
            return float(parts[-1].replace(",", ""))
        return float(s.replace(",", ""))
    except Exception:
        return np.nan


def resolve_file(path: str) -> str:
    if os.path.exists(path):
        return path
    if path.lower().endswith(".csv"):
        alt = path[:-4] + "_UPDATED.csv"
        if os.path.exists(alt):
            return alt
    raise FileNotFoundError(f"Input file not found: {path}")


def load_and_extract(filename: str, model_idx: int, expected_cols: List[str]) -> pd.DataFrame:
    path = resolve_file(filename)
    raw = pd.read_csv(path, dtype=str)
    mapping: Dict[str, Optional[str]] = {}

    for exp in expected_cols:
        expn = _norm_col(exp)
        found = None
        for c in raw.columns:
            if _norm_col(c) == expn:
                found = c
                break
        if found is None:
            for c in raw.columns:
                if exp.lower() in str(c).lower():
                    found = c
                    break
        mapping[exp] = found

    out = pd.DataFrame(index=range(len(raw)))
    for exp, real_col in mapping.items():
        if real_col is None:
            out[exp] = np.nan
        else:
            out[exp] = raw[real_col].apply(lambda v: parse_cell_value(v, model_idx))

    required = ["New_Meningitis", "New_Pneumo_Sev", "New_Pneumo_NonSev", "New_Sepsis", "New_AOM", "New_Deaths_Pneumo"]
    missing = [c for c in required if c not in out.columns or out[c].isna().all()]
    if missing:
        raise ValueError(f"{filename}: missing required columns or all-NaN after parsing: {missing}")

    return out


# =========================
# Hardcoded parameters from your Excel table
# =========================
RAW_PARAM_TABLE = {
    # Core model settings from the parameter table
    "birth_cohort_2023": "newborns=1384610",
    "discount_rate": 0.03,
    "coverage_pcv13_epi": 0.9,
    "coverage_pcv10_epi": 0.9,

    # Direct medical costs (table section labeled 2024)
    "cost_medical_2024": {
        "inpatient_pneumonia": "mean=221.09",
        "outpatient_pneumonia": "mean=16.09",
        "inpatient_meningitis": "mean=574.50",
        "inpatient_sepsis": "mean=589.60",
        "outpatient_aom": "mean=13.31",
        "cognitive_annual": "mean=73.36",
        "hearing_case_range": "680-1,360",
        "epilepsy_case": 329.06,
        "hemiplegia_case_range": "606.61-1390.34",
        "cochlear_implant_case_range": "20691.18-28967.65",
    },

    # Direct non-medical costs (table section labeled 2025)
    "cost_nonmedical_2025": {
        "inpatient_pneumonia": "mean=142.33",
        "outpatient_pneumonia": "mean=18.98",
        "inpatient_meningitis": "mean=297.95",
        "inpatient_sepsis": "mean=305.24",
        "outpatient_aom": "mean=21.21",
    },

    # Caregiver absenteeism costs (table section labeled 2025)
    "cost_caregiver_2025": {
        "inpatient_pneumonia": "mean=184.25",
        "outpatient_pneumonia": "mean=124.49",
        "inpatient_meningitis": "mean=297.82",
        "inpatient_sepsis": "mean=193.55",
        "outpatient_aom": "mean=139.12",
    },

    # Special education annual cost (same 2025 section)
    "special_education_annual_range_2025": "1638.74-3421.49",

    # Vaccine / program costs
    "vaccine_price_pcv10_2024": "2 (median price)",
    "vaccine_price_pcv13_unknown_year": 4.22,
    "service_delivery_per_dose_2024": 2.1,
    "severe_ae_cost_range_unknown_year": "4717.59-4883.12",
    "severe_ae_incidence_per_100k_doses": 4.01,
    "wastage_rate": 0.05,

    # Utilities
    "utilities": {
        "inpatient_pneumonia": 0.9941,
        "outpatient_pneumonia": 0.9963,
        "meningitis": 0.9768,
        "sepsis": 0.9921,
        "aom": 0.9950,
        "cognitive": 0.62,
        "hearing": 0.91,
        "epilepsy": 0.83,
        "hemiplegia": 0.3903,
    },

    # Productivity discount sheet, column "Discounted productivity loss per capita"
    # Ages 15-60 only, exactly as requested.
    "prod_discounted_age_contrib_15_60": {
        15: 1589.287938641466781,
        16: 1542.997998681036051,
        17: 1498.056309399064048,
        18: 1454.423601358314500,
        19: 1412.061748891567504,
        20: 1370.933736787929774,
        21: 1331.003627949446582,
        22: 1292.236531989753757,
        23: 1254.598574747333714,
        24: 1218.056868686731832,
        25: 1263.746972649093777,
        26: 1226.938808397178263,
        27: 1191.202726599202379,
        28: 1156.507501552623808,
        29: 1122.822817041382450,
        30: 1090.119239846002301,
        31: 1058.368194025244748,
        32: 1027.541935946839658,
        33: 997.613530045475500,
        34: 968.556825286869525,
        35: 940.346432317348899,
        36: 912.957701278979584,
        37: 886.366700270854039,
        38: 860.550194437722439,
        39: 835.485625667691579,
        40: 811.151092881254044,
        41: 787.525332894421354,
        42: 764.587701839244005,
        43: 742.318157125479615,
        44: 720.697239927650230,
        45: 699.706058182184734,
        46: 679.326270079790902,
        47: 659.540068038631944,
        48: 640.330163144302901,
        49: 621.679770043012581,
        50: 604.738027158226828,
        51: 587.124298211870723,
        52: 570.023590496961901,
        53: 553.420961647535933,
        54: 537.301904512170722,
        55: 521.652334477835552,
        56: 506.458577162947279,
        57: 491.707356468880789,
        58: 477.385782979495900,
        59: 463.481342698539720,
        60: 449.981886115087150
    },
}


def build_parameters_from_hardcoded_table(cfg: dict) -> dict:
    base_year = cfg["BASE_COHORT_YEAR"]
    r_align = cfg["PRICE_YEAR_ALIGNMENT_RATE"]

    def align_money(value: float, source_year: Optional[int]) -> float:
        if source_year is None or r_align == 0:
            return float(value)
        delta = source_year - base_year
        if delta <= 0:
            return float(value)
        return float(value) / ((1.0 + r_align) ** delta)

    raw = RAW_PARAM_TABLE

    # Direct medical costs (2024)
    m24 = raw["cost_medical_2024"]
    costs_medical = {
        "inpatient_pneumonia": align_money(parse_numeric_any(m24["inpatient_pneumonia"]), 2024),
        "outpatient_pneumonia": align_money(parse_numeric_any(m24["outpatient_pneumonia"]), 2024),
        "inpatient_meningitis": align_money(parse_numeric_any(m24["inpatient_meningitis"]), 2024),
        "inpatient_sepsis": align_money(parse_numeric_any(m24["inpatient_sepsis"]), 2024),
        "outpatient_aom": align_money(parse_numeric_any(m24["outpatient_aom"]), 2024),
    }

    costs_sequelae = {
        "cognitive_annual": align_money(parse_numeric_any(m24["cognitive_annual"]), 2024),
        "hearing_case": align_money(parse_range_midpoint(m24["hearing_case_range"]), 2024),
        "epilepsy_case": align_money(parse_numeric_any(m24["epilepsy_case"]), 2024),
        "hemiplegia_case": align_money(parse_range_midpoint(m24["hemiplegia_case_range"]), 2024),
        "cochlear_implant_case": align_money(parse_range_midpoint(m24["cochlear_implant_case_range"]), 2024),
    }

    # Direct non-medical and caregiver costs (2025)
    n25 = raw["cost_nonmedical_2025"]
    c25 = raw["cost_caregiver_2025"]
    costs_nonmedical = {
        "inpatient_pneumonia": align_money(parse_numeric_any(n25["inpatient_pneumonia"]), 2025),
        "outpatient_pneumonia": align_money(parse_numeric_any(n25["outpatient_pneumonia"]), 2025),
        "inpatient_meningitis": align_money(parse_numeric_any(n25["inpatient_meningitis"]), 2025),
        "inpatient_sepsis": align_money(parse_numeric_any(n25["inpatient_sepsis"]), 2025),
        "outpatient_aom": align_money(parse_numeric_any(n25["outpatient_aom"]), 2025),
    }
    costs_caregiver = {
        "inpatient_pneumonia": align_money(parse_numeric_any(c25["inpatient_pneumonia"]), 2025),
        "outpatient_pneumonia": align_money(parse_numeric_any(c25["outpatient_pneumonia"]), 2025),
        "inpatient_meningitis": align_money(parse_numeric_any(c25["inpatient_meningitis"]), 2025),
        "inpatient_sepsis": align_money(parse_numeric_any(c25["inpatient_sepsis"]), 2025),
        "outpatient_aom": align_money(parse_numeric_any(c25["outpatient_aom"]), 2025),
    }
    costs_sequelae["special_education_annual"] = align_money(parse_range_midpoint(raw["special_education_annual_range_2025"]), 2025)

    # Vaccine/program costs
    vaccine_prices = {
        "pcv10": align_money(parse_numeric_any(raw["vaccine_price_pcv10_2024"]), 2024),
        "pcv13": align_money(parse_numeric_any(raw["vaccine_price_pcv13_unknown_year"]), None),
        "pcv10sii": align_money(parse_numeric_any(raw["vaccine_price_pcv10_2024"]), 2024),  # no separate table row
    }
    service_delivery_per_dose = align_money(parse_numeric_any(raw["service_delivery_per_dose_2024"]), 2024)
    ae_cost_per_case = align_money(parse_range_midpoint(raw["severe_ae_cost_range_unknown_year"]), None)
    ae_incidence_per_100k_doses = parse_numeric_any(raw["severe_ae_incidence_per_100k_doses"])
    wastage_rate = parse_numeric_any(raw["wastage_rate"])

    prod_age_contrib = {
        int(k): align_money(float(v), 2024) for k, v in raw["prod_discounted_age_contrib_15_60"].items()
    }
    prod_loss_per_death_age0 = sum(prod_age_contrib.values())

    params = {
        "birth_cohort": parse_numeric_any(raw["birth_cohort_2023"]),
        "coverage_rate": {
            "pcv10": parse_numeric_any(raw["coverage_pcv10_epi"]),
            "pcv10sii": parse_numeric_any(raw["coverage_pcv10_epi"]),
            "pcv13": parse_numeric_any(raw["coverage_pcv13_epi"]),
        },
        "doses_per_child": cfg["DOSES_PER_CHILD"],
        "costs_medical": costs_medical,
        "costs_nonmedical": costs_nonmedical,
        "costs_caregiver": costs_caregiver,
        "costs_sequelae": costs_sequelae,
        "vaccine_prices": vaccine_prices,
        "service_delivery_per_dose": service_delivery_per_dose,
        "wastage_rate": wastage_rate,
        "ae_incidence_per_100k_doses": ae_incidence_per_100k_doses,
        "ae_cost_per_case": ae_cost_per_case,
        "utilities": raw["utilities"].copy(),
        "life_expectancy": cfg["LIFE_EXPECTANCY_YEARS"],
        "discount_rate": cfg["DISCOUNT_RATE"],
        "discount_rate_align": cfg["PRICE_YEAR_ALIGNMENT_RATE"],
        "prod_discounted_age_contrib": prod_age_contrib,
        "prod_loss_per_death_age0": prod_loss_per_death_age0,
    }
    return params


# =========================
# Economics helpers
# =========================
def unit_cost_per_episode(params: dict) -> Dict[str, float]:
    uc = {}
    for k in ["inpatient_pneumonia", "outpatient_pneumonia", "inpatient_meningitis", "inpatient_sepsis", "outpatient_aom"]:
        uc[k] = params["costs_medical"][k] + params["costs_nonmedical"][k] + params["costs_caregiver"][k]
    return uc


def discount_factor(age_years_from_birth: float, r: float) -> float:
    if r == 0:
        return 1.0
    return 1.0 / ((1.0 + r) ** max(age_years_from_birth, 0.0))


def pv_future_qaly_from_death_age(age_at_death: float, life_expectancy: float, r: float) -> float:
    if age_at_death >= life_expectancy:
        return 0.0
    total = 0.0
    start_year = int(np.floor(age_at_death))
    end_year = int(np.ceil(life_expectancy))
    for y in range(start_year, end_year):
        seg_start = max(age_at_death, float(y))
        seg_end = min(life_expectancy, float(y + 1))
        if seg_end <= seg_start:
            continue
        duration = seg_end - seg_start
        mid_age = (seg_start + seg_end) / 2.0
        total += duration * discount_factor(mid_age, r)
    return total


def pv_prod_loss_from_death_age(age_at_death: float, prod_age_contrib: Dict[int, float]) -> float:
    cutoff = int(np.floor(age_at_death))
    total = 0.0
    for age, contrib in prod_age_contrib.items():
        if age > cutoff:
            total += contrib
    return total


def positive_increment(curr: float, prev: float) -> float:
    curr = 0.0 if (curr is None or np.isnan(curr)) else float(curr)
    prev = 0.0 if (prev is None or np.isnan(prev)) else float(prev)
    return max(curr - prev, 0.0)


def scenario_vaccine_key(file_name: str) -> Optional[str]:
    fname = os.path.basename(file_name).lower()
    if "baseline" in fname:
        return CONFIG["BASELINE_VACCINE_PROGRAM"]
    if "pcv13" in fname:
        return "pcv13"
    if "pcv10sii" in fname or "sii" in fname:
        return "pcv10sii"
    if "pcv10" in fname:
        return "pcv10"
    return None


def add_vaccine_program_pv(yearly_rows: List[dict], params: dict, vaccine_key: Optional[str], cfg: dict) -> None:
    if vaccine_key is None:
        return
    if vaccine_key not in params["vaccine_prices"]:
        raise ValueError(f"Unknown vaccine key: {vaccine_key}")

    cohort = params["birth_cohort"]
    cov = params["coverage_rate"][vaccine_key]
    doses = params["doses_per_child"]
    wastage = params["wastage_rate"]

    total_doses = cohort * cov * doses * (1.0 + wastage)
    vaccine_component = total_doses * (params["vaccine_prices"][vaccine_key] + params["service_delivery_per_dose"])
    ae_cases = total_doses * (params["ae_incidence_per_100k_doses"] / 100000.0)
    ae_component = ae_cases * params["ae_cost_per_case"]
    raw_total = vaccine_component + ae_component

    event_age = cfg["VACCINE_EVENT_AGE_YEARS"]
    pv_total = raw_total * discount_factor(event_age, params["discount_rate"])

    if not yearly_rows:
        return
    yearly_rows[0]["raw_vaccine_cost"] += raw_total
    yearly_rows[0]["pv_vaccine_cost"] += pv_total
    yearly_rows[0]["pv_total_cost"] += pv_total


def compute_yearly_pv(df: pd.DataFrame, params: dict, cfg: dict) -> pd.DataFrame:
    uc = unit_cost_per_episode(params)
    util = params["utilities"]
    los = cfg["LOS_DAYS"]
    r = params["discount_rate"]
    life_exp = params["life_expectancy"]
    prod_age_contrib = params["prod_discounted_age_contrib"]

    rows: List[dict] = []
    prev_stocks = {
        "Stock_Seq_Hearing": 0.0,
        "Stock_Seq_Epilepsy": 0.0,
        "Stock_Seq_Hemiplegia": 0.0,
        "Stock_Seq_Hearing_Implant": 0.0,
    }

    for i in range(len(df)):
        year_idx = i + 1
        row = df.iloc[i]
        cycle_mid_age = year_idx - 0.5
        cycle_df = discount_factor(cycle_mid_age, r)
        age_at_death = cfg["AVG_DEATH_AGE_BY_YEAR"].get(year_idx, cycle_mid_age)

        new_mening = 0.0 if np.isnan(row.get("New_Meningitis", np.nan)) else float(row.get("New_Meningitis", 0.0))
        new_psev = 0.0 if np.isnan(row.get("New_Pneumo_Sev", np.nan)) else float(row.get("New_Pneumo_Sev", 0.0))
        new_pnon = 0.0 if np.isnan(row.get("New_Pneumo_NonSev", np.nan)) else float(row.get("New_Pneumo_NonSev", 0.0))
        new_seps = 0.0 if np.isnan(row.get("New_Sepsis", np.nan)) else float(row.get("New_Sepsis", 0.0))
        new_aom = 0.0 if np.isnan(row.get("New_AOM", np.nan)) else float(row.get("New_AOM", 0.0))
        deaths = 0.0 if np.isnan(row.get("New_Deaths_Pneumo", np.nan)) else float(row.get("New_Deaths_Pneumo", 0.0))

        stock_cog = 0.0 if np.isnan(row.get("Stock_Seq_Cognitive", np.nan)) else float(row.get("Stock_Seq_Cognitive", 0.0))
        stock_hear = 0.0 if np.isnan(row.get("Stock_Seq_Hearing", np.nan)) else float(row.get("Stock_Seq_Hearing", 0.0))
        stock_epi = 0.0 if np.isnan(row.get("Stock_Seq_Epilepsy", np.nan)) else float(row.get("Stock_Seq_Epilepsy", 0.0))
        stock_hemi = 0.0 if np.isnan(row.get("Stock_Seq_Hemiplegia", np.nan)) else float(row.get("Stock_Seq_Hemiplegia", 0.0))
        stock_impl = 0.0 if np.isnan(row.get("Stock_Seq_Hearing_Implant", np.nan)) else float(row.get("Stock_Seq_Hearing_Implant", 0.0))
        stock_sped = 0.0 if np.isnan(row.get("Stock_Seq_Hearing_SpEd", np.nan)) else float(row.get("Stock_Seq_Hearing_SpEd", 0.0))

        new_hear = positive_increment(stock_hear, prev_stocks["Stock_Seq_Hearing"])
        new_epi = positive_increment(stock_epi, prev_stocks["Stock_Seq_Epilepsy"])
        new_hemi = positive_increment(stock_hemi, prev_stocks["Stock_Seq_Hemiplegia"])
        new_impl = positive_increment(stock_impl, prev_stocks["Stock_Seq_Hearing_Implant"])

        prev_stocks["Stock_Seq_Hearing"] = stock_hear
        prev_stocks["Stock_Seq_Epilepsy"] = stock_epi
        prev_stocks["Stock_Seq_Hemiplegia"] = stock_hemi
        prev_stocks["Stock_Seq_Hearing_Implant"] = stock_impl

        raw_disease_cost = (
            new_psev * uc["inpatient_pneumonia"]
            + new_pnon * uc["outpatient_pneumonia"]
            + new_mening * uc["inpatient_meningitis"]
            + new_seps * uc["inpatient_sepsis"]
            + new_aom * uc["outpatient_aom"]
        )

        s = params["costs_sequelae"]
        raw_sequela_cost_annual = (
            stock_cog * s["cognitive_annual"]
            + stock_sped * s["special_education_annual"]
        )
        raw_sequela_cost_onetime = (
            new_hear * s["hearing_case"]
            + new_epi * s["epilepsy_case"]
            + new_hemi * s["hemiplegia_case"]
            + new_impl * s["cochlear_implant_case"]
        )
        raw_sequela_cost = raw_sequela_cost_annual + raw_sequela_cost_onetime

        raw_death_prod_cost = deaths * pv_prod_loss_from_death_age(age_at_death, prod_age_contrib)
        pv_death_prod_cost = raw_death_prod_cost

        def acute_qaly(cases: float, los_days: float, utility: float) -> float:
            if cases <= 0:
                return 0.0
            return cases * (los_days / 365.0) * (1.0 - utility)

        raw_qaly_acute = (
            acute_qaly(new_psev, los["inpatient_pneumonia"], util["inpatient_pneumonia"])
            + acute_qaly(new_pnon, los["outpatient_pneumonia"], util["outpatient_pneumonia"])
            + acute_qaly(new_mening, los["meningitis"], util["meningitis"])
            + acute_qaly(new_seps, los["sepsis"], util["sepsis"])
            + acute_qaly(new_aom, los["aom"], util["aom"])
        )

        raw_qaly_sequela = (
            stock_cog * (1.0 - util["cognitive"])
            + stock_hear * (1.0 - util["hearing"])
            + stock_epi * (1.0 - util["epilepsy"])
            + stock_hemi * (1.0 - util["hemiplegia"])
        )

        qaly_loss_per_death_pv = pv_future_qaly_from_death_age(age_at_death, life_exp, r)
        raw_qaly_death = deaths * qaly_loss_per_death_pv
        pv_qaly_death = raw_qaly_death

        pv_disease_cost = raw_disease_cost * cycle_df
        pv_sequela_cost = raw_sequela_cost * cycle_df
        pv_qaly_acute = raw_qaly_acute * cycle_df
        pv_qaly_sequela = raw_qaly_sequela * cycle_df

        pv_total_cost = pv_disease_cost + pv_sequela_cost + pv_death_prod_cost
        pv_total_qaly_loss = pv_qaly_acute + pv_qaly_sequela + pv_qaly_death

        rows.append({
            "year_index": year_idx,
            "cycle_mid_age": cycle_mid_age,
            "avg_age_at_death_used": age_at_death,

            "New_Meningitis": new_mening,
            "New_Pneumo_Sev": new_psev,
            "New_Pneumo_NonSev": new_pnon,
            "New_Sepsis": new_seps,
            "New_AOM": new_aom,
            "New_Deaths_Pneumo": deaths,

            "raw_disease_cost": raw_disease_cost,
            "raw_sequela_cost": raw_sequela_cost,
            "raw_death_prod_cost_pvbase": raw_death_prod_cost,
            "raw_vaccine_cost": 0.0,

            "pv_disease_cost": pv_disease_cost,
            "pv_sequela_cost": pv_sequela_cost,
            "pv_death_prod_cost": pv_death_prod_cost,
            "pv_vaccine_cost": 0.0,
            "pv_total_cost": pv_total_cost,

            "raw_qaly_acute": raw_qaly_acute,
            "raw_qaly_sequela": raw_qaly_sequela,
            "raw_qaly_death_pvbase": raw_qaly_death,

            "pv_qaly_acute": pv_qaly_acute,
            "pv_qaly_sequela": pv_qaly_sequela,
            "pv_qaly_death": pv_qaly_death,
            "pv_total_qaly_loss": pv_total_qaly_loss,

            "cases_total": new_mening + new_psev + new_pnon + new_seps + new_aom,
            "deaths_total": deaths,

            "new_seq_hearing_proxy": new_hear,
            "new_seq_epilepsy_proxy": new_epi,
            "new_seq_hemiplegia_proxy": new_hemi,
            "new_seq_implant_proxy": new_impl,
        })

    return pd.DataFrame(rows)


def summarize_results(resdf: pd.DataFrame) -> dict:
    return {
        "PV_total_cost": float(resdf["pv_total_cost"].sum()),
        "PV_total_qaly_loss": float(resdf["pv_total_qaly_loss"].sum()),
        "Total_cases": float(resdf["cases_total"].sum()),
        "Total_deaths": float(resdf["deaths_total"].sum()),
    }


def process_all(files: List[str], params: dict, cfg: dict, model_name: str) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    model_idx = cfg["MODEL_INDEX"][model_name]
    summaries: List[dict] = []
    details: Dict[str, pd.DataFrame] = {}

    f0 = files[0]
    df0 = load_and_extract(f0, model_idx, EXPECTED_COLS)
    res0 = compute_yearly_pv(df0, params, cfg)
    rows0 = res0.to_dict("records")
    add_vaccine_program_pv(rows0, params, cfg["BASELINE_VACCINE_PROGRAM"], cfg)
    res0 = pd.DataFrame(rows0)
    s0 = summarize_results(res0)
    scen0 = os.path.basename(resolve_file(f0)).replace(".csv", "")
    details[scen0] = res0
    summaries.append({
        "file": os.path.basename(resolve_file(f0)),
        "scenario": scen0,
        **s0,
        "Delta_cost": np.nan,
        "Cases_averted": np.nan,
        "Deaths_averted": np.nan,
        "QALY_gained": np.nan,
        "USD_per_case": np.nan,
        "USD_per_death": np.nan,
        "USD_per_qaly": np.nan,
        "vaccine_used": cfg["BASELINE_VACCINE_PROGRAM"],
    })

    for f in files[1:]:
        df = load_and_extract(f, model_idx, EXPECTED_COLS)
        res = compute_yearly_pv(df, params, cfg)
        rows = res.to_dict("records")
        vkey = scenario_vaccine_key(f)
        add_vaccine_program_pv(rows, params, vkey, cfg)
        res = pd.DataFrame(rows)

        summ = summarize_results(res)
        delta_cost = summ["PV_total_cost"] - s0["PV_total_cost"]
        cases_averted = s0["Total_cases"] - summ["Total_cases"]
        deaths_averted = s0["Total_deaths"] - summ["Total_deaths"]
        qaly_gained = s0["PV_total_qaly_loss"] - summ["PV_total_qaly_loss"]

        summaries.append({
            "file": os.path.basename(resolve_file(f)),
            "scenario": os.path.basename(resolve_file(f)).replace(".csv", ""),
            **summ,
            "Delta_cost": delta_cost,
            "Cases_averted": cases_averted,
            "Deaths_averted": deaths_averted,
            "QALY_gained": qaly_gained,
            "USD_per_case": (delta_cost / cases_averted) if cases_averted != 0 else np.nan,
            "USD_per_death": (delta_cost / deaths_averted) if deaths_averted != 0 else np.nan,
            "USD_per_qaly": (delta_cost / qaly_gained) if qaly_gained != 0 else np.nan,
            "vaccine_used": vkey,
        })
        details[os.path.basename(resolve_file(f)).replace(".csv", "")] = res

    return pd.DataFrame(summaries), details


def save_outputs(summary_df: pd.DataFrame, details: Dict[str, pd.DataFrame], prefix: str) -> None:
    summary_df.to_csv(f"results_summary_{prefix}.csv", index=False)

    for scenario, df in details.items():
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", scenario).strip("_")
        df.to_csv(f"yearly_details_{prefix}_{safe}.csv", index=False)

    with pd.ExcelWriter(f"results_details_{prefix}.xlsx", engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary", index=False)
        for scenario, df in details.items():
            sname = re.sub(r"[^A-Za-z0-9_]+", "_", scenario)[:31] or "sheet"
            df.to_excel(writer, sheet_name=sname, index=False)


def print_parameter_audit(params: dict, cfg: dict) -> None:
    print("\n=== Parameter audit (hardcoded from your Excel table) ===")
    print(f"Birth cohort: {params['birth_cohort']}")
    print(f"Discount rate: {params['discount_rate']}")
    print(f"Life expectancy (QALY): {cfg['LIFE_EXPECTANCY_YEARS']}")
    print("Average death age by year:", cfg["AVG_DEATH_AGE_BY_YEAR"])
    print("Working ages for productivity loss: 15-60")
    print(f"Productivity loss per death at age 0 (PV@birth, aligned to 2023): {params['prod_loss_per_death_age0']:.2f}")
    print(f"Price-year alignment rate: {cfg['PRICE_YEAR_ALIGNMENT_RATE']}")
    print("Coverage:", params["coverage_rate"])
    print("Vaccine prices (aligned):", params["vaccine_prices"])
    print("LOS days:", cfg["LOS_DAYS"])


def main():
    resolved_files = []
    missing = []
    for f in CONFIG["FILES"]:
        try:
            resolved_files.append(resolve_file(f))
        except FileNotFoundError:
            missing.append(f)

    if missing:
        print("Missing input CSVs:")
        for m in missing:
            print("  -", m, "(or *_UPDATED.csv fallback)")
        print("Please place the files next to this script or edit CONFIG['FILES'].")
        return

    params = build_parameters_from_hardcoded_table(CONFIG)
    print_parameter_audit(params, CONFIG)

    for model_name in ["ODE", "Markov"]:
        print(f"\nProcessing {model_name} ...")
        summary_df, details = process_all(resolved_files, params, CONFIG, model_name)
        save_outputs(summary_df, details, model_name)
        print(summary_df[["scenario", "Delta_cost", "Cases_averted", "Deaths_averted", "QALY_gained", "USD_per_qaly"]])

    print("\nDone.")
    print("Generated:")
    print("  - results_summary_ODE.csv")
    print("  - results_summary_Markov.csv")
    print("  - results_details_ODE.xlsx")
    print("  - results_details_Markov.xlsx")
    print("  - yearly_details_*.csv (per scenario)")


if __name__ == "__main__":
    main()
