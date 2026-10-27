"""Reconstruct parameter provenance without treating a citation identity as a numeric validation.

This script reads frozen inputs only and writes inside its own sources directory.
The internal comparison document is a consistency checkpoint, never a bibliography source.
The supplementary text is maintained separately in source_audit.tex.
"""
from pathlib import Path
import ast, csv, hashlib, json, math, re

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
PROJECT = CODE / "inputs/Ablation/inputs/project"
RUN = CODE / "inputs/Final/outputs/runs/pcv13_D1W1G1"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def literal(path, name):
    tree = ast.parse(path.read_text(encoding="utf-8-sig"))
    return next(ast.literal_eval(n.value) for n in tree.body
        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in n.targets))

def save_csv(name, rows):
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def numeric(value):
    if isinstance(value,(float,int)):
        return float(value)
    nums=[float(x) for x in re.findall(r"\d+(?:\.\d+)?",str(value).replace(",",""))]
    if not nums:
        raise ValueError(value)
    return nums[0]

def midpoint(value):
    if isinstance(value,(float,int)):
        return float(value)
    nums=[float(x) for x in re.findall(r"\d+(?:\.\d+)?",str(value).replace(",",""))]
    return sum(nums[:2])/min(2,len(nums))

def flatten(obj, prefix=""):
    if isinstance(obj, dict):
        for key, val in obj.items():
            yield from flatten(val, prefix + "/" + str(key) if prefix else str(key))
    elif isinstance(obj, (list, tuple)):
        for i, val in enumerate(obj):
            yield from flatten(val, prefix + "/" + str(i))
    else:
        yield prefix, obj

def epi_source(key):
    root = key.split("/")[0]
    if root == "birth_cohort":
        return "UNWPP", "https://population.un.org/wpp/", "Frozen demographic scenario", "Target-year country table has not been independently recovered."
    if root in ("u5_mortality_rate", "neonatal_mortality_rate"):
        return "WorldBankMortality", "https://data.worldbank.org/indicator/SH.DYN.MORT", "Frozen demographic scenario", "Numeric series vintage not independently recovered; all-cause mortality must not be added twice."
    if root in ("inc_pneumonia", "inc_meningitis", "prop_pneumo_severe"):
        return "Wahl2018", "https://doi.org/10.1016/S2214-109X(18)30247-X", "Transferred scenario input", "Source identity verified; exact country/age/year/table value not independently recovered."
    if root in ("inc_sepsis", "inc_aom_all"):
        return "PhilippinesHTA", "https://hta.dost.gov.ph/wp-content/uploads/2024/02/PCV-Full-HTA-Report-FINAL-.pdf", "Foreign transferred scenario input", "Primary assessment inspected; exact rate and age-band correspondence remain unverified."
    if root == "prop_aom_spn":
        return "Do2025", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11826232/", "Unverified attribution proxy", "Frozen theta=0.303; selected hospital culture proportion=0.301 and combined detection=0.411 do not identify a national attribution fraction."
    if root.startswith("cfr_"):
        return "Dilok2019", "https://doi.org/10.1016/j.vaccine.2019.06.015", "Transferred Thai scenario input", "Identity verified; exact retained CFR has not been independently matched to primary table."
    if root in ("prob_cog", "prob_epi", "prob_hemi", "prob_hear"):
        return "Edmond2010", "https://doi.org/10.1016/S1473-3099(10)70048-7", "Transferred sequela probability", "Identity verified; exact subtype and conditioning denominator unverified."
    if root == "prob_cochlear_implant":
        return "", "", "Unverified treatment-use scenario", "No primary denominator supporting 0.025 retrieved; this is conditional on incident hearing sequela."
    if root == "prob_special_education":
        return "", "", "Unverified service-use scenario", "No primary extract supporting 0.2574 retrieved; this is conditional on incident hearing sequela."
    if root in ("pcv10_sii_relative_ve_multiplier", "pcv10_sii_relative_ve_assumption") or (root == "vaccine_direct_effectiveness" and "/PCV10-SII/" in key):
        return "Adigweme2023", "https://doi.org/10.1016/S1473-3099(22)00734-4", "Retained summary multiplier / calculated endpoint profile", "Published post-booster immunogenicity motivates the bridge; original serotype weights and GMC extraction underlying 0.8203164085517026 are unavailable."
    if root == "vaccine_irrs":
        return "Bennett2025;Adigweme2023", "https://doi.org/10.1016/S1473-3099(24)00665-0", "Target population-effect scenario", "Different product-specific histories are imposed, including a constant mature SII series; exact target extraction unverified."
    if root == "vaccine_direct_effectiveness":
        return "", "", "Illustrative endpoint effectiveness profile", "Frozen profile is not a matched head-to-head clinical efficacy estimate."
    return "", "", "Structural / computational assumption", "Declared scenario coefficient; no independent empirical calibration claimed."

def econ_source(key):
    root = key.split("/")[0]
    if root in ("costs_medical", "costs_nonmedical", "costs_caregiver"):
        return "", "", "Unverified aggregate scenario constant", "Original hospital records / caregiver questionnaires were not available for independent inspection; no surrogate publication is assigned."
    if root == "costs_sequelae":
        if "cognitive_annual" in key:
            return "Pham2013", "https://doi.org/10.2185/jrm.2876", "Transferred costing scenario", "General childhood disability household cost; exact retained value, period and cognitive-sequela mapping unverified."
        return "", "", "Unverified transferred costing scenario", "Retained range midpoint or service cost; original primary tariff / study extract and currency conversion have not been independently verified."
    if root == "vaccine_prices":
        return "", "", "Unverified procurement scenario", "GSK and SII share the same assumed price; PCV13 4.22 is retained without a confirmed base-year conversion."
    if root == "service_delivery_per_dose":
        return "Mai2022", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9053159/", "Transferred delivery scenario", "Primary study concerns TT/Td delivery; its fiscal facility/outreach/school estimates are 1.76/3.86/1.49 USD in 2017, not a validated 2.10 USD infant-PCV unit cost."
    if root == "utilities":
        if key.endswith("/aom"):
            refs, url = "Oh1996", "https://pubmed.ncbi.nlm.nih.gov/8851461/"
        elif key.endswith("/outpatient_pneumonia"):
            refs, url = "vanHoek2012", "https://doi.org/10.1016/j.vaccine.2012.10.017"
        elif any(key.endswith("/"+s) for s in ("cognitive", "hearing", "epilepsy")):
            refs, url = "Oostenbrink2002", "https://pubmed.ncbi.nlm.nih.gov/12384194/"
        else:
            refs, url = "Bennett2000", "https://doi.org/10.1001/archpedi.154.1.43"
        bennett_map={"inpatient_pneumonia":"local infection", "sepsis":"hospitalization", "meningitis":"meningitis with recovery", "hemiplegia":"severe brain damage"}
        endpoint=key.split("/")[-1]
        if refs == "Bennett2000" and endpoint in bennett_map:
            return refs, "https://pubmed.ncbi.nlm.nih.gov/10632249/", "Numeric source match / transferred state valuation", "Retained number appears in primary abstract; source state is "+bennett_map[endpoint]+". Its endpoint and episode-duration transfer remain assumptions."
        return refs, url, "Transferred valuation scenario", "Reference identity verified; exact utility, elicitation method and duration mapping remain unverified."
    if root.startswith("prod_"):
        return "", "", "Computationally reconstructed valuation scenario", "46 stored annual contributions reconstruct the stated employment bands and an implied unrounded income; national primary series and currency conversion remain unverified."
    if root in ("ae_cost_per_case", "ae_incidence_per_100k_doses"):
        return "", "", "Unverified adverse-event scenario", "No original local cost record or vaccine-specific surveillance rate denominator recovered."
    return "", "", "Structural / computational assumption", "Declared costing setting; no independent empirical validation claimed."

def main():
    HERE.mkdir(parents=True, exist_ok=True)
    p = json.loads((RUN / "parameters.json").read_text())
    e = json.loads((RUN / "economic_parameters.json").read_text())
    raw = literal(PROJECT / "QALY_new.py", "RAW_PARAM_TABLE")
    cfg = literal(PROJECT / "Refined/Week 7_qaly_refinement_impact.py", "CALENDAR_CONFIG")
    rows = []
    for group, values, classifier in (("epidemiology/programme", p, epi_source), ("economics", e, econ_source)):
        for key, value in flatten(values):
            reference, url, status, note = classifier(key)
            rows.append(dict(group=group, parameter=key, value=value,
                reference_key=reference, primary_url=url, numeric_evidence_status=status,
                source_identity_is_parameter_validation=False, qualification=note))
    for key, value in flatten(cfg):
        rows.append(dict(group="hidden economic configuration", parameter=key, value=value,
            reference_key="", primary_url="", numeric_evidence_status="Structural / computational assumption",
            source_identity_is_parameter_validation=False,
            qualification="LOS duration is frozen and no primary duration extract was recovered." if key.startswith("LOS_DAYS") else "Explicit modelling setting."))
    for disease, bands in {
        "meningitis": [(0,11,.65),(12,23,.18),(24,35,.07),(36,59,.10)],
        "pneumonia/sepsis": [(0,11,.50),(12,23,.25),(24,35,.15),(36,59,.10)],
        "AOM": [(0,11,.1869),(12,23,.1889),(24,59,.6242)]}.items():
        for left,right,mass in bands:
            rows.append(dict(group="hidden disease-age allocation", parameter=f"{disease}/{left}-{right} months mass", value=mass,
                reference_key="", primary_url="", numeric_evidence_status="Unverified age-allocation scenario",
                source_identity_is_parameter_validation=False, qualification=f"Uniform monthly weight {mass/(right-left+1):.17g} within band; all bands sum to one. No country/cause/measure query extract recovered."))
    save_csv("parameter_sources.csv", rows)
    # Reproduce the retained endpoint multiplication only. This does not recover
    # the original serotype weights or GMC extraction behind the summary proxy.
    sii_multiplier = float(p["pcv10_sii_relative_ve_multiplier"])
    sii_bridge = []
    for endpoint in ("ipd", "pneumonia", "aom"):
        pcv13_profile = float(p["vaccine_direct_effectiveness"]["PCV13"][endpoint])
        derived_sii_profile = pcv13_profile * sii_multiplier
        frozen_sii_profile = float(p["vaccine_direct_effectiveness"]["PCV10-SII"][endpoint])
        assert math.isclose(derived_sii_profile, frozen_sii_profile,
                            rel_tol=1e-14, abs_tol=1e-15), endpoint
        sii_bridge.append(dict(endpoint=endpoint, pcv13_profile=pcv13_profile,
            summary_multiplier=sii_multiplier, derived_sii_profile=derived_sii_profile,
            verification_scope="calculation_only",
            original_derivation_status="unrecovered"))
    save_csv("sii_bridge.csv", sii_bridge)
    # Numeric compatibility uses the frozen table values, not table citations as scientific evidence.
    comp=[]
    def compare(name, raw_value, comparator, runtime, source_year, note=""):
        comp.append(dict(parameter=name, raw_constant=raw_value, comparison_table_value=comparator,
            runtime_value=runtime, labelled_source_year=source_year,
            raw_difference=float(raw_value)-float(comparator), agreement=math.isclose(float(raw_value),float(comparator),rel_tol=1e-10,abs_tol=1e-9),
            comparison_document="internal manuscript Table 1, pp 21-24; consistency check only",
            publication_citation="Original source or scenario declaration only", note=note))
    for key,val in {"birth_cohort":1384610,"u5_mortality_rate":20,"neonatal_mortality_rate":10.3,
        "inc_pneumonia":1771,"inc_meningitis":12,"inc_sepsis":124,"inc_aom_all":9570,
        "prop_pneumo_severe":.389,"prop_aom_spn":.303,"cfr_pneumo_severe":.014241,"cfr_meningitis":.032258,"cfr_sepsis":.08,
        "prob_cog":.016,"prob_epi":.022,"prob_hemi":.032,"prob_hear":.046,"prob_cochlear_implant":.025,"prob_special_education":.2574,"vaccine_coverage":.9}.items():
        compare(key,p[key],val,p[key],"2023" if key in ("birth_cohort","u5_mortality_rate","neonatal_mortality_rate") else "mixed / unspecified")
    for category,values,year in (("costs_medical",[221.09,16.09,574.5,589.6,13.31],2024),
        ("costs_nonmedical",[142.33,18.98,297.95,305.24,21.21],2025),
        ("costs_caregiver",[184.25,124.49,297.82,193.55,139.12],2025)):
        for endpoint,val in zip(("inpatient_pneumonia","outpatient_pneumonia","inpatient_meningitis","inpatient_sepsis","outpatient_aom"),values):
            note="Internal table footnote says outpatient pneumonia equals AOM; displayed values differ. Frozen displayed values are retained as separate scenarios." if endpoint=="outpatient_pneumonia" and category in ("costs_nonmedical","costs_caregiver") else "Original records unavailable; 3% price-year conversion is an illustrative proxy."
            rawkey={"costs_medical":"cost_medical_2024", "costs_nonmedical":"cost_nonmedical_2025", "costs_caregiver":"cost_caregiver_2025"}[category]
            compare(category+"/"+endpoint,numeric(raw[rawkey][endpoint]),val,e[category][endpoint],year,note)
    for endpoint,val,year in (("cognitive_annual",73.36,2024),("hearing_case",1020,2024),("epilepsy_case",329.06,2024),("hemiplegia_case",998.475,2024),("cochlear_implant_case",24829.415,2024),("special_education_annual",2530.115,2025)):
        if endpoint=="special_education_annual":
            observed=midpoint(raw["special_education_annual_range_2025"])
        else:
            rawkey={"hearing_case":"hearing_case_range", "hemiplegia_case":"hemiplegia_case_range", "cochlear_implant_case":"cochlear_implant_case_range"}.get(endpoint,endpoint)
            observed=midpoint(raw["cost_medical_2024"][rawkey])
        compare("costs_sequelae/"+endpoint,observed,val,e["costs_sequelae"][endpoint],year,"For ranges the raw value is the arithmetic midpoint; source-year labels follow the frozen code, not a verified price index.")
    for endpoint,val in raw["utilities"].items():
        compare("utilities/"+endpoint,val,val,e["utilities"][endpoint],"not applicable","Dimensionless transferred utility; exact primary table extraction remains unverified.")
    compare("vaccine_prices/pcv10sii",numeric(raw["vaccine_price_pcv10_2024"]),2,e["vaccine_prices"]["pcv10sii"],2024,"Price identity/procurement terms not independently verified.")
    compare("vaccine_prices/pcv10",numeric(raw["vaccine_price_pcv10_2024"]),2,e["vaccine_prices"]["pcv10"],2024,"GSK price is the same retained proxy as SII; no separate original price row.")
    compare("vaccine_prices/pcv13",numeric(raw["vaccine_price_pcv13_unknown_year"]),4.22,e["vaccine_prices"]["pcv13"],"unknown in code; table labels NIHE2025","Retained unchanged because source year was coded unknown; the whole ledger uses scenario USD.")
    compare("service_delivery_per_dose",numeric(raw["service_delivery_per_dose_2024"]),2.1,e["service_delivery_per_dose"],2024,"Mai2022 supplies delivery context, not verified infant-PCV 2.10.")
    compare("ae_cost_per_case",midpoint(raw["severe_ae_cost_range_unknown_year"]),4800.355,e["ae_cost_per_case"],"unknown","Arithmetic midpoint of 4717.59 and 4883.12; local original record unavailable.")
    compare("ae_incidence_per_100k_doses",numeric(raw["severe_ae_incidence_per_100k_doses"]),4.01,e["ae_incidence_per_100k_doses"],"2022 surveillance cited in table","Exact vaccine-specific denominator has not been recovered.")
    compare("wastage_rate",numeric(raw["wastage_rate"]),.05,e["wastage_rate"],"structural","Operational mark-up convention uses procurement cost times 1+w; not 1/(1-w).")
    save_csv("input_comparison.csv",comp)
    # Hidden productivity ledger reconstructed independently from the frozen array.
    income=224.08583333333333; discount=.03; wage_rows=[]
    for age,stored in raw["prod_discounted_age_contrib_15_60"].items():
        unemployment=.0792 if age<25 else .016 if age<50 else .0141
        reconstructed=12*income*(1-unemployment)/(1+discount)**age
        wage_rows.append(dict(age=age,unemployment=unemployment,implied_monthly_income=income,
            stored_birth_discounted_2024=stored,reconstructed_birth_discounted_2024=reconstructed,
            difference=stored-reconstructed,proxy_aligned_2023=e["prod_discounted_age_contrib"][str(age)]))
    save_csv("productivity_age_weights.csv",wage_rows)
    maxerr=max(abs(x["difference"]) for x in wage_rows)
    raw_sum=sum(x["stored_birth_discounted_2024"] for x in wage_rows)
    event=float(e["prod_loss_per_death_age0"])*(1.03**2.5)
    qaly=sum(1.03**-(k+.5) for k in range(71))+.5*(1.03**-71.25)
    assert maxerr<2e-9
    assert all(math.isclose(x["proxy_aligned_2023"],x["stored_birth_discounted_2024"]/1.03,rel_tol=1e-13) for x in wage_rows)
    result=dict(date="2026-10-12",flattened_parameter_rows=len(rows),raw_table_comparison_rows=len(comp),
        raw_numeric_matches=sum(bool(x["agreement"]) for x in comp),productivity_reconstruction_max_error=maxerr,
        raw_productivity_sum=raw_sum,proxy_aligned_productivity_sum=e["prod_loss_per_death_age0"],
        event_relative_productivity_at_2p5=event,event_relative_qaly_at_2p5=qaly,
        sii_bridge_endpoint_rows=len(sii_bridge),
        sii_bridge_max_difference_to_frozen_profile=max(abs(x["derived_sii_profile"]-float(p["vaccine_direct_effectiveness"]["PCV10-SII"][x["endpoint"]])) for x in sii_bridge),
        sii_bridge_evidence_status="calculation_only; original serotype-weight/GMC derivation unrecovered",
        using_rounded_income224p09_productivity=sum(12*224.09*(1-x["unemployment"])/1.03**x["age"] for x in wage_rows),
        numeric_validation_scope="Code/table arithmetic only; source population and exact primary values remain qualified.",
        public_manuscript_internal_comparison_citations=False)
    (HERE/"source_verification.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    inputs=[RUN/"parameters.json",RUN/"economic_parameters.json",PROJECT/"QALY_new.py",PROJECT/"Refined/Week 7_qaly_refinement_impact.py",PROJECT/"Refined/Week 4_calendar_original.py"]
    (HERE/"provenance.json").write_text(json.dumps([dict(path=str(x.relative_to(CODE)),sha256=sha(x)) for x in inputs],indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))

if __name__ == "__main__":
    main()
