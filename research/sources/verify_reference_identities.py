"""Cross-check retained article identities against exact DOI records in Europe PMC.

Does not infer numeric parameter support from a bibliographic match and does not
rewrite the manuscript bibliography. Metadata corrections remain reviewable.
"""
from pathlib import Path
import concurrent.futures, csv, html, json, re, requests

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent / "audit"
OLD = AUDIT / "references.json"

def normalized_title(title):
    title=re.sub(r"<[^>]+>","",html.unescape(title))
    return re.sub(r"[^a-z0-9]","",title.lower())

def query(record):
    doi = record.get("requested_doi") or record.get("metadata", {}).get("DOI", "")
    out = dict(key=record["key"], doi=doi, checked_date="2026-10-12", exact_match=False)
    try:
        r = requests.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={"query":'DOI:"'+doi+'"', "format":"json", "resultType":"core"}, timeout=28)
        r.raise_for_status()
        hits = [x for x in r.json().get("resultList",{}).get("result",[]) if x.get("doi", "").lower()==doi.lower()]
        if not hits:
            out["status"]="Exact DOI not found in queried index"; return out
        x=hits[0]; journal=x.get("journalInfo",{})
        out.update(exact_match=True,status="Exact DOI identity match",pmid=x.get("id"),title=x.get("title"),
            author_list=x.get("authorList",{}).get("author",[]),
            journal=journal.get("journal",{}).get("title"),year=x.get("pubYear"),
            volume=journal.get("volume"),issue=journal.get("issue"),pages=x.get("pageInfo"),
            primary_index_url="https://europepmc.org/article/MED/"+x.get("id",""))
        old=record.get("metadata",{})
        out["old_title"]=old.get("title",[""])[0]
        out["old_author_count"]=len(old.get("author",[]))
        out["new_author_count"]=len(out["author_list"])
        out["author_count_agrees"]=out["old_author_count"]==out["new_author_count"]
        out["title_agrees_ignoring_format"]=normalized_title(out["old_title"])==normalized_title(out["title"])
    except Exception as exc:
        out.update(status="Retrieval unavailable",error=str(exc))
    return out

def main():
    old=json.loads(OLD.read_text(encoding="utf-8"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        records=list(pool.map(query,old))
    manual={r["reference_key"]:r for r in json.loads((HERE/"primary_crosschecks.json").read_text())["records"]}
    for record in records:
        if not record["exact_match"] and record["key"] in ("Bennett2000","Tian2006"):
            record["manual_primary_crosscheck"]=manual[record["key"]]
            record["identity_status"]="Verified against manually inspected publisher/author primary record"
        else:
            record["identity_status"]=record["status"]
    (HERE/"reference_identities.json").write_text(json.dumps(records,indent=2,ensure_ascii=False),encoding="utf-8")
    keys=["key","doi","checked_date","exact_match","status","identity_status","pmid","year","volume","issue","pages","old_author_count","new_author_count","author_count_agrees","title_agrees_ignoring_format","primary_index_url"]
    with (HERE/"reference_identities.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows([{k:r.get(k,"") for k in keys} for r in records])
    print(json.dumps(dict(records=len(records),exact_identity_matches=sum(r["exact_match"] for r in records),
        manual_primary_identity_crosschecks=sum("manual_primary_crosscheck" in r for r in records),
        metadata_differences=[r["key"] for r in records if r["exact_match"] and (not r.get("author_count_agrees") or not r.get("title_agrees_ignoring_format"))],
        unavailable=[r["key"] for r in records if not r["exact_match"]],parameter_values_verified_by_identity_check=False),indent=2))

if __name__ == "__main__":
    main()
