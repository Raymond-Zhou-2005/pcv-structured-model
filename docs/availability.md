# Data and code availability

## Manuscript statement

All aggregate model inputs, run-level simulated states, disease-event and
economic outputs, source and transformation ledgers, figure and table source
data, and author-generated analysis and verification code underlying the
reported results are publicly available at
https://github.com/Raymond-Zhou-2005/pcv-structured-model, release **v1.0.0**
(2026-10-28). The fixed snapshot is
https://github.com/Raymond-Zhou-2005/pcv-structured-model/tree/v1.0.0.
The repository includes the manuscript source, dependency versions, seed rules
and instructions for reproducing the analyses. It contains aggregate scenario
constants and simulated outputs; no individual patient records are included.
Original published sources and input transformations are identified in the
parameter ledger and manuscript bibliography. Original code is licensed under
MIT, and author-owned original data and prose under CC BY 4.0; third-party
rights and notices are retained as described in `LICENSES.md`. Unavailable
original local cost records are represented only by disclosed scenario
constants and are not claimed as newly collected or independently accessed data.

No DOI has been assigned to this software/data snapshot. The code citation
record is `CITATION.cff`; it does not represent an accepted journal article.

## What supports the results

| Required material | Location |
|---|---|
| Figure/table-to-source-data and generation-script mapping | `docs/artifact_map.csv` |
| Units, denominators, array definitions, seeds and estimands | `docs/data_dictionary.md` |
| Complete parameter/source/transformation ledger and qualifications | `research/sources/parameter_sources.csv` |
| Current 144 configurations and run-level records | `research/inputs/Final/outputs/` |
| Mathematical and historical controlled checks | `research/outputs/`, `research/coverage/outputs/`, `revision/outputs/` |
| Complete analysis and figure regeneration | `python reproduce.py` after installing `requirements.txt` with Python 3.12 |
| Retained-output validation | `python reproduce.py --stage check` |
| Main manuscript and publication figures | `paper/` and release download assets |
| Current file integrity | `release_manifest.json` |

The repository records computational reproducibility and the evidence limits of
the scenario inputs. It does not establish independent national calibration,
validate every transferred cost or utility constant, or provide the missing
joint calendar-model parameter PSA.

## Sharing status and remaining submission tasks

The depositor approved public upload and the specified licences, and confirmed
permission to share inherited aggregate inputs on 2026-10-28. That
representation is distinct from an independent legal audit. The original source
rights remain applicable. Third-party publication identities are documented;
full-text articles and original unavailable administrative records are not
represented as deposited datasets.

The [PLOS data policy](https://journals.plos.org/globalpublichealth/s/data-availability)
requires the data and metadata needed to replicate reported findings and a
specific availability statement. The
[PLOS code policy](https://journals.plos.org/globalpublichealth/s/materials-software-and-code-sharing)
expects reusable, documented author-generated code and encourages an open
licence and persistent identifier. The paths above provide the sharing map;
the GitHub release fixes a version but is not a DOI archive.

For stronger long-term preservation, deposit the exact release in a repository
such as Zenodo and add its actual DOI to the manuscript and citation metadata.
No archive record or DOI is fabricated here. This is a preservation improvement,
not evidence of PLOS acceptance. Complete authorship, affiliations, funding,
competing interests and acknowledgment permissions still require the authors'
confirmation in the manuscript and submission form.

The review manuscript has page numbers, continuous line numbers and double
spacing. The retained TIFF figures are single-page RGB, LZW-compressed,
600 dpi and below 10 MB each. Some original canvas widths and small plot labels
still need adjustment to the final production limits (maximum 7.5-inch width,
8–12-point Arial/Times/Symbol text). PLOS waives these formatting requirements
before provisional Editorial Accept; this repository does not certify those
figures as production-approved. See the [PLOS figure specifications](https://journals.plos.org/globalpublichealth/s/figures).
