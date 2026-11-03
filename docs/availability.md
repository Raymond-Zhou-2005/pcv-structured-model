# Data and code availability

## Version map

The aggregate inputs, 144 application configurations and retained numerical
results are publicly available in
[pcv-structured-model, v1.0.0](https://github.com/Raymond-Zhou-2005/pcv-structured-model/tree/v1.0.0).
The public main branch also includes later platform diagnostics.

The **v1.0.1** release incorporates the 2026-11-03-r1 manuscript and figure
revision, source data, saved-trajectory CFR translation audit, dependency records
and file hashes. Its fixed public archive is
[v1.0.1](https://github.com/Raymond-Zhou-2005/pcv-structured-model/releases/tag/v1.0.1).
The simulation inputs and outcomes remain unchanged from v1.0.0.
The current `release_manifest.json` identifies the exact files. The original
Overleaf PDF is retained as `manuscript_overleaf.pdf`; the publication PDF only
updates its availability paragraph to this release.

## Reuse and scope

Author-owned software is licensed under MIT; original data, figures and prose
are licensed under CC BY 4.0. Third-party rights and notices are retained in
`LICENSES.md`. The shared numerical materials contain aggregate scenario
constants and simulated outputs. Published sources and transformations are
identified in the parameter ledger and manuscript bibliography.

Original unavailable local administrative cost records are represented by
explicit scenario constants. The repository distinguishes source identity,
arithmetic reconstruction and independent parameter validation. The source
ledger retains the evidence status of every flattened input.

## Reproduction map

| Material | Location |
|---|---|
| Main manuscript and six publication figures | `paper/` |
| Figure/table data and script mapping | `docs/artifact_map.csv` |
| Units, denominator, array layout and seeds | `docs/data_dictionary.md` |
| Input, source and transformation ledger | `research/sources/parameter_sources.csv` |
| 144 saved configurations and outputs | `research/inputs/Final/outputs/` |
| Mathematical and historical controlled checks | `research/outputs/`, `research/coverage/outputs/`, `revision/outputs/` |
| CFR translation and revised-asset checks | `checks/revision_20261103/` |
| Numerical reference checks | `checks/verification.json` |
| File integrity | `release_manifest.json` |

After installing `requirements-lock.txt`, run `python verify_release.py` to
verify the downloaded archive and `python reproduce.py --stage check` to check
retained outputs and reconstruct the CFR diagnostic. `python reproduce.py`
reruns the complete default numerical and plotting pipeline. The optional
`--with-cohort` regenerates the historical closed-cohort repetitions.

The recorded reference environment is Windows x86-64 / CPython 3.12.14.
Historical Windows and Linux results and their limits are retained in
`checks/portability/README.md`. A new manuscript or figure revision does not
change the status of those earlier platform tests.

The software citation record is `CITATION.cff`. A persistent archive DOI can be
added after an actual deposit; none has been assigned to this release.
