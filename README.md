# ProsodySemanticsIntegration — Scripts & Data

Behavioral and ROI analysis scripts and data for the prosody–semantics integration
(irony/sarcasm) studies.

## Structure

| Folder | Contents |
|---|---|
| `Study_1/` | `s1_analyses.R`, behavioral `data/`, `plots/` |
| `Study_2/` | `s2_analyses.R`, behavioral `data/`, `plots/` |
| `Study_3/` | `s3_behavioral_analyses.R`, `s3_ROI_analyses.R`, per-participant `data/` (p01–p50), `Betas/`, `plots/` |

## Notes

- Analyses are run in R. Each study's `.R` script reads from its `data/` folder and writes to `plots/`.
- R workspace (`.RData`) and history (`.Rhistory`) files are not tracked — results are regenerable from the scripts and data.
