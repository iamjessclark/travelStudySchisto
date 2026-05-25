# Session Notes — 2026-05-25
**Focus:** Expert feedback from Simon Babayan on sensitivity analyses; language revisions to manuscript, supplementary, and reviewer response letter

---

## 1. Simon's feedback (received via Slack)

Simon reviewed the revised sensitivity analyses and provided the following:

**General**: Both analyses (spillover bounding and measurement error bounding) should be described consistently as bounding or sensitivity analyses — not as if a transmission network was fitted or true reporting fractions were measured from the data. The fact that results are "unsurprising" is not a problem; that is the point of these analyses.

**Table S12 (E-values)**: Strong, no changes needed.

**Table S14 (DAG structural)**: Strong; directly addresses Reviewer 3's concern; supports cautious language ("consistent with mediation").

**Table S11 / Figure S8 (spillover bounding)**:
- Rename to "environmental benefit transfer" or "shared-water bounding scenarios" — avoid implying epidemiological spillover was estimated from data
- State explicitly that 10–30% fractions are illustrative bounds
- In response letter and briefly in limitations: note spillover scenarios are optimistic and substitution scenarios are pessimistic — together they bracket the baseline estimate
- Check table reference in reviewer response (Simon flagged possible S10/S11 confusion)

**Table S13 (measurement error)**:
- Keep in discussion and as a table — appropriate for Reviewer 3
- Trim specific percentages and OR values from main text limitations — keep full detail in supplement

---

## 2. Changes made to manuscript (Word documents — pasted in by RML)

### Reviewer response letter — Comment 23
- Fixed wrong table reference: `Table S10` → `Table S11` (S10 is the duration/activity table; S11 is the spillover bounding table)
- Replaced "The spillover analysis..." with "Shared-water bounding scenarios (Table S11 and Figure S8)..."
- Added that fractions are illustrative and labelled scenarios as optimistic/pessimistic
- Last sentence reworded by RML to be appropriately hedged: "together confirm that the conclusion would likely hold or have greater benefit across plausible violations of the independence assumption"

### Main manuscript — Methods paragraph [138]
- Replaced "potential spillover effects" with "shared-water environmental benefit"
- Labelled fractions as "illustrative"
- Added behavioural substitution bounding to the same paragraph

### Main manuscript — Results paragraph [186]
- Replaced "Modelling spillover fractions of 10–30%" with "Shared-water bounding scenarios applying illustrative benefit-transfer fractions of 10–30%..."
- Added brief description of what the fractions represent (cercarial contamination reduction assumed to benefit non-travellers)

### Main manuscript — Robustness results paragraph [210]
- Removed: specific reporting fractions (0–50%), ORs at 60 minutes (1.34–1.82), and 90th percentile reference
- Retained: directional claim (bias toward null → conservative lower bounds) and pointer to Table S13
- New text: "A bounding analysis confirmed this, with full results in Table S13."

### Main manuscript — Limitations paragraph [230]
- Added explicit optimistic/pessimistic framing
- Added statement that fractions are illustrative bounds, not estimated from data
- Replaced "sensitivity analyses bounding the impact of these differences" with "bounding analyses that bracket the baseline estimate from both directions"

### Table S11 caption (supplementary)
- Renamed to "Shared-water bounding analysis for environmental benefit transfer and behavioural substitution effects..."
- Added sentence stating fractions are illustrative bounds, not data-derived
- Added optimistic/pessimistic direction for each scenario type

---

## 3. Changes made to Sensitivity_analyses.R

Four reader-visible label changes in the Sensitivity 3 section (lines ~449–500):

| Location | Old | New |
|---|---|---|
| `scenario` column (CSV / Table S11) | `"% spillover benefit to non-daily-travellers"` | `"% shared-water benefit transfer to non-travellers"` |
| `type` column (plot legend / Fig S8) | `"Spillover to non-travellers"` | `"Shared-water benefit transfer"` |
| Plot x-axis label | `"Assumed fraction (spillover or substitution)"` | `"Illustrative assumed fraction"` |
| Plot title | `"Sensitivity 3: spillover and substitution bounding"` | `"Sensitivity 3: shared-water bounding and behavioural substitution"` |

Internal variable names (`spillover_fracs`, `spillover_results`, etc.) were left unchanged as they are not reader-visible.

**Action needed**: re-run `Sensitivity_analyses.R` to regenerate `Figures/spillover_bounding.png` and `Figures/spillover_bounding.csv` with updated labels before final submission.

---

## 4. Files changed this session
- `Sensitivity_analyses.R` — Sensitivity 3 plot/CSV labels updated
- `Session_notes_2026-05-25.md` — this file
- Manuscript Word documents edited by RML (not tracked in this repo)
