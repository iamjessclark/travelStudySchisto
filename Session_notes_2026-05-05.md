# Session Notes — 2026-05-05
**Focus:** Reviewer responses — sensitivity analyses, discussion text, methods text

---

## 1. Sensitivity 2a: Measurement error / underreporting bounding (DurMin)

### Key conceptual change from previous session
The analysis was originally framed as "classical measurement error correction" using a reliability ratio λ. We revised this framing because:
- Classical ME assumes random, mean-zero noise — not applicable here
- The actual mechanism is **proportional systematic underreporting** due to social desirability bias (people want to please researchers by reporting less time in the water)
- The formula β_true = β_observed / λ is mathematically identical, but λ now represents a **reporting fraction** not a reliability ratio
- This distinction matters because proportional underreporting is a more honest and accurate description of the mechanism

### Critical correction to output table
The previous output showed **per-minute ORs (1.005–1.010)** which appear negligible but are misleading. At realistic visit durations the effect is substantial:

| Reporting fraction (λ) | Underreporting | OR at 30 min | OR at 60 min | OR at 120 min |
|---|---|---|---|---|
| 1.0 | 0% | 1.16 | 1.35 | 1.82 |
| 0.9 | 10% | 1.18 | 1.39 | 1.94 |
| 0.8 | 20% | 1.21 | 1.45 | 2.11 |
| 0.7 | 30% | 1.24 | 1.53 | 2.35 |
| 0.6 | 40% | 1.28 | 1.65 | 2.71 |
| 0.5 | 50% | 1.35 | 1.82 | 3.31 |

90th percentile of DurMin = 60 minutes. The conclusion is NOT "negligible" — the true effect of duration may be substantially larger than observed.

### Key structural point
Duration is a **mediator** on the Travel → Activity → Duration → Infection pathway. Measurement error in Duration does NOT affect the primary travel estimate — Duration is not in the total effect model.

### Non-differential assumption
The correction assumes underreporting is non-differential with respect to infection status. This is plausible because participants did not know their infection status at interview.

### Changes made to Sensitivity_analyses.R
- Header comment rewritten: "proportional underreporting" replaces "classical measurement error"
- Tibble now computes cumulative ORs at 30, 60, 120 minutes (not per-minute ORs)
- Includes 95% CrI at 60 min
- Variable names: `reporting_fraction`, `pct_underreporting`, `or_at_30min`, `or_at_60min_med/lo/hi`, `or_at_120min`
- Output CSV updated accordingly

### Changes made to session notes (Sensitivity_analyses_session_notes.md)
- Technical methods: reframed as proportional underreporting, non-differential assumption stated, mediator point explicit
- Layman's description: corrected "per-minute effect stays tiny" (this was wrong) — replaced with correct cumulative OR interpretation
- Results table: replaced per-minute OR table with cumulative OR table at 30/60/120 min

---

## 2. MDA social desirability bias

### Direction of bias
- Duration: underreporting → attenuation toward null (observed effect too small)
- MDA: **overreporting** → dilution of protective effect toward null (some untreated people say they were treated, dragging the "treated" group's infection rate up)
- The observed null MDA effect (OR 0.87, CrI 0.57–1.34) may therefore understate MDA's true protective effect
- This is a limitation, not a reassurance
- However: MDA's causal placement does not affect the travel estimate (Sensitivity 4 Model B: OR 1.74, CrI above null)

---

## 3. Reviewer comment mapping

### Comment 1: Temporal ordering (cross-sectional design)
- No analysis addresses this
- Answer in **Discussion/Limitations** using biological plausibility:
  - Schistosomiasis largely asymptomatic → participants don't know their status → can't modify behaviour
  - Water exposure is a biological precondition for infection → travel necessarily upstream
  - Any bias would be toward null (heavily infected = more morbidity = travel less)

### Comment 2: Unmeasured confounders
- Addressed by **Sensitivity 1 (E-values)**: E = 2.20 for daily travel → infection
- Brief mention in Discussion

### Comment 3: Social desirability — duration and MDA effect on duration-infection pathway
- Duration: addressed by **Sensitivity 2a** (underreporting bounding)
- MDA: addressed by overreporting argument (direction toward null) + Sensitivity 4 Model B
- Key point: total travel effect unaffected (Duration is mediator, not in total effect model)
- Addressed in **Limitations** section

### Comment 4: MDA participants differ systematically (confounding)
- Addressed by **Sensitivity 4 Model B**: treating MDA as confounder vs mediator leaves travel estimate unchanged (OR 1.74, CrI above null)
- E-value for MDA association shows it is fragile (E = 1.41) but the MDA effect is already null so this doesn't change conclusions
- Brief mention in Discussion

### Comment 5: Counterfactual assumptions (independent manipulation, behavioural adaptation, economic necessity)
- Independent manipulation → addressed by **Sensitivity 3 spillover analysis** ✓
- Behavioural substitution → addressed by **Sensitivity 3 substitution bounding** ✓
- Economic/social necessity → needs **Discussion text**: counterfactuals are idealised contrasts to quantify causal contribution, not prescriptions; practical interventions would be targeted (alternative water sources, protective equipment)

### Comment 6: DAG not empirically verifiable, mediation not established, language too strong
- Addressed by **Sensitivity 4** (all 8 structural combinations)
- Key finding: point estimate stable (OR 1.73–1.74) across all structures; only models not conditioning on Activity or Duration give CrI above null — consistent with over-adjustment when mediators incorrectly treated as covariates
- Needs **language change** throughout: "consistent with" / "compatible with" not "indicates" or "confirms" mediation
- DAG was constructed from domain knowledge, confirmed by published research, and validated using implied independencies in the data

---

## 4. Discussion text — what was drafted today

### Revised robustness results sentence (Results section)
Old (wrong): "produced negligible changes to the per-minute effect estimate"
New: "the direction of bias was toward the null, meaning observed duration–infection estimates are conservative lower bounds of the true effect. Proportional underreporting correction across reporting fractions of 0–50% indicated that corrected cumulative odds ratios at 60 minutes ranged from 1.35 under perfect reporting to 1.82 under 50% underreporting."

### Counterfactuals paragraph addition (end of paragraph 2)
> These counterfactuals represent idealised contrasts to quantify the causal contribution of travel rather than prescribe realistic interventions; spillover and substitution bounding confirmed that population-level reductions persisted even when independence assumptions were relaxed (Supplementary Table X).

### Limitations section — temporal ordering + DAG (combined, deduplicated)
Full draft written and agreed. Key structure:
1. Cross-sectional design → temporal ordering cannot be established
2. Reverse causation biologically implausible (asymptomatic, water exposure upstream, bias toward null if anything)
3. Real threat is unmeasured confounding — SES as candidate (does double duty: confounding threat AND alternative DAG structure threat)
4. E-value evidence (RR ≥ 2.20 needed) + geographic/socioeconomic homogeneity of study villages limits this
5. DAG grounded in domain knowledge + validated against implied independencies
6. Sensitivity 4: direction and magnitude preserved across all 8 structures; precision changes consistent with over-adjustment
7. Findings "consistent with, rather than conclusive evidence of" proposed mediation

### Limitations section — self-report / MDA (drafted)
- Duration most susceptible (continuous vs categorical); underreporting → conservative lower bound; total effect unaffected
- MDA overreporting → dilutes protective effect toward null; null MDA result may understate true benefit; structural sensitivity confirms travel estimate unaffected

---

## 5. Methods text for Sensitivity 2a (drafted, needs one gap filled)

A full methods paragraph was drafted. **Gap to fill**: what did the activity plausibility check find? (The check used granular pre-grouped activity data to assess whether reported durations were consistent with the activity type — e.g., fishing could logically involve hours in water, water collection could not. The result of this check needs to be stated.)

---

## 6. Files changed this session
- `Sensitivity_analyses.R` — Sensitivity 2a section rewritten (lines ~339–383)
- `Sensitivity_analyses_session_notes.md` — Methods (technical + layman) and Results table for Sensitivity 2a updated
