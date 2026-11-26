# Waiting-time-code

Reproducible numerical pipeline for the manuscript:

> **"Can Mutation and Selection Assemble Human Innovations in the Time Available? A Quantitative Re-examination of the Waiting-Time Problem"**
> Atikur Rahman (2026)

This repository is the complete deposited code and output for the paper's quantitative core — corresponding to **Supplementary File S1** (the analysis script) and **Supplementary File S2** (the complete output log) of the manuscript, mirrored on Zenodo ([DOI: 10.5281/zenodo.24096351](https://doi.org/10.5281/zenodo.24096351)).

## What it does

The script `waiting.py` implements every numerical claim in the manuscript:

- Each proposed resolution mechanism of the evolutionary *waiting-time problem* — single substitutions, stochastic tunneling, duplication-based accretion, standing variation, de novo ORF supply — is priced against an empirically constrained hominin time budget.
- Arrival times (de novo emergence of a specified genotype) are separated from establishment times (fixation).
- Published numbers quoted from other papers (Durrett & Schmidt 2008; Sanford et al. 2015; Lynch 2005; Kimura & Ohta 1969, and others) are reproduced to the precision shown in those papers and **labelled QUOTED**. All other numbers are **COMPUTED** by the stated formulas.

Running the script writes the complete plain-text log `numerical_outputs.txt` — every waiting time, table cell, scenario estimate and appendix number in the manuscript is emitted by this one script.

## Repository contents

| File | Description |
| --- | --- |
| `waiting.py` | The single numerical pipeline (pure Python standard library, no dependencies) |
| `numerical_outputs.txt` | Complete output log of the script (= manuscript Supplementary File S2) |
| `LICENSE` | MIT License |

## Fixed conventions

| Symbol | Meaning | Central value |
| --- | --- | --- |
| `Ne` | effective population size, diploid individuals | 1.0 × 10⁴ |
| `2Ne` | gene copies producing mutations each generation | 2.0 × 10⁴ |
| `mu` | per-site mutation rate per generation | 1.25 × 10⁻⁸ |
| `B` | time budget | 2.4 × 10⁵ generations (6 Myr at 25 yr/generation) |
| `B_gen` | generous budget | 4.5 × 10⁵ generations (9 Myr at 20 yr/generation) |
| `omega_k` | functional fraction of k-mer space, C_k / 4^k | scenario-dependent |
| `T_arr(k)` | de novo first-arrival time = 1 / (2Ne · mu^k · omega_k) | computed per scenario |

## Reproduce

```bash
python3 waiting.py
```

Requirements: **Python 3.6+** (standard library only — no packages to install).
The script prints the log to stdout and writes `numerical_outputs.txt` in its own directory; the committed log was produced by this exact code.

## What the log contains

`numerical_outputs.txt` (191 lines) is organized into labelled sections:

- **0. Parameters** — Table 1 of the manuscript
- **A.1 Class I** — single substitutions: arrival + fixation, favored (s = 0.01) and neutral cases
- **A.3** — segregating probability of a specified k-tuple (Eq. 4)
- **Table 3** — first-arrival times T = 4^k / (2Ne · mu^k · C_k) across Ne = 10⁴–10¹⁰, k = 2–8, C_k = 1 and 10³
- **Per-position multiplier** — 4/mu (C_k fixed) and 1/mu (omega_k fixed) per added specified position
- **Scenario 1** — 8-nucleotide micro-RNA seed match, presite-rich vs presite-free regimes
- **Scenario 2** — Durrett–Schmidt two-mutation regulatory case (quoted) + census-sensitivity axis
- **A.2** — polynomial tunneling correction at the DS parameters (~10⁴ factor at k = 2)
- **Scenario 5** — composite best case: all mechanism families operating jointly, presite-free 8-mer
- **Scenario 3 / A.5** — Behe–Snoke two-residue feature via duplication, extrapolated to hominin Ne
- **Scenario 4 / A.4** — de novo beneficial ORF supply (Keefe–Szostak and Axe benchmark densities)
- **Table 2** — Class IV conversion under the generous trial ceiling (θ = 1 novel chain/genome/generation)
- **Sanford et al. 2015** — Mendel's Accountant hominin run (quoted)
- **Figure 2 / Figure 3 data** — the plotted values for the manuscript figures

## Headline numbers (central parameters)

- Presite-bearing regulatory-motif turnover **fits the budget**: an 8-mer completed from a 7-of-8 presite takes 3.75 × 10⁵ years ≈ 1.5 × 10⁴ generations ≈ 0.06 × B (Durrett & Schmidt 2008, presite-rich regime).
- The one measured rate of de novo gene birth is **marginally compatible** with mutation-limited supply at the generous census reading.
- A **conditional residual** remains for presite-free specified strings: even when every mechanism family operates jointly at best case, ≥ 10¹⁸ × the budget at k = 5, rising to ≥ 10⁴¹ × at k = 8.
- Novel structured folds retain ≥ 10⁶⁷ × the budget if they require Axe-class specificity.

All figures are conditional on the stated parameters; the manuscript examines the extrapolation scope of each anchor. See the manuscript for full framing and derivations.

## Related

- **Manuscript**: "Can Mutation and Selection Assemble Human Innovations in the Time Available? A Quantitative Re-examination of the Waiting-Time Problem" (under submission, BioSystems — transferred from the Journal of Theoretical Biology).
- **Zenodo deposit** (script + output, Version 1.0.0): [DOI: 10.5281/zenodo.24096351](https://doi.org/10.5281/zenodo.24096351)

## Author

**Atikur Rahman** — Independent Researcher, Dhaka, Bangladesh

## License

[MIT](LICENSE) © 2026 Atikur Rahman
