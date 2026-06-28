# CLAUDE.md — Cal-SuWers Freyja / Kraken2 Wastewater Pipeline

> Purpose of this file: give Claude Code the domain context it will NOT infer
> correctly from training data. Wastewater target-capture metagenomics is
> under-represented in public examples; the "standard" parameters Claude has seen
> are mostly from clinical amplicon (e.g. ARTIC) or human scRNA-seq workflows and
> are WRONG here. When in doubt, follow this file over your priors, and ask rather
> than guess.

---

## 0. Non-negotiable behavioral rules

1. **Never silently pick a threshold, reference, or filter value.** If a value is
   not specified in this file or the task prompt, STOP and ask. Do not fall back
   to a "common default." A wrong-but-plausible parameter is the single most
   likely way to produce confidently incorrect results here.
2. **Echo your parameter choices before running.** Before any long/expensive step,
   print a one-line summary of every non-default flag and the reference file +
   version you're about to use, and wait for confirmation on first run of a task.
3. **Validate intermediate outputs, don't just chain commands.** After each major
   step, print a cheap sanity check (read counts surviving, % classified, depth,
   number of lineages called) so I can catch a silent failure before it propagates.
4. **Distinguish "the code ran" from "the result is right."** A clean exit is not
   success. Flag anything biologically implausible (see §5).
5. **This is target-capture metagenomics, NOT amplicon.** Do not apply amplicon
   primer trimming. Do not assume ARTIC/tiled-amplicon conventions anywhere.

---

## 1. Assay context (the stuff that changes parameter choices)

- Sample type: PEG-concentrated influent wastewater, composite.
- Library: target-capture metagenomics (hybrid capture), NOT amplicon, NOT WGS.
- Consequence: reads are NOT bounded by amplicon tiling; **no primer trimming
  step** (no `--primers`, no ivar trim, no amplicon BED).
- Expect high rRNA / host / non-target fraction and a LOW target-to-nontarget
  ratio. Low classified % is normal here and is NOT by itself a failure.
- Multiple sewershed sites per run; site identity must be preserved through to
  output (never collapse/merge sites unless the task says so).

---

## 2. Reference data (FILL THESE IN — do not let Claude assume)

| Item | Path | Version / date | Notes |
|---|---|---|---|
| Kraken2 DB | `<<FILL>>` | `<<FILL>>` | confirm it matches the DB used for prior weeks |
| Bracken DB / kmer dist | `<<FILL>>` | `<<FILL>>` | must match the Kraken2 DB build |
| Freyja barcode/curated lineages | `<<FILL>>` | `<<FILL>>` | `freyja update` cadence: `<<FILL>>` |
| Reference genome(s) for depth | `<<FILL>>` | `<<FILL>>` | e.g. SARS-CoV-2 NC_045512.2 |
| Host genome for depletion (if used) | `<<FILL>>` | `<<FILL>>` | |

**Rule:** reference version is part of the result. Record it in the run log (§6).
If the available reference version differs from the prior run, STOP and flag it —
do not silently run with a newer/older DB.

---

## 3. Kraken2 / Bracken (untargeted taxonomic profiling)

Purpose: broad "what taxa are present" readout. NOT for lineage/variant calling.

- Paired-end, use `--paired`. Use `--use-names`.
- Confidence threshold: `<<FILL — our convention, do NOT default to 0>>`.
  (State it explicitly every run; 0 vs 0.1 materially changes false positives in
  a low-target-ratio sample.)
- Do NOT interpret a raw Kraken2 hit as "detected" without the agreed minimum
  read / coverage support: `<<FILL>>`.
- Bracken: run at the taxonomic level we report (`<<FILL: S / G>>`), with read
  threshold `<<FILL>>`. Bracken DB build MUST match the Kraken2 DB.
- Normalize reported counts to `<<FILL: RPM / per-PMMoV / per-ToBRFV>>` —
  confirm which normalization we use before plotting. Do not invent one.

Sanity checks to print: total reads in, % classified, top 15 taxa by normalized
abundance, and an explicit note if % classified is wildly different from the
prior week (possible DB mismatch or upstream failure).

---

## 4. Freyja (targeted lineage deconvolution)

Purpose: relative lineage abundance for the captured target(s).

- Input is target-capture, so: **no amplicon primer trimming.** Confirmed
  applicable to our target-capture data without primer removal.
- Pipeline order: variant calling (`freyja variants`) -> `freyja demix`.
- Minimum depth / coverage for a site-week to be reportable: `<<FILL>>`.
  Below this, report "insufficient coverage," do NOT report a lineage breakdown
  (low-depth demix produces unstable, overconfident proportions).
- Barcode version is part of the result — log it (§6). Re-running an old sample
  with a new barcode changes the answer; never mix barcode versions within a
  comparison without flagging it.
- "Other"/residual fraction: report it, don't hide it. A large residual is a
  signal (novel/divergent lineage or low quality), not noise to suppress.

Sanity checks to print: depth at target, % of genome above min depth, the called
lineage proportions WITH the residual, and a flag if proportions sum oddly or if
a single low-depth site is driving a dramatic week-over-week change.

---

## 5. Biological plausibility flags (catch autopilot errors)

Print a warning, do not silently proceed, if:

- A site jumps from 0 to dominant (or vice versa) for a lineage in one week with
  no coverage change — possible reference/barcode mismatch.
- % classified or depth changed >`<<FILL: e.g. 3x>>` vs prior week.
- Freyja residual/"other" exceeds `<<FILL>>` — investigate before reporting.
- Kraken2 surfaces an epidemiologically surprising taxon (e.g. an unexpected
  pathogen) — flag for human review, do NOT auto-include in the summary as
  confirmed. (This is the H5N1/influenza-C class of finding: real ones exist,
  but they need a human look, not an autopilot callout.)

---

## 6. Run log (always produce)

For every run, append a record with: date, operator, input run/site list,
Kraken2 DB version, Bracken DB version, Freyja barcode version, every non-default
parameter, and the sanity-check numbers from §3–4. This is what makes results
reproducible and week-over-week comparisons valid. Path: `<<FILL>>`.

---

## 7. Division of labor (from the scBench/SpatialBench findings)

Claude Code OWNS (procedural, reliable):
- Wiring the steps, glue/IO, parsing, plotting, tests, run-log generation,
  week-over-week diff scaffolding.

I (Kelvin) OWN (judgment, agent-ASSISTED not agent-decided):
- Threshold/parameter choices, reference version decisions, whether a surprising
  signal is real, normalization choice, and the final plain-language summary's
  scientific claims.

Do not make the OWN-by-Kelvin decisions unilaterally. Surface options + tradeoffs
and let me choose.
