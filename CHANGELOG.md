# Changelog

## v1.1.0 — 2026-05-15

### Changed
- **`MAFFT_ScoreNGo.py` is now an alignment generator only.** The simple BLOSUM62-based scoring previously embedded in this script has been removed.
- The script now supports **batch processing** of multiple FASTA files in a single run, creating separate output directories per input.
- At the end of a run, users are pointed to `AlignmentScorerWithJalview.py` for alignment evaluation.

### Added
- **`AlignmentScorerWithJalview.py`** — new standalone scoring tool implementing:
  - Exact Jalview conservation algorithm (per-column 0–11 scores, 10 physicochemical property classes, AMAS-derived)
  - Exact Jalview quality algorithm (BLOSUM62-based)
  - Sum-of-Pairs score (BLOSUM62)
  - Hydrophobic position conservation
  - Gap coherence
  - Parsimony-informative sites
  - Pairwise identity
  - Three ranking methods (Weighted Average, TOPSIS, Borda Count)
  - Consensus best-alignment detection
  - Parallel processing with chunked execution for large alignment sets
- Added `numpy`, `scipy` and `tqdm` to dependencies (for the new scorer).

### Removed
- `evaluate_alignment()` and related ranking output from `MAFFT_ScoreNGo.py`.
- `mafft_results_summary.txt` is no longer produced by `MAFFT_ScoreNGo.py`; equivalent (and more detailed) ranking output is now produced by the scorer.
- `MAFFT_ScoreNGo_Documentation.pdf` removed (the README and PARAMETERS.md now serve as the canonical documentation).

### Notes
- Users who relied on the previous in-script scoring should migrate to the new two-step workflow described in the README.
- The French legacy version `MAFFT_ScoreNGo_Fr.py` retains the original single-file behaviour and is kept unchanged.

---

## v1.0.0 — 2024-07-16

Initial public release.
- Interactive screening of MAFFT parameters across three predefined levels.
- Embedded simple scoring (conservation, gap penalty, complexity).
- French version `MAFFT_ScoreNGo_Fr.py` provided alongside the English version.
