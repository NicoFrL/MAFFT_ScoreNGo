# Parameters Reference

This document describes the MAFFT parameter combinations tested by `MAFFT_ScoreNGo.py` and the scoring metrics computed by `AlignmentScorerWithJalview.py`.

---

## Part 1 — MAFFT parameters tested by `MAFFT_ScoreNGo.py`

For a comprehensive description of each MAFFT option, see the official MAFFT documentation: [https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)

### Alignment strategies (always tested)

| Flag | Description |
| --- | --- |
| `--genafpair` | G-INS-i: generalized affine gap costs, suitable for sequences with conserved regions |
| `--localpair` | L-INS-i: local pairwise alignment, suitable for sequences with one alignable domain |
| `--globalpair` | G-INS-1/2: global pairwise alignment, suitable for sequences of similar length |

### Substitution matrices (always tested)

| Flag | Description |
| --- | --- |
| *(none)* | BLOSUM62 (MAFFT default) |
| `--bl 80` | BLOSUM80 (more stringent, for closer homologs) |

### Screening levels

Each level defines which gap, weight, and tree parameters are explored.

#### Light screening (`1`) — quick exploration

| Parameter group | Values tested |
| --- | --- |
| Gap open (`--op`) | 1.53 |
| Gap extension (`--ep`) | 0, 0.123 |
| Local opening (`--lop`, `--lexp`) | (default only) |
| Large gap (`--LOP`) | (default only) |
| Weight (`--weighti`) | (default only) |
| Tree iterations (`--retree`) | 2 |

Typical run: **13 combinations**.

#### Standard screening (`2`) — balanced (recommended)

| Parameter group | Values tested |
| --- | --- |
| Gap open (`--op`) | 1.53, 3.0 |
| Gap extension (`--ep`) | 0, 0.123, 0.5 |
| Local opening (`--lop`, `--lexp`) | default, `-1.0 / -0.5` |
| Large gap (`--LOP`) | default, `-8.00` |
| Weight (`--weighti`) | default, 4.0 |
| Tree iterations (`--retree`) | 2, 3 |

Typical run: **~337 combinations**. This is the level used in published analyses (see Lipp et al. 2026).

#### Aggressive screening (`3`) — thorough exploration

| Parameter group | Values tested |
| --- | --- |
| Gap open (`--op`) | 1.0, 1.53, 2.0, 3.0 |
| Gap extension (`--ep`) | 0, 0.123, 0.5, 1.0 |
| Local opening (`--lop`, `--lep`, `--lexp`) | default, `-1.0 / 0.0 / -0.1`, `-3.00 / 0.2 / -0.2` |
| Large gap (`--LOP`, `--LEXP`) | default, `-6.00 / 0.00`, `-8.00 / 0.1` |
| Weight (`--weighti`) | default, 1.0, 4.0 |
| Tree iterations (`--retree`) | 2, 3 |

Typical run: **>1000 combinations**. Expect significantly longer runtimes; recommended only for refined exploration once Standard has been run.

### Always-included parameters

Every combination is run with:
- `--maxiterate 1000` (iterative refinement)
- `--thread -1` (use all available cores)

An explicit E-INS-i combination (`--ep 0 --genafpair --maxiterate 1000 --thread -1`) is always added at the end.

### Custom parameters

After choosing a screening level, users may supply additional MAFFT parameter strings interactively. These are appended to the list of combinations to evaluate.

---

## Part 2 — Scoring metrics computed by `AlignmentScorerWithJalview.py`

The scorer evaluates every `alignment_*.fasta` file produced by `MAFFT_ScoreNGo.py` and computes the following metrics.

### Per-alignment metrics

| Metric | Description | Source |
| --- | --- | --- |
| **Sum-of-Pairs (SP)** | Total BLOSUM62 substitution score over all aligned pairs, excluding gap-gap | Classical metric |
| **Jalview Conservation (total)** | Sum of per-column conservation scores (0–11) | Livingstone & Barton (1993); Jalview |
| **Jalview Quality (total)** | Sum of per-column BLOSUM62-based quality scores | Jalview |
| **Hydrophobic Conservation** | Fraction of conserved hydrophobic positions across putatively hydrophobic columns | Custom |
| **Gap Coherence** | Measure of whether gaps in each sequence form coherent blocks rather than scattered insertions | Custom |
| **Parsimony-Informative Sites** | Fraction of columns containing at least two amino acids each present in at least two sequences (excluding gaps) | Phylogenetics standard |
| **Pairwise Identity** | Mean pairwise identity (%) across all sequence pairs, excluding gaps | Classical metric |

### Per-column outputs (in `per_position_scores.txt`)

For each alignment, the scorer also reports per-column **conservation** (0–11 scale) and **quality** values, allowing detailed inspection of which regions drive the alignment's score.

### Ranking methods

The seven metrics above are normalised and combined into three aggregate ranking methods:

#### 1. Weighted Average
Each normalised metric is multiplied by a phylogenetically-informed weight and summed:

| Metric | Weight |
| --- | --- |
| Conservation | 0.25 |
| Sum-of-Pairs | 0.20 |
| Quality | 0.20 |
| Parsimony-informative | 0.20 |
| Gap Coherence | 0.10 |
| Hydrophobic | 0.05 |

#### 2. TOPSIS
Technique for Order of Preference by Similarity to Ideal Solution. Each alignment is scored by its Euclidean distance to the ideal (best on every metric) and anti-ideal (worst on every metric) points in the normalised metric space.

#### 3. Borda Count
Each metric independently ranks all alignments; the per-metric ranks are summed. The alignment with the lowest total rank wins.

### Consensus detection

If at least two of the three ranking methods agree on the top alignment, the result is reported as a **consensus best**. Otherwise the output flags that the methods disagree, prompting the user to inspect the per-metric winners manually.

### Conservation algorithm (technical detail)

The Jalview conservation algorithm classifies each column based on agreement across 10 physicochemical property classes:

1. Hydrophobic
2. Small
3. Positive
4. Negative
5. Charged
6. Aromatic
7. Aliphatic
8. Tiny
9. Proline
10. Polar

For each column, residues present above a minimum count threshold contribute to determining whether each property is conserved (all residues share the property value) or not. The score is the number of conserved properties (0–10), with 11 reserved for full identity. Columns exceeding a gap threshold (default 25%) are scored as 0.

This is a direct re-implementation of the algorithm described in Livingstone & Barton (1993) and used by the Jalview alignment viewer.
