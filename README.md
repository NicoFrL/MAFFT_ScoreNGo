# MAFFT_ScoreNGo

A two-step pipeline to screen MAFFT alignment parameters and identify the optimal multiple sequence alignment for a given protein dataset.

Created by Nicolas-Frédéric Lipp, PhD.

---

## Overview

`MAFFT_ScoreNGo` performs systematic exploration of MAFFT alignment parameters to find the best alignment strategy for your dataset. The workflow is split into two complementary scripts:

1. **`MAFFT_ScoreNGo.py`** — generates alignments across many MAFFT parameter combinations, with support for batch processing of multiple FASTA files.
2. **`MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py`** — re-scores those alignments using rigorous, multi-metric evaluation (exact Jalview conservation and quality algorithms, Sum-of-Pairs, parsimony-informative sites, and more) to identify the optimal alignment.

This separation allows fast iteration on alignment generation while applying more sophisticated scoring methods as a distinct, reproducible step.

---

## Features

### MAFFT_ScoreNGo.py — Alignment generator
- Automated screening of MAFFT parameter combinations across three pre-defined levels (Light, Standard, Aggressive)
- Support for batch processing of multiple FASTA files in a single run
- Custom parameter injection
- Separate output directories per input file
- Comprehensive logging of MAFFT commands and execution

### MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py — Rigorous alignment scorer
- **Exact Jalview conservation** (per-column 0–11 scores based on 10 physicochemical property classes, AMAS-derived scheme)
- **Exact Jalview quality** (per-column BLOSUM62-based quality)
- Sum-of-Pairs score (BLOSUM62)
- Hydrophobic position conservation
- Gap coherence
- Parsimony-informative sites
- Pairwise identity
- Three ranking methods aggregating all metrics: Weighted Average, TOPSIS, Borda Count
- Consensus best-alignment detection across ranking methods
- Parallel processing for large alignment sets

***MAFFT stands for Multiple sequence Alignment using Fast Fourier Transform. More documentation can be found at [mafft.cbrc.jp](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html), [Katoh et al. (Nucleic Acids Res., 2002)](https://doi.org/10.1093%2Fnar%2Fgkf436), and [Katoh et al. (Brief. Bioinform., 2017)](https://doi.org/10.1093/bib/bbx108).*

---

## Parameters Overview

`MAFFT_ScoreNGo` tests various combinations of the following MAFFT parameters:

- Alignment strategies (`--genafpair`, `--localpair`, `--globalpair`)
- Substitution matrices (BLOSUM62, BLOSUM80)
- Gap opening penalties
- Gap extension penalties
- Large gap penalties
- Tree iteration count

For a detailed explanation of the parameters tested and the screening levels, see [PARAMETERS.md](./PARAMETERS.md).

---

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/NicoFrL/MAFFT_ScoreNGo.git
   cd MAFFT_ScoreNGo
   ```

2. Install the required Python packages:
   ```
   pip3 install -r requirements.txt
   ```

3. Ensure MAFFT is installed and accessible from your command line (see Dependencies section).

---

## Usage

### Step 1 — Generate alignments

```
python3 MAFFT_ScoreNGo.py
```

Follow the prompts to:
1. Select one or more input FASTA files (multi-selection supported)
2. Choose the screening level:
   - **Light** (`1`): Quick screening with fewer parameter combinations
   - **Standard** (`2`): Balanced screening (recommended for most use cases)
   - **Aggressive** (`3`): Thorough screening with many combinations
3. (Optional) Add custom MAFFT parameters to be tested alongside the predefined set
4. Confirm and let the script run

For each input file, a folder named `mafft_results_<basename>/` is created next to the input, containing:
- `alignment_<n>.fasta` — one file per parameter combination
- `mafft_commands.txt` — the exact MAFFT commands used
- `debug_logs.txt` — execution times and stderr output

### Step 2 — Score alignments and identify the best

```
python3 MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py
```

When prompted:
- Select the folder containing the alignments produced by Step 1 (or a parent folder containing multiple such subfolders — the script will detect them automatically).

The scorer will evaluate every `alignment_*.fasta` file and produce:
- `alignment_scores_EXACT_JALVIEW.txt` — ranked results with all metrics and consensus best alignment
- `per_position_scores.txt` — per-column conservation and quality scores for each alignment

---

## Output Files

### From `MAFFT_ScoreNGo.py`
| File | Description |
| --- | --- |
| `alignment_<n>.fasta` | One alignment per MAFFT parameter combination |
| `mafft_commands.txt` | All MAFFT commands executed |
| `debug_logs.txt` | Execution logs and MAFFT stderr output |

### From `MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py`
| File | Description |
| --- | --- |
| `alignment_scores_EXACT_JALVIEW.txt` | Full scoring results, ranking by three methods, and consensus best alignment |
| `per_position_scores.txt` | Per-column conservation (0–11) and quality scores |

---

## Dependencies

- Python 3.7 or later (3.9+ recommended)
- biopython >= 1.81
- numpy >= 1.24
- scipy >= 1.10
- tqdm >= 4.65
- tkinter (usually included with Python)
- MAFFT (must be installed separately and available in your system PATH)

### Verifying MAFFT installation

```
mafft --version
```

If the command is not found, install MAFFT:

- **macOS** (with [Homebrew](https://brew.sh)):
  ```
  brew install mafft
  ```
- **Ubuntu/Debian Linux**:
  ```
  sudo apt-get install mafft
  ```
- **Other systems**: see [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)

### Tkinter installation (if missing)

- **macOS** (Homebrew Python):
  ```
  brew install python-tk
  ```
- **Ubuntu/Debian**:
  ```
  sudo apt-get install python3-tk
  ```

---

## Performance Notes

- `MAFFT_ScoreNGo.py` runtime scales linearly with the number of parameter combinations. Light screening typically completes in under a minute for ~70 sequences of ~400 residues; Aggressive can take significantly longer.
- `MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py` uses parallel processing (one process per Performance core on Apple Silicon; bounded on other systems) and processes alignments in chunks of 50 to keep memory usage stable.
- For very large datasets (>200 sequences), consider running on a server or workstation rather than a laptop.
- Don't forget to "[caffeinate](https://www.theapplegeek.co.uk/blog/caffeinate)" your Mac during long runs (or use [systemd-inhibit](https://evanhahn.com/systemd-inhibit-alternative-to-macos-caffeinate/) on Linux).

---

## Examples

Example input files are provided in the `examples/` directory. To reproduce the example workflow:

1. Run `MAFFT_ScoreNGo.py` and select `examples/sample_input.fasta`.
2. Choose **Light** screening (fastest).
3. After completion, run `MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py` and point it to the resulting `mafft_results_sample_input/` folder.

---

## French Version

A French interactive version of the original (single-file, with simple scoring) is available as `MAFFT_ScoreNGo_Fr.py`. It is kept for legacy reasons; new users are encouraged to use the two-step workflow described above.

---

## License

This project is licensed under the MIT License — see the [LICENSE](./LICENSE) file for details.

---

## Contributing

Contributions are welcome. Please open an issue or submit a pull request on the GitHub repository.

---

## Support

If you encounter any problems or have questions, please open an issue on the [GitHub repository](https://github.com/NicoFrL/MAFFT_ScoreNGo/issues).

---

## Author

Nicolas-Frédéric Lipp, PhD
<https://github.com/NicoFrL>

---

## Acknowledgements

The exact Jalview conservation and quality algorithms in `MAFFT_AlignmentScorer_EXACT_JALVIEW-2.py` are re-implementations of the algorithms described in:

- Livingstone, C. D. & Barton, G. J. (1993) *Protein sequence alignments: a strategy for the hierarchical analysis of residue conservation.* Comput. Appl. Biosci. 9, 745–756.
- Waterhouse, A. M., Procter, J. B., Martin, D. M. A., Clamp, M. & Barton, G. J. (2009) *Jalview Version 2 — a multiple sequence alignment editor and analysis workbench.* Bioinformatics 25, 1189–1191.

This project was developed with assistance from AI language models for code structure and documentation. The scientific approach and core algorithm were designed and implemented by the author.
