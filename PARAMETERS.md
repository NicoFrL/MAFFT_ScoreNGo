# MAFFT_ScoreNGo Parameters

This document provides detailed information about the parameters tested by MAFFT_ScoreNGo and the scoring algorithms used to evaluate alignments.

## MAFFT Parameters Tested

MAFFT_ScoreNGo tests various combinations of the following MAFFT parameters:

1. Alignment strategies:
   - `--genafpair`: FFT-NS-2 (Fast but rough)
   - `--localpair`: L-INS-i (Accurate but slow)
   - `--globalpair`: G-INS-i (Accurate but slow)

   Note: The combination of `--genafpair` with `--ep 0` and `--maxiterate 1000` corresponds to the E-INS-i strategy.

2. Substitution matrices:
   - BLOSUM62 (default)
   - BLOSUM80 (`--bl 80`)

3. Gap opening penalties (`--op`):
   - Light screening: 1.53
   - Standard screening: 1.53, 3.0
   - Aggressive screening: 1.0, 1.53, 2.0, 3.0

4. Gap extension penalties (`--ep`):
   - Light screening: 0, 0.123
   - Standard screening: 0, 0.123, 0.5
   - Aggressive screening: 0, 0.123, 0.5, 1.0

5. Local parameters (`--lop`, `--lep`, `--lexp`):
   - Light screening: None
   - Standard screening: None, `--lop -1.0 --lexp -0.5`
   - Aggressive screening: None, `--lop -1.0 --lep 0.0 --lexp -0.1`, `--lop -3.00 --lep 0.2 --lexp -0.2`

6. Large gap penalties (`--LOP`, `--LEXP`):
   - Light screening: None
   - Standard screening: None, `--LOP -8.00`
   - Aggressive screening: None, `--LOP -6.00 --LEXP 0.00`, `--LOP -8.00 --LEXP 0.1`

7. Weighting factor (`--weighti`):
   - Light screening: None
   - Standard and Aggressive screening: None, 4.0
   - Aggressive screening additionally: 1.0

8. Guide tree rebuilding (`--retree`):
   - Light screening: 2
   - Standard and Aggressive screening: 2, 3

Additional parameters:
- `--maxiterate 1000`: Maximum number of iterative refinement cycles (applied to all combinations)
- `--thread -1`: Use all available threads

## Screening Levels

MAFFT_ScoreNGo offers three screening levels:

1. Light: Quick screening with fewer parameter combinations (13 combinations).
2. Standard: Balanced screening with a moderate number of combinations (337 combinations).
3. Aggressive: Thorough screening with many parameter combinations (2497 combinations).

The specific parameters tested depend on the chosen screening level.

## Scoring Algorithm

MAFFT_ScoreNGo evaluates each alignment using a custom scoring algorithm that considers the following factors:

1. Weighted Conservation Score (WCS)
2. Gap Penalty (GP)
3. Complexity Score (CX)

### 1. Weighted Conservation Score (WCS)

The Weighted Conservation Score is calculated as follows:

    WCS = (1/N) * Σ(i=1 to N) [WCS'(i) / (L * (L-1))]

Where:
- WCS is the final Weighted Conservation Score for the entire alignment
- N is the total number of columns in the alignment
- L is the number of sequences in the alignment
- WCS'(i) is the raw Weighted Conservation Score for column i, calculated as:

    WCS'(i) = Σ(a,b) [count(a) * count(b) * BLOSUM62(a,b)]

Where:
- a and b are amino acids in column i
- count(x) is the number of occurrences of amino acid x in column i
- BLOSUM62(a,b) is the substitution score between a and b in the BLOSUM62 matrix

### 2. Gap Penalty (GP)

The Gap Penalty Score is calculated as follows:

    GP = (1/N) * Σ(i=1 to N) [(g_i / L) + 0.5 * I(i)]

Where:
- GP is the final Gap Penalty Score
- N is the total number of columns in the alignment
- g_i is the number of gaps in column i
- L is the number of sequences in the alignment
- I(i) is an indicator function that equals 1 if the gap in column i is isolated (no gap in the previous column), and 0 otherwise

### 3. Complexity Score (CX)

The Complexity Score is calculated as follows:

    CX = (1/N) * Σ(i=1 to N) (u_i / L)

Where:
- CX is the final Complexity Score
- N is the total number of columns in the alignment
- u_i is the number of unique amino acids in column i
- L is the number of sequences in the alignment

### Final Score

The Final Score is calculated as follows:

    Final Score = WCS - GP + CX

For each set of MAFFT parameters, MAFFT_ScoreNGo runs the alignment and calculates this score. The parameter set that produces the highest Final Score is considered the optimal alignment strategy for the given dataset.

## Interpretation of Results

The output includes:
- Individual alignment files for each parameter combination
- A summary file (mafft_results_summary.txt) with scores for each combination
- A list of MAFFT commands used (mafft_commands.txt)
- Debug logs (debug_logs.txt) for troubleshooting

The combination with the highest final score is recommended as the best alignment strategy for your specific dataset. The top 13 alignments are ranked and presented in the results summary.
