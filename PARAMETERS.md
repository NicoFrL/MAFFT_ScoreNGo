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

1. Light: Quick screening with fewer parameter combinations.
2. Standard: Balanced screening with a moderate number of combinations.
3. Aggressive: Thorough screening with many parameter combinations.

The specific parameters tested depend on the chosen screening level.

## Scoring Algorithm

MAFFT_ScoreNGo evaluates each alignment using a custom scoring algorithm that considers the following factors:

1. Conservation Score: Measures the similarity of amino acids in each column of the alignment using the BLOSUM62 substitution matrix.
2. Gap Penalty: Penalizes the presence of gaps, with additional penalties for isolated gaps.
3. Complexity Score: Measures the diversity of amino acids in each column.

The final score is calculated as:
Final Score = Conservation Score - Gap Penalty + Complexity Score

For each set of MAFFT parameters, MAFFT_ScoreNGo runs the alignment and calculates this score. The parameter set that produces the highest final score is considered the optimal alignment strategy for the given dataset.

## Interpretation of Results

The output includes:
- Individual alignment files for each parameter combination
- A summary file (mafft_results_summary.txt) with scores for each combination
- A list of MAFFT commands used (mafft_commands.txt)
- Debug logs (debug_logs.txt) for troubleshooting

The combination with the highest final score is recommended as the best alignment strategy for your specific dataset. The top 13 alignments are ranked and presented in the results summary.
