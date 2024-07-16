# MAFFT_ScoreNGo Examples

This directory contains sample input and output files for the MAFFT_ScoreNGo tool. These files demonstrate how to use the tool and what kind of results to expect.

## Files in this directory

1. `sample_input.fasta`: An example input file containing unaligned sequences.
2. `sample_results_summary.txt`: Summary of the alignment results, including scores and rankings.
3. `sample_commands.txt`: The MAFFT commands used for each alignment attempt.
4. `sample_debug_logs.txt`: Debugging information from the tool's execution.
5. `sample_best_alignment.fasta`: The highest-scoring alignment produced by the tool.
6. `sample_second_best_alignment.fasta`: The second-highest scoring alignment.

## How to use these examples

1. To see how MAFFT_ScoreNGo works, you can run the tool on the sample input file:
```python3 MAFFT_ScoreNGo.py```
When prompted, select the `sample_input.fasta` file from this directory.

2. Compare your results with the sample output files provided here. Your results may vary slightly depending on your system and MAFFT version.

3. Examine the `sample_results_summary.txt` to see how different MAFFT parameters performed on this dataset.

4. Look at `sample_commands.txt` to see the exact MAFFT commands used for each alignment attempt.

5. If you encounter any issues, refer to `sample_debug_logs.txt` to see what kind of debugging information the tool provides.

Note: A full run of MAFFT_ScoreNGo typically generates multiple alignment files (one for each parameter combination tested). For space reasons, we've only included the top two alignments here. When you run the tool, expect to see more output files in your results directory.

## Using your own data

To use MAFFT_ScoreNGo with your own data, simply replace `sample_input.fasta` with your own FASTA file of unaligned sequences. Follow the same process as described above.

Remember, for best results, it's recommended to use a dataset of 200 sequences or less. If your dataset is larger, consider using the "Rand_NSamp_MyFasta.py" script (available soon) to extract a subset of sequences.

