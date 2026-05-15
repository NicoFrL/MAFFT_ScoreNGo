#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import time
import tkinter as tk
from tkinter import filedialog

def select_files():
    """Select one or multiple FASTA files"""
    root = tk.Tk()
    root.withdraw()
    file_paths = filedialog.askopenfilenames(
        title="Select input FASTA file(s)",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn")]
    )
    return [os.path.abspath(fp) for fp in file_paths]

def run_mafft(input_file, output_file, params):
    """Execute MAFFT with given parameters"""
    start_time = time.time()
    cmd = f"mafft {params} \"{input_file}\" > \"{output_file}\""
    try:
        result = subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE, text=True)
        end_time = time.time()
        return end_time - start_time, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error executing MAFFT: {e}")
        return None, e.stderr

def get_mafft_combinations(screening_level):
    """Generate MAFFT parameter combinations based on screening level"""
    strategies = ["--genafpair", "--localpair", "--globalpair"]
    matrices = ["", "--bl 80"]  # "" for BLOSUM62 (default), "--bl 80" for BLOSUM80

    if screening_level == "light":
        gap_opens = ["--op 1.53"]
        gap_extensions = ["--ep 0", "--ep 0.123"]
        local_params = [""]
        large_gaps = [""]
        weights = [""]
        retrees = ["--retree 2"]
    elif screening_level == "standard":
        gap_opens = ["--op 1.53", "--op 3.0"]
        gap_extensions = ["--ep 0", "--ep 0.123", "--ep 0.5"]
        local_params = ["", "--lop -1.0 --lexp -0.5"]
        large_gaps = ["", "--LOP -8.00"]
        weights = ["", "--weighti 4.0"]
        retrees = ["--retree 2", "--retree 3"]
    elif screening_level == "aggressive":
        gap_opens = ["--op 1.0", "--op 1.53", "--op 2.0", "--op 3.0"]
        gap_extensions = ["--ep 0", "--ep 0.123", "--ep 0.5", "--ep 1.0"]
        local_params = ["", "--lop -1.0 --lep 0.0 --lexp -0.1", "--lop -3.00 --lep 0.2 --lexp -0.2"]
        large_gaps = ["", "--LOP -6.00 --LEXP 0.00", "--LOP -8.00 --LEXP 0.1"]
        weights = ["", "--weighti 1.0", "--weighti 4.0"]
        retrees = ["--retree 2", "--retree 3"]

    combinations = []
    for strategy in strategies:
        for matrix in matrices:
            for gap_open in gap_opens:
                for gap_ext in gap_extensions:
                    for weight in weights:
                        for retree in retrees:
                            base_params = f"{strategy} {matrix} {gap_open} {gap_ext} {weight} {retree} --maxiterate 1000 --thread -1"

                            if strategy in ["--genafpair", "--localpair"]:
                                for local_param in local_params:
                                    if strategy == "--genafpair":
                                        for large_gap in large_gaps:
                                            params = f"{base_params} {local_param} {large_gap}"
                                            combinations.append(params)
                                    else:
                                        params = f"{base_params} {local_param}"
                                        combinations.append(params)
                            else:
                                combinations.append(base_params)

    # Explicitly add E-INS-i
    combinations.append("--ep 0 --genafpair --maxiterate 1000 --thread -1")

    return combinations

def process_single_file(input_file, level, custom_params, skip_confirmation=False):
    """Process a single FASTA file with all parameter combinations"""
    combinations = get_mafft_combinations(level)

    if custom_params:
        combinations.append(custom_params)
        print(f"Custom parameters added: {custom_params}")

    print(f"\nTotal number of combinations for {os.path.basename(input_file)}: {len(combinations)}")

    # Ask for confirmation only for the first file
    if not skip_confirmation:
        while True:
            confirmation = input("Confirm this number of combinations? (Y/N): ").strip().lower()
            if confirmation == 'y':
                break
            elif confirmation == 'n':
                print("Processing cancelled for this file.")
                return False
            else:
                print("Invalid response. Please answer 'Y' for yes or 'N' for no.")

    # Create output directory with input file name
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = os.path.join(os.path.dirname(input_file), f"mafft_results_{base_name}")
    os.makedirs(output_dir, exist_ok=True)

    mafft_commands = {}
    debug_logs = []

    for i, params in enumerate(combinations, 1):
        output_file = os.path.join(output_dir, f"alignment_{i}.fasta")

        print(f"\nRunning combination {i}/{len(combinations)}:")
        print(f"Parameters: {params}")

        debug_logs.append(f"\nRunning combination {i}/{len(combinations)}:")
        debug_logs.append(f"Parameters: {params}")

        mafft_commands[i] = f"mafft {params} \"{input_file}\" > \"{output_file}\""

        execution_time, mafft_output = run_mafft(input_file, output_file, params)
        if execution_time is None:
            debug_logs.append(f"Execution failed for combination {i}")
            debug_logs.append(mafft_output)
            continue

        debug_logs.append(f"Execution time: {execution_time:.2f} seconds")
        debug_logs.append(f"MAFFT output: {mafft_output}")

    # Save MAFFT commands
    mafft_commands_file = os.path.join(output_dir, "mafft_commands.txt")
    with open(mafft_commands_file, 'w', encoding='utf-8') as f:
        f.write(f"MAFFT commands for: {os.path.basename(input_file)}\n")
        f.write("=" * 50 + "\n\n")
        for i, cmd in mafft_commands.items():
            f.write(f"Alignment {i}:\n{cmd}\n\n")

    print(f"\nMAFFT commands saved to: {mafft_commands_file}")

    # Save debug logs
    debug_log_file = os.path.join(output_dir, "debug_logs.txt")
    with open(debug_log_file, 'w', encoding='utf-8') as f:
        f.write(f"Debug logs for: {os.path.basename(input_file)}\n")
        f.write("=" * 50 + "\n")
        f.write("\n".join(debug_logs))

    print(f"Debug logs saved to: {debug_log_file}")
    print(f"\nAlignments saved to: {output_dir}")
    print("➡  To score and identify the best alignment, use:")
    print("   AlignmentScorerWithJalview.py")
    print("   (available in this repository: https://github.com/NicoFrL/MAFFT_ScoreNGo)")

    return True

def main():
    """Main function to coordinate batch processing"""
    input_files = select_files()
    if not input_files:
        print("No file selected. Exiting.")
        return

    print(f"\n{len(input_files)} file(s) selected:")
    for file in input_files:
        print(f"  - {os.path.basename(file)}")

    # Ask for screening level once for all files
    print("\nChoose screening level:")
    print("1. Light")
    print("2. Standard")
    print("3. Aggressive")

    choice = input("Enter your choice (1, 2 or 3): ")

    if choice == "1":
        level = "light"
    elif choice == "2":
        level = "standard"
    elif choice == "3":
        level = "aggressive"
    else:
        print("Invalid choice. Using standard level.")
        level = "standard"

    # Ask for custom parameters once
    print("\nWould you like to add custom parameters to compare?")
    print("(Enter your parameters and press Enter. If nothing happens, press Enter a second time.)")
    print("(Leave empty and press Enter to skip.)")

    custom_params = []
    while True:
        line = input("Custom parameters: ").strip()
        if line:
            custom_params.append(line)
        else:
            break

    custom_params = " ".join(custom_params)

    # Process each file
    successful_files = []
    failed_files = []

    for file_index, input_file in enumerate(input_files, 1):
        print(f"\n{'='*60}")
        print(f"Processing file {file_index}/{len(input_files)}: {os.path.basename(input_file)}")
        print(f"{'='*60}")

        # Only ask for confirmation on the first file
        skip_confirmation = (file_index > 1)

        success = process_single_file(input_file, level, custom_params, skip_confirmation)

        if success:
            successful_files.append(input_file)
        else:
            failed_files.append(input_file)

    # Final summary
    print(f"\n{'='*60}")
    print("BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    print(f"\nFiles processed successfully ({len(successful_files)}):")
    for file in successful_files:
        print(f"  ✓ {os.path.basename(file)}")

    if failed_files:
        print(f"\nFailed or cancelled files ({len(failed_files)}):")
        for file in failed_files:
            print(f"  ✗ {os.path.basename(file)}")

    print(f"\nProcessing complete!")
    print(f"Total: {len(successful_files)}/{len(input_files)} files processed successfully")
    print(f"\n{'='*60}")
    print("NEXT STEP: ALIGNMENT SCORING")
    print(f"{'='*60}")
    print("To identify the best alignment among those generated,")
    print("use AlignmentScorerWithJalview.py,")
    print("available in this repository:")
    print("https://github.com/NicoFrL/MAFFT_ScoreNGo")

if __name__ == "__main__":
    main()
