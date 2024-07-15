import os
import subprocess
import time
from Bio import AlignIO
from Bio.Align import substitution_matrices
import itertools
import tkinter as tk
from tkinter import filedialog
from collections import Counter

def select_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select input FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn")]
    )
    return os.path.abspath(file_path)

def run_mafft(input_file, output_file, params):
    start_time = time.time()
    cmd = f"mafft {params} \"{input_file}\" > \"{output_file}\""
    try:
        result = subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE, text=True)
        end_time = time.time()
        return end_time - start_time, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error while executing MAFFT: {e}")
        return None, e.stderr

def evaluate_alignment(alignment_file):
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        
        conservation_score = 0
        gap_penalty = 0
        complexity_score = 0
        
        blosum = substitution_matrices.load("BLOSUM62")
        
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            char_counts = Counter(column)
            
            # Weighted conservation score
            col_score = 0
            for a in char_counts:
                for b in char_counts:
                    if a != '-' and b != '-' and a.isalpha() and b.isalpha():
                        if (a, b) in blosum:
                            col_score += char_counts[a] * char_counts[b] * blosum[(a, b)]
                        elif (b, a) in blosum:
                            col_score += char_counts[a] * char_counts[b] * blosum[(b, a)]
            conservation_score += col_score / (len(alignment) * (len(alignment) - 1))
            
            # Gap penalty
            gap_count = char_counts.get('-', 0)
            if gap_count > 0:
                gap_penalty += gap_count / len(alignment)
                if i > 0 and '-' not in alignment[:, i-1]:
                    gap_penalty += 0.5  # Additional penalty for isolated gap
            
            # Column complexity
            complexity_score += len(char_counts) / len(alignment)
        
        conservation_score /= alignment.get_alignment_length()
        gap_penalty /= alignment.get_alignment_length()
        complexity_score /= alignment.get_alignment_length()
        
        # Final score
        final_score = conservation_score - gap_penalty + complexity_score
        
        return final_score, conservation_score, gap_penalty, complexity_score
    except Exception as e:
        print(f"Error while evaluating the alignment: {e}")
        return None, None, None, None

def get_mafft_combinations(screening_level):
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
    
    # Explicit addition of E-INS-i
    combinations.append("--ep 0 --genafpair --maxiterate 1000 --thread -1")
    
    return combinations

def main():
    input_file = select_file()
    if not input_file:
        print("No file selected. The program is terminating.")
        return

    print("Choose the screening level:")
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
    
    combinations = get_mafft_combinations(level)
    
    # Request for custom parameter
    print("\nDo you want to add a specific parameter line to compare?")
    print("(Enter your parameters and press Enter. If nothing happens, press Enter a second time.)")
    print("(Leave empty and press Enter to ignore.)")
    
    custom_params = []
    while True:
        line = input("Custom parameters: ").strip()
        if line:
            custom_params.append(line)
        else:
            break
    
    custom_params = " ".join(custom_params)
    if custom_params:
        combinations.append(custom_params)
        print(f"Custom parameters added: {custom_params}")
    else:
        print("No custom parameters added.")
    
    print(f"\nTotal number of combinations: {len(combinations)}")
    
    while True:
        confirmation = input("Do you confirm this number of combinations? (Y/N): ").strip().lower()
        if confirmation == 'y':
            break
        elif confirmation == 'n':
            print("Operation cancelled by the user.")
            return
        else:
            print("Invalid response. Please answer 'Y' for yes or 'N' for no.")

    results = {}
    mafft_commands = {}
    debug_logs = []

    output_dir = os.path.join(os.path.dirname(input_file), "mafft_results")
    os.makedirs(output_dir, exist_ok=True)

    for i, params in enumerate(combinations, 1):
        output_file = os.path.join(output_dir, f"alignment_{i}.fasta")
        
        print(f"\nExecuting combination {i}/{len(combinations)}:")
        print(f"Parameters: {params}")
        
        debug_logs.append(f"\nExecuting combination {i}/{len(combinations)}:")
        debug_logs.append(f"Parameters: {params}")
        
        mafft_commands[i] = f"mafft {params} \"{input_file}\" > \"{output_file}\""
        
        execution_time, mafft_output = run_mafft(input_file, output_file, params)
        if execution_time is None:
            debug_logs.append(f"Execution failed for combination {i}")
            debug_logs.append(mafft_output)
            continue
        
        debug_logs.append(f"Execution time: {execution_time:.2f} seconds")
        debug_logs.append(f"MAFFT output: {mafft_output}")
        
        final_score, conservation_score, gap_penalty, complexity_score = evaluate_alignment(output_file)
        if final_score is None:
            continue
        
        results[i] = {
            "params": params,
            "execution_time": execution_time,
            "final_score": final_score,
            "conservation_score": conservation_score,
            "gap_penalty": gap_penalty,
            "complexity_score": complexity_score
        }

    if not results:
        print("No results were obtained. Check the errors above.")
        return

    # Sort results by descending final score
    sorted_results = sorted(results.items(), key=lambda x: x[1]['final_score'], reverse=True)

    print("\nRanking of the 13 best alignments:")
    print("=" * 50)
    for i, (combo, res) in enumerate(sorted_results[:13], 1):
        print(f"\nRank {i}:")
        print(f"Combination {combo}:")
        print(f"Parameters: {res['params']}")
        print(f"Final score: {res['final_score']:.4f}")
        print(f"Execution time: {res['execution_time']:.2f} seconds")
        print(f"Conservation score: {res['conservation_score']:.4f}")
        print(f"Gap penalty: {res['gap_penalty']:.4f}")
        print(f"Complexity score: {res['complexity_score']:.4f}")

    best_combo = sorted_results[0][0]
    print(f"\nThe combination with the best final score is: {best_combo}")
    print(f"Parameters: {results[best_combo]['params']}")

    results_file = os.path.join(output_dir, "mafft_results_summary.txt")
    with open(results_file, 'w') as f:
        f.write("Ranking of the 13 best alignments:\n")
        f.write("=" * 50 + "\n")
        for i, (combo, res) in enumerate(sorted_results[:13], 1):
            f.write(f"\nRank {i}:\n")
            f.write(f"Combination {combo}:\n")
            f.write(f"Parameters: {res['params']}\n")
            f.write(f"Final score: {res['final_score']:.4f}\n")
            f.write(f"Execution time: {res['execution_time']:.2f} seconds\n")
            f.write(f"Conservation score: {res['conservation_score']:.4f}\n")
            f.write(f"Gap penalty: {res['gap_penalty']:.4f}\n")
            f.write(f"Complexity score: {res['complexity_score']:.4f}\n")
        f.write(f"\nBest combination: {best_combo}\n")
        f.write(f"Parameters: {results[best_combo]['params']}\n")

    print(f"\nDetailed results have been saved in: {results_file}")

    # Write summary file of MAFFT commands
    mafft_commands_file = os.path.join(output_dir, "mafft_commands.txt")
    with open(mafft_commands_file, 'w') as f:
        for i, cmd in mafft_commands.items():
            f.write(f"Alignment {i}:\n{cmd}\n\n")

    print(f"\nThe MAFFT commands used have been saved in: {mafft_commands_file}")

    # Write debug log file
    debug_log_file = os.path.join(output_dir, "debug_logs.txt")
    with open(debug_log_file, 'w') as f:
        f.write("\n".join(debug_logs))

    print(f"\nDebug logs have been saved in: {debug_log_file}")

if __name__ == "__main__":
    main()
