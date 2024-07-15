#This is a preliminary version of the code, it will run 108 alignments, large gap options are not coherent yet, and more options will be add later, notably to include EinsI methods.
#The script is currently in French, I'll translate it later.
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
        title="Sélectionnez le fichier FASTA d'entrée",
        filetypes=[("Fichiers FASTA", "*.fasta *.fa *.fna *.ffn *.faa *.frn")]
    )
    return file_path

def run_mafft(input_file, output_file, params):
    start_time = time.time()
    cmd = f"mafft {params} {input_file} > {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Erreur lors de l'exécution de MAFFT: {e}")
        print(f"Sortie d'erreur: {e.stderr}")
        return None
    end_time = time.time()
    return end_time - start_time

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
            
            # Score de conservation pondéré
            col_score = 0
            for a in char_counts:
                for b in char_counts:
                    if a != '-' and b != '-':
                        if (a, b) in blosum:
                            col_score += char_counts[a] * char_counts[b] * blosum[(a, b)]
                        elif (b, a) in blosum:
                            col_score += char_counts[a] * char_counts[b] * blosum[(b, a)]
            conservation_score += col_score / (len(alignment) * (len(alignment) - 1))
            
            # Pénalité de gap
            gap_count = char_counts.get('-', 0)
            if gap_count > 0:
                gap_penalty += gap_count / len(alignment)
                if i > 0 and '-' not in alignment[:, i-1]:
                    gap_penalty += 0.5  # Pénalité supplémentaire pour gap isolé
            
            # Complexité de colonne
            complexity_score += len(char_counts) / len(alignment)
        
        conservation_score /= alignment.get_alignment_length()
        gap_penalty /= alignment.get_alignment_length()
        complexity_score /= alignment.get_alignment_length()
        
        # Score final
        final_score = conservation_score - gap_penalty + complexity_score
        
        return final_score, conservation_score, gap_penalty, complexity_score
    except Exception as e:
        print(f"Erreur lors de l'évaluation de l'alignement: {e}")
        return None, None, None, None

def main():
    input_file = select_file()
    if not input_file:
        print("Aucun fichier sélectionné. Le programme se termine.")
        return

    strategies = ["--genafpair", "--localpair", "--globalpair"]
    matrices = ["", "--bl 80"]  # "" pour BLOSUM62 (défaut), "--bl 80" pour BLOSUM80
    gap_opens = ["--op 1.53", "--op 2.0", "--op 3.0"]
    gap_extensions = ["--ep 0.123", "--ep 0.5", "--ep 1.0"]
    large_gaps = ["", "--lop -2.00 --lep -0.1"]

    results = {}
    mafft_commands = {}

    combinations = list(itertools.product(strategies, matrices, gap_opens, gap_extensions, large_gaps))

    output_dir = os.path.join(os.path.dirname(input_file), "mafft_results")
    os.makedirs(output_dir, exist_ok=True)

    for i, (strategy, matrix, gap_open, gap_ext, large_gap) in enumerate(combinations, 1):
        params = f"{strategy} {matrix} {gap_open} {gap_ext} {large_gap} --maxiterate 1000 --thread -1"
        output_file = os.path.join(output_dir, f"alignment_{i}.fasta")
        
        print(f"\nExécution de la combinaison {i}/{len(combinations)}:")
        print(f"Paramètres : {params}")
        
        mafft_commands[i] = f"mafft {params} {input_file} > {output_file}"
        
        execution_time = run_mafft(input_file, output_file, params)
        if execution_time is None:
            continue
        
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
        print("Aucun résultat n'a été obtenu. Vérifiez les erreurs ci-dessus.")
        return

    print("\nRésultats :")
    print("=" * 50)
    for i, res in results.items():
        print(f"\nCombination {i}:")
        print(f"Paramètres : {res['params']}")
        print(f"Temps d'exécution : {res['execution_time']:.2f} secondes")
        print(f"Score final : {res['final_score']:.4f}")
        print(f"Score de conservation : {res['conservation_score']:.4f}")
        print(f"Pénalité de gap : {res['gap_penalty']:.4f}")
        print(f"Score de complexité : {res['complexity_score']:.4f}")

    best_combo = max(results, key=lambda x: results[x]['final_score'])
    print(f"\nLa combinaison avec le meilleur score final est : {best_combo}")
    print(f"Paramètres : {results[best_combo]['params']}")

    results_file = os.path.join(output_dir, "mafft_results_summary.txt")
    with open(results_file, 'w') as f:
        for i, res in results.items():
            f.write(f"Combination {i}:\n")
            f.write(f"Paramètres : {res['params']}\n")
            f.write(f"Temps d'exécution : {res['execution_time']:.2f} secondes\n")
            f.write(f"Score final : {res['final_score']:.4f}\n")
            f.write(f"Score de conservation : {res['conservation_score']:.4f}\n")
            f.write(f"Pénalité de gap : {res['gap_penalty']:.4f}\n")
            f.write(f"Score de complexité : {res['complexity_score']:.4f}\n\n")
        f.write(f"Meilleure combinaison : {best_combo}\n")
        f.write(f"Paramètres : {results[best_combo]['params']}\n")

    print(f"\nLes résultats détaillés ont été sauvegardés dans : {results_file}")

    # Écriture du fichier récapitulatif des commandes MAFFT
    mafft_commands_file = os.path.join(output_dir, "mafft_commands.txt")
    with open(mafft_commands_file, 'w') as f:
        for i, cmd in mafft_commands.items():
            f.write(f"Alignement {i}:\n{cmd}\n\n")

    print(f"\nLes commandes MAFFT utilisées ont été sauvegardées dans : {mafft_commands_file}")

if __name__ == "__main__":
    main()
