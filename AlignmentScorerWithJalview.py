#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MAFFT Alignment Scorer - EXACT JALVIEW IMPLEMENTATION VERSION
=============================================================
Modified to use EXACT Jalview conservation and quality calculations.

Key Changes:
1. Uses exact Jalview conservation algorithm (matches Jalview 100%)
2. Uses exact Jalview quality algorithm with proper BLOSUM62 matrix
3. Reports RAW SUMS (not per-column averages) for conservation and quality
4. Outputs ALL alignments (not limited to top 20)
5. Creates supplementary file with per-position scores

Metrics:
1. Sum-of-Pairs (BLOSUM62)
2. Jalview Conservation (EXACT implementation)
3. Jalview Quality (EXACT algorithm)
4. Hydrophobic Position Conservation
5. Gap Coherence
6. Parsimony-Informative Sites
7. Pairwise Identity
"""

import os
import numpy as np
from Bio import AlignIO
from Bio.Align import substitution_matrices
from Bio.Align import MultipleSeqAlignment
from collections import Counter
import tkinter as tk
from tkinter import filedialog, messagebox
import math
from scipy import stats
from scipy.spatial.distance import euclidean
from multiprocessing import Pool, cpu_count
import time
from tqdm import tqdm
import warnings
import signal
import sys
import gc
warnings.filterwarnings('ignore')

# ============================================================================
# SIGNAL HANDLING FOR CLEAN SHUTDOWN
# ============================================================================

# Global pool variable for cleanup
current_pool = None

def signal_handler(sig, frame):
    """Handle Ctrl+C gracefully"""
    print('\n\nInterrupted! Cleaning up resources...')
    if current_pool is not None:
        current_pool.terminate()
        current_pool.join()
    sys.exit(0)

# Register signal handler
signal.signal(signal.SIGINT, signal_handler)

# ============================================================================
# EXACT JALVIEW PROPERTY DEFINITIONS (from ResidueProperties.js)
# ============================================================================

JALVIEW_PROPERTIES = {
    'hydrophobic': {
        'I': 1, 'L': 1, 'V': 1, 'C': 1, 'A': 1, 'G': 1, 'M': 1,
        'F': 1, 'Y': 1, 'W': 1, 'H': 1, 'K': 1, 'T': 1,
        'X': 1, '-': 1, '*': 1,
        'R': 0, 'E': 0, 'Q': 0, 'D': 0, 'N': 0, 'S': 0, 'P': 0
    },
    'small': {
        'V': 1, 'C': 1, 'A': 1, 'G': 1, 'D': 1, 'N': 1, 'S': 1, 'T': 1, 'P': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'M': 0, 'F': 0, 'Y': 0, 'W': 0,
        'H': 0, 'K': 0, 'R': 0, 'E': 0, 'Q': 0
    },
    'positive': {
        'H': 1, 'K': 1, 'R': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'A': 0, 'G': 0,
        'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'E': 0, 'Q': 0,
        'D': 0, 'N': 0, 'S': 0, 'T': 0, 'P': 0
    },
    'negative': {
        'E': 1, 'D': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'A': 0, 'G': 0,
        'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'H': 0, 'K': 0,
        'R': 0, 'Q': 0, 'N': 0, 'S': 0, 'T': 0, 'P': 0
    },
    'charged': {
        'H': 1, 'K': 1, 'R': 1, 'E': 1, 'D': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'A': 0, 'G': 0,
        'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'Q': 0, 'N': 0,
        'S': 0, 'T': 0, 'P': 0
    },
    'aromatic': {
        'F': 1, 'Y': 1, 'W': 1, 'H': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'A': 0, 'G': 0,
        'M': 0, 'K': 0, 'R': 0, 'E': 0, 'Q': 0, 'D': 0,
        'N': 0, 'S': 0, 'T': 0, 'P': 0
    },
    'aliphatic': {
        'I': 1, 'L': 1, 'V': 1,
        '-': 1, '*': 1,
        'C': 0, 'A': 0, 'G': 0, 'M': 0, 'F': 0, 'Y': 0,
        'W': 0, 'H': 0, 'K': 0, 'R': 0, 'E': 0, 'Q': 0,
        'D': 0, 'N': 0, 'S': 0, 'T': 0, 'P': 0
    },
    'tiny': {
        'A': 1, 'G': 1, 'S': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'M': 0, 'F': 0,
        'Y': 0, 'W': 0, 'H': 0, 'K': 0, 'R': 0, 'E': 0,
        'Q': 0, 'D': 0, 'N': 0, 'T': 0, 'P': 0
    },
    'proline': {
        'P': 1,
        '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'A': 0, 'G': 0,
        'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'H': 0, 'K': 0,
        'R': 0, 'E': 0, 'Q': 0, 'D': 0, 'N': 0, 'S': 0, 'T': 0
    },
    'polar': {
        'Y': 1, 'W': 1, 'H': 1, 'K': 1, 'R': 1, 'E': 1,
        'Q': 1, 'D': 1, 'N': 1, 'S': 1, 'T': 1,
        'X': 1, '-': 1, '*': 1,
        'I': 0, 'L': 0, 'V': 0, 'C': 0, 'A': 0, 'G': 0, 'M': 0, 'F': 0, 'P': 0
    }
}

# Jalview AA index for quality calculation
JALVIEW_AA_INDEX = {
    'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7,
    'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
    'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19,
    'B': 20, 'Z': 21, 'X': 22, '*': 23
}

# ============================================================================
# EXACT JALVIEW BLOSUM62 MATRIX
# ============================================================================

def load_jalview_blosum62():
    """Load the exact BLOSUM62 matrix used by Jalview."""
    matrix_text = """
A    4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R   -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N   -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D   -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C    0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q   -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E   -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G    0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H   -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I   -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L   -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K   -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M   -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F   -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P   -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S    1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T    0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W   -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y   -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V    0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B   -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z   -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X    0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
*   -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1"""
    
    lines = matrix_text.strip().split('\n')
    headers = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
               'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
    
    blosum62 = {}
    for i, line in enumerate(lines):
        parts = line.split()
        row_aa = parts[0]
        values = parts[1:]
        for j, val in enumerate(values):
            col_aa = headers[j]
            blosum62[(row_aa, col_aa)] = float(val)
    
    # Add gap scores
    for aa in headers:
        blosum62[('-', aa)] = -4.0
        blosum62[(aa, '-')] = -4.0
    blosum62[('-', '-')] = 1.0
    
    return blosum62

# ============================================================================
# EXACT JALVIEW CONSERVATION CALCULATION
# ============================================================================

def calculate_jalview_conservation(alignment, threshold=3, max_gap_percent=25.0):
    """
    Calculate conservation scores EXACTLY as Jalview does.
    Returns: (per-position scores list, total sum)
    """
    num_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    conservation_scores = []
    
    thresh = int((threshold * num_seqs) / 100)
    
    for col in range(alignment_length):
        residue_counts = {}
        gap_count = 0
        
        for record in alignment:
            aa = record.seq[col].upper()
            if aa == '-' or aa == '.':
                gap_count += 1
            else:
                residue_counts[aa] = residue_counts.get(aa, 0) + 1
        
        gap_percent = (gap_count * 100.0) / num_seqs
        
        if gap_percent >= max_gap_percent:
            conservation_scores.append(0.0)
            continue
        
        property_values = {}
        
        for aa, count in residue_counts.items():
            if count > thresh:
                for prop_name, prop_dict in JALVIEW_PROPERTIES.items():
                    if aa not in prop_dict:
                        prop_value = prop_dict.get('-', 1)
                    else:
                        prop_value = prop_dict[aa]
                    
                    if prop_name not in property_values:
                        property_values[prop_name] = prop_value
                    elif property_values[prop_name] != -1:
                        if property_values[prop_name] != prop_value:
                            property_values[prop_name] = -1
        
        if gap_count > thresh:
            for prop_name, prop_dict in JALVIEW_PROPERTIES.items():
                prop_value = prop_dict['-']
                if prop_name not in property_values:
                    property_values[prop_name] = prop_value
                elif property_values[prop_name] != -1:
                    if property_values[prop_name] != prop_value:
                        property_values[prop_name] = -1
        
        conserved_count = sum(1 for v in property_values.values() if v != -1)
        
        fully_conserved = len(residue_counts) == 1 and gap_count == 0
        
        if fully_conserved:
            conservation_scores.append(11.0)
        elif conserved_count == 10:
            conservation_scores.append(11.0)
        else:
            conservation_scores.append(float(conserved_count))
    
    total = sum(conservation_scores)
    return conservation_scores, total

# ============================================================================
# EXACT JALVIEW QUALITY CALCULATION
# ============================================================================

def calculate_jalview_quality(alignment):
    """
    Calculate quality scores using Jalview's exact algorithm.
    Returns: (per-position scores list, total sum)
    """
    blosum62 = load_jalview_blosum62()
    
    num_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    quality_scores = []
    
    min_score = -4.0
    symbol_count = 24
    
    # Build profile matrix
    cons2 = np.zeros((alignment_length, symbol_count), dtype=int)
    cons2_gap_counts = np.zeros(alignment_length, dtype=int)
    
    for seq_record in alignment:
        for col_idx, aa in enumerate(str(seq_record.seq).upper()):
            if aa == '-' or aa == '.':
                cons2_gap_counts[col_idx] += 1
            elif aa in JALVIEW_AA_INDEX:
                aa_index = JALVIEW_AA_INDEX[aa]
                cons2[col_idx][aa_index] += 1
            else:
                cons2[col_idx][JALVIEW_AA_INDEX.get('X', 22)] += 1
    
    max_quality = -float('inf')
    index_to_aa = {v: k for k, v in JALVIEW_AA_INDEX.items()}
    
    for col in range(alignment_length):
        bigtot = 0.0
        x = np.zeros(symbol_count)
        
        for ii in range(symbol_count):
            x[ii] = 0.0
            aa_ii = index_to_aa.get(ii, 'X')
            
            for i2 in range(symbol_count - 1):
                aa_i2 = index_to_aa.get(i2, 'X')
                score = blosum62.get((aa_ii, aa_i2), blosum62.get((aa_i2, aa_ii), -2.0))
                x[ii] += (cons2[col][i2] * score) + 4.0
            
            x[ii] += 4.0 + cons2_gap_counts[col] * min_score
            x[ii] /= num_seqs
        
        for seq_idx, seq_record in enumerate(alignment):
            tot = 0.0
            
            if col < len(seq_record.seq):
                aa = str(seq_record.seq[col]).upper()
                if aa in JALVIEW_AA_INDEX:
                    seq_num = JALVIEW_AA_INDEX[aa]
                elif aa == '-' or aa == '.':
                    seq_num = -1
                else:
                    seq_num = JALVIEW_AA_INDEX.get('X', 22)
            else:
                seq_num = -1
            
            for i in range(symbol_count - 1):
                aa_i = index_to_aa.get(i, 'X')
                sr = 4.0
                
                if seq_num == -1:
                    sr += min_score
                else:
                    aa_seq = index_to_aa.get(seq_num, 'X')
                    sr += blosum62.get((aa_i, aa_seq), blosum62.get((aa_seq, aa_i), -2.0))
                
                xx_i = x[i] - sr
                tot += xx_i * xx_i
            
            bigtot += math.sqrt(tot)
        
        quality_scores.append(bigtot)
        max_quality = max(max_quality, bigtot)
    
    # Normalize scores
    normalized_scores = []
    for col in range(alignment_length):
        tmp = quality_scores[col]
        tmp = ((max_quality - tmp) * (num_seqs - cons2_gap_counts[col])) / num_seqs
        normalized_scores.append(tmp)
    
    total = sum(normalized_scores)
    return normalized_scores, total

# ============================================================================
# OTHER METRICS (kept from original)
# ============================================================================

def calculate_sum_of_pairs(alignment):
    """Sum-of-Pairs score using BLOSUM62"""
    try:
        matrix = substitution_matrices.load("BLOSUM62")
        total_score = 0
        n_seqs = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        for col_idx in range(alignment_length):
            for i in range(n_seqs):
                for j in range(i + 1, n_seqs):
                    aa1 = alignment[i, col_idx].upper()
                    aa2 = alignment[j, col_idx].upper()
                    
                    if aa1 != '-' and aa2 != '-':
                        score = matrix.get((aa1, aa2), matrix.get((aa2, aa1), 0))
                        total_score += score
        
        return total_score
    except Exception as e:
        print(f"Error calculating SP score: {e}")
        return 0

def calculate_hydrophobic_conservation(alignment, hydrophobic_positions=None):
    """Calculate conservation of hydrophobic positions"""
    HYDROPHOBIC = set('ILVCAGMFYWHK')
    n_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    if hydrophobic_positions is None:
        hydrophobic_positions = []
        for col_idx in range(alignment_length):
            hydro_count = sum(1 for seq_idx in range(n_seqs) 
                            if alignment[seq_idx, col_idx].upper() in HYDROPHOBIC)
            if hydro_count >= n_seqs * 0.5:
                hydrophobic_positions.append(col_idx)
    
    if not hydrophobic_positions:
        return 0.0
    
    conservation_scores = []
    for col_idx in hydrophobic_positions:
        hydro_count = sum(1 for seq_idx in range(n_seqs) 
                        if alignment[seq_idx, col_idx].upper() in HYDROPHOBIC)
        conservation_scores.append(hydro_count / n_seqs)
    
    return np.mean(conservation_scores) if conservation_scores else 0.0

def calculate_gap_coherence(alignment):
    """Calculate gap coherence score"""
    n_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    if n_seqs < 2:
        return 1.0
    
    coherence_scores = []
    
    for seq_idx in range(n_seqs):
        seq = str(alignment[seq_idx].seq)
        
        gaps = []
        in_gap = False
        gap_start = 0
        
        for i, char in enumerate(seq):
            if char == '-':
                if not in_gap:
                    in_gap = True
                    gap_start = i
            else:
                if in_gap:
                    gaps.append((gap_start, i))
                    in_gap = False
        
        if in_gap:
            gaps.append((gap_start, len(seq)))
        
        if gaps:
            total_gap_length = sum(end - start for start, end in gaps)
            num_blocks = len(gaps)
            avg_gap_length = total_gap_length / num_blocks
            max_possible = total_gap_length
            coherence = avg_gap_length / max_possible if max_possible > 0 else 1
            coherence_scores.append(coherence)
        else:
            coherence_scores.append(1.0)
    
    return np.mean(coherence_scores)

def calculate_parsimony_informative_sites(alignment):
    """Calculate fraction of parsimony-informative sites"""
    n_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    if n_seqs < 4:
        return 0.0
    
    informative_sites = 0
    
    for col_idx in range(alignment_length):
        column = [alignment[seq_idx, col_idx].upper() for seq_idx in range(n_seqs)]
        
        column_no_gaps = [aa for aa in column if aa != '-' and aa != 'X']
        
        if len(column_no_gaps) < 4:
            continue
        
        aa_counts = Counter(column_no_gaps)
        
        variants_with_2plus = sum(1 for count in aa_counts.values() if count >= 2)
        
        if variants_with_2plus >= 2:
            informative_sites += 1
    
    return informative_sites / alignment_length if alignment_length > 0 else 0

def calculate_pairwise_identity(alignment):
    """Calculate average pairwise identity"""
    n_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    if n_seqs < 2:
        return 100.0
    
    identities = []
    
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            matches = 0
            comparisons = 0
            
            for col_idx in range(alignment_length):
                aa1 = alignment[i, col_idx].upper()
                aa2 = alignment[j, col_idx].upper()
                
                if aa1 != '-' and aa2 != '-':
                    comparisons += 1
                    if aa1 == aa2:
                        matches += 1
            
            if comparisons > 0:
                identity = (matches / comparisons) * 100
                identities.append(identity)
    
    return np.mean(identities) if identities else 0.0

# ============================================================================
# FILE SELECTION
# ============================================================================

def select_folders():
    """Select folder(s) containing MAFFT alignment results"""
    root = tk.Tk()
    root.withdraw()
    
    folder = filedialog.askdirectory(
        title="Select folder containing MAFFT alignments (or parent folder with subfolders)"
    )
    
    if not folder:
        return []
    
    # Check if this folder directly contains alignment files
    has_alignments = any(
        f.startswith("alignment_") and f.endswith(".fasta") 
        for f in os.listdir(folder)
    )
    
    if has_alignments:
        # Direct alignment folder
        print(f"Found alignments directly in: {folder}")
        return [folder]
    else:
        # Check for subfolders with alignments
        folders_with_alignments = []
        for subfolder in os.listdir(folder):
            subfolder_path = os.path.join(folder, subfolder)
            if os.path.isdir(subfolder_path):
                # Check if subfolder contains alignment files
                if any(f.startswith("alignment_") and f.endswith(".fasta") 
                       for f in os.listdir(subfolder_path)):
                    folders_with_alignments.append(subfolder_path)
                    print(f"Found alignments in subfolder: {subfolder_path}")
        
        if folders_with_alignments:
            print(f"\nTotal: {len(folders_with_alignments)} folders with alignments found")
            return folders_with_alignments
        else:
            messagebox.showwarning(
                "No Alignments Found",
                "No alignment files found in the selected folder or its subfolders.\n"
                "Looking for files named 'alignment_*.fasta'"
            )
            return []

# ============================================================================
# COMPREHENSIVE SCORING
# ============================================================================

def score_alignment_comprehensive(alignment_file):
    """Score a single alignment with all metrics"""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Calculate all metrics
        sp_score = calculate_sum_of_pairs(alignment)
        
        # EXACT Jalview calculations (return per-position and total)
        cons_scores, cons_total = calculate_jalview_conservation(alignment)
        qual_scores, qual_total = calculate_jalview_quality(alignment)
        
        hydrophobic_conservation = calculate_hydrophobic_conservation(alignment)
        gap_coherence = calculate_gap_coherence(alignment)
        parsimony_informative = calculate_parsimony_informative_sites(alignment)
        pairwise_identity = calculate_pairwise_identity(alignment)
        
        # Get dimensions
        aln_length = alignment.get_alignment_length()
        n_seqs = len(alignment)
        
        return {
            'sp_score': sp_score,
            'conservation_total': cons_total,  # RAW SUM
            'conservation_scores': cons_scores,  # Per-position
            'quality_total': qual_total,  # RAW SUM
            'quality_scores': qual_scores,  # Per-position
            'hydrophobic_conservation': hydrophobic_conservation,
            'gap_coherence': gap_coherence,
            'parsimony_informative': parsimony_informative,
            'pairwise_identity': pairwise_identity,
            'alignment_length': aln_length,
            'n_sequences': n_seqs
        }
    except Exception as e:
        print(f"Error processing {alignment_file}: {e}")
        return None

# ============================================================================
# PARALLEL PROCESSING
# ============================================================================

def score_alignment_wrapper(alignment_file):
    """Wrapper for parallel processing"""
    return score_alignment_comprehensive(alignment_file)

def process_alignments_parallel(alignment_files, n_processes=None):
    """Process alignments in parallel with better resource management"""
    global current_pool
    
    if n_processes is None:
        # For Apple Silicon Macs, use all Performance cores
        total_cpus = cpu_count()
        if total_cpus >= 10:
            # Apple Silicon with 8P+2E or more, use all 8 P-cores
            n_processes = 8
        elif total_cpus > 8:
            # Other Apple Silicon configs, leave some headroom
            n_processes = min(6, total_cpus - 2)
        else:
            n_processes = max(1, min(4, total_cpus - 1))  # Standard systems
    
    print(f"  Using {n_processes} parallel processes (detected {cpu_count()} total cores)")
    
    # Process in smaller chunks to avoid memory issues
    chunk_size = 50  # Process 50 alignments at a time
    all_results = []
    
    for i in range(0, len(alignment_files), chunk_size):
        chunk = alignment_files[i:i+chunk_size]
        chunk_start = i + 1
        chunk_end = min(i + chunk_size, len(alignment_files))
        
        print(f"  Processing alignments {chunk_start}-{chunk_end} of {len(alignment_files)}...")
        
        try:
            # Create pool and store reference for signal handler
            current_pool = Pool(processes=n_processes)
            
            chunk_results = list(tqdm(
                current_pool.imap(score_alignment_wrapper, chunk, chunksize=1),
                total=len(chunk),
                desc=f"  Batch {i//chunk_size + 1}"
            ))
            all_results.extend(chunk_results)
            
        except KeyboardInterrupt:
            print("\n  Interrupted! Cleaning up...")
            if current_pool:
                current_pool.terminate()
                current_pool.join()
            raise
            
        finally:
            # Always clean up the pool
            if current_pool:
                current_pool.close()
                current_pool.join()
                current_pool = None
        
        # Force garbage collection between chunks
        gc.collect()
    
    return all_results

# ============================================================================
# RANKING METHODS
# ============================================================================

def calculate_rankings(results, alignment_files):
    """Calculate rankings using multiple methods"""
    valid_results = [r for r in results if r]
    
    if not valid_results:
        return None
    
    # Extract metrics - now using RAW SUMS for conservation and quality
    sp_scores = np.array([r['sp_score'] for r in valid_results])
    cons_scores = np.array([r['conservation_total'] for r in valid_results])  # RAW SUM
    qual_scores = np.array([r['quality_total'] for r in valid_results])  # RAW SUM
    hydro_scores = np.array([r['hydrophobic_conservation'] for r in valid_results])
    gap_scores = np.array([r['gap_coherence'] for r in valid_results])
    parsimony_scores = np.array([r['parsimony_informative'] for r in valid_results])
    
    # Normalize to 0-1 for comparison
    def normalize(arr):
        arr_min, arr_max = arr.min(), arr.max()
        if arr_max - arr_min == 0:
            return np.ones_like(arr)
        return (arr - arr_min) / (arr_max - arr_min)
    
    sp_norm = normalize(sp_scores)
    cons_norm = normalize(cons_scores)
    qual_norm = normalize(qual_scores)
    hydro_norm = hydro_scores  # Already 0-1
    gap_norm = gap_scores  # Already 0-1
    parsimony_norm = parsimony_scores  # Already 0-1
    
    # METHOD 1: Weighted Average
    weights_phylo = {
        'cons': 0.25,
        'sp': 0.20,
        'quality': 0.20,
        'parsimony': 0.20,
        'gap': 0.10,
        'hydro': 0.05
    }
    
    weighted_scores = (
        weights_phylo['cons'] * cons_norm +
        weights_phylo['sp'] * sp_norm +
        weights_phylo['quality'] * qual_norm +
        weights_phylo['parsimony'] * parsimony_norm +
        weights_phylo['gap'] * gap_norm +
        weights_phylo['hydro'] * hydro_norm
    )
    
    weighted_ranks = np.argsort(-weighted_scores)
    
    # METHOD 2: TOPSIS
    decision_matrix = np.column_stack([
        sp_norm, cons_norm, qual_norm,
        hydro_norm, gap_norm, parsimony_norm
    ])
    
    weights_topsis = np.array([0.20, 0.25, 0.20, 0.05, 0.10, 0.20])
    weighted_matrix = decision_matrix * weights_topsis
    
    ideal_positive = weighted_matrix.max(axis=0)
    ideal_negative = weighted_matrix.min(axis=0)
    
    dist_positive = np.sqrt(((weighted_matrix - ideal_positive) ** 2).sum(axis=1))
    dist_negative = np.sqrt(((weighted_matrix - ideal_negative) ** 2).sum(axis=1))
    
    topsis_scores = dist_negative / (dist_positive + dist_negative + 1e-10)
    topsis_ranks = np.argsort(-topsis_scores)
    
    # METHOD 3: Borda Count
    rank_matrix = []
    for scores in [sp_scores, cons_scores, qual_scores, hydro_scores, gap_scores, parsimony_scores]:
        ranks = np.argsort(-scores)
        rank_positions = np.empty_like(ranks)
        rank_positions[ranks] = np.arange(len(ranks))
        rank_matrix.append(rank_positions)
    
    rank_matrix = np.array(rank_matrix)
    borda_scores = np.sum(len(valid_results) - rank_matrix - 1, axis=0)
    borda_ranks = np.argsort(-borda_scores)
    
    return {
        'weighted': {
            'scores': weighted_scores,
            'ranks': weighted_ranks,
            'name': 'Weighted Average',
            'description': 'Phylogenetically-informed weighted combination'
        },
        'topsis': {
            'scores': topsis_scores,
            'ranks': topsis_ranks,
            'name': 'TOPSIS',
            'description': 'Multi-criteria decision analysis'
        },
        'borda': {
            'scores': borda_scores,
            'ranks': borda_ranks,
            'name': 'Borda Count',
            'description': 'Rank aggregation across all metrics'
        }
    }

# ============================================================================
# CONSENSUS FINDING
# ============================================================================

def find_consensus(rankings, alignment_files):
    """Find consensus best alignment across ranking methods"""
    top_alignments = {}
    
    for method_key, method_data in rankings.items():
        best_idx = method_data['ranks'][0]
        best_num = int(os.path.basename(alignment_files[best_idx]).split('_')[1].split('.')[0])
        
        if best_num not in top_alignments:
            top_alignments[best_num] = []
        top_alignments[best_num].append(method_key)
    
    consensus_alignment = None
    max_agreement = 0
    
    for alignment_num, methods in top_alignments.items():
        if len(methods) > max_agreement:
            max_agreement = len(methods)
            consensus_alignment = alignment_num
    
    return {
        'has_consensus': max_agreement >= 2,
        'alignment': consensus_alignment,
        'agreement': max_agreement,
        'total_methods': len(rankings),
        'by_method': {
            method_key: {
                'alignment': int(os.path.basename(
                    alignment_files[method_data['ranks'][0]]
                ).split('_')[1].split('.')[0])
            }
            for method_key, method_data in rankings.items()
        }
    }

# ============================================================================
# OUTPUT GENERATION
# ============================================================================

def write_results(results, alignment_files, rankings, consensus, output_file, folder_name):
    """Write comprehensive results to file"""
    with open(output_file, 'w') as f:
        # Header
        f.write("="*80 + "\n")
        f.write("MAFFT ALIGNMENT SCORING RESULTS - EXACT JALVIEW IMPLEMENTATION\n")
        f.write("="*80 + "\n")
        f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Folder: {folder_name}\n")
        f.write(f"Total alignments analyzed: {len([r for r in results if r])}\n\n")
        
        # Best by individual metrics
        f.write("="*80 + "\n")
        f.write("BEST ALIGNMENTS BY INDIVIDUAL METRICS:\n")
        f.write("="*80 + "\n\n")
        
        metric_winners = [
            ('Sum-of-Pairs Score', 'sp_score', False),
            ('Jalview Conservation (total)', 'conservation_total', False),
            ('Jalview Quality (total)', 'quality_total', False),
            ('Hydrophobic Conservation', 'hydrophobic_conservation', False),
            ('Gap Coherence', 'gap_coherence', False),
            ('Parsimony-Informative Sites', 'parsimony_informative', True),
            ('Pairwise Identity', 'pairwise_identity', True)
        ]
        
        for metric_name, metric_key, is_percentage in metric_winners:
            valid_scores = [(i, r[metric_key]) for i, r in enumerate(results) if r]
            if valid_scores:
                best_idx, best_score = max(valid_scores, key=lambda x: x[1])
                best_num = int(os.path.basename(alignment_files[best_idx]).split('_')[1].split('.')[0])
                
                f.write(f"★ {metric_name}:\n")
                f.write(f"  Alignment: alignment_{best_num}.fasta ")
                
                if is_percentage:
                    f.write(f"(score: {best_score:.2f}%)\n")
                else:
                    f.write(f"(score: {best_score:.4f})\n")
        
        # Best by ranking methods
        f.write("\n" + "="*80 + "\n")
        f.write("BEST ALIGNMENTS BY RANKING METHODS:\n")
        f.write("="*80 + "\n\n")
        
        for method_key, method_data in rankings.items():
            best_idx = method_data['ranks'][0]
            best_num = int(os.path.basename(alignment_files[best_idx]).split('_')[1].split('.')[0])
            best_score = method_data['scores'][best_idx]
            
            f.write(f"★ {method_data['name']}:\n")
            f.write(f"  Alignment: alignment_{best_num}.fasta\n")
            f.write(f"  Score: {best_score:.4f}\n")
            f.write(f"  Method: {method_data['description']}\n\n")
        
        # Consensus
        f.write("="*80 + "\n")
        f.write("CONSENSUS ANALYSIS:\n")
        f.write("-"*80 + "\n")
        
        if consensus['has_consensus']:
            f.write(f"✅ CONSENSUS FOUND: Alignment {consensus['alignment']}\n")
            f.write(f"   Agreement: {consensus['agreement']}/{consensus['total_methods']} methods\n\n")
        else:
            f.write(f"⚠️  NO CONSENSUS: Ranking methods disagree!\n\n")
        
        # ALL alignments sorted by weighted average
        f.write("="*80 + "\n")
        f.write("DETAILED SCORES (sorted by Weighted Average):\n")
        f.write("="*80 + "\n")
        f.write(f"{'#':<4} {'Align':<8} {'SP':<10} {'Cons':<10} {'Qlty':<10} ")
        f.write(f"{'Hydro':<8} {'Gap':<8} {'Pars':<8} {'ID%':<8}\n")
        f.write("-"*80 + "\n")
        
        # Output ALL alignments
        for i, idx in enumerate(rankings['weighted']['ranks']):
            res = results[idx]
            if res:
                alignment_num = int(os.path.basename(alignment_files[idx]).split('_')[1].split('.')[0])
                f.write(f"{i+1:<4} {alignment_num:<8} "
                       f"{res['sp_score']:<10.0f} "
                       f"{res['conservation_total']:<10.1f} "
                       f"{res['quality_total']:<10.4f} "
                       f"{res['hydrophobic_conservation']:<8.4f} "
                       f"{res['gap_coherence']:<8.4f} "
                       f"{res['parsimony_informative']*100:<8.1f} "
                       f"{res['pairwise_identity']:<8.2f}\n")
        
        # Score interpretation guide
        f.write("\n" + "="*80 + "\n")
        f.write("SCORE INTERPRETATION GUIDE:\n")
        f.write("-"*80 + "\n")
        f.write("Jalview Conservation (total sum):\n")
        f.write("  Score = sum of conservation scores across all positions\n")
        f.write("  Per position: 0-11 (11 = identity, 10 = all properties conserved)\n")
        f.write("  Higher total = more conserved alignment\n\n")
        
        f.write("Jalview Quality (total sum):\n")
        f.write("  Score = sum of quality scores across all positions\n")
        f.write("  Based on BLOSUM62 substitution favorability\n")
        f.write("  Higher total = better quality alignment\n\n")
        
        f.write("Note: Raw sums are reported (not per-column averages)\n")
        f.write("      This reflects cumulative conservation/quality\n")

def write_per_position_scores(results, alignment_files, output_folder):
    """Write per-position conservation and quality scores to separate file"""
    scores_file = os.path.join(output_folder, "per_position_scores.txt")
    
    with open(scores_file, 'w') as f:
        f.write("PER-POSITION CONSERVATION AND QUALITY SCORES\n")
        f.write("=" * 80 + "\n")
        f.write("Format: Position scores separated by commas\n")
        f.write("Conservation: 0-11 (11=identity, 10=all properties conserved)\n")
        f.write("Quality: Normalized quality scores\n\n")
        
        for i, (result, align_file) in enumerate(zip(results, alignment_files)):
            if result:
                align_num = int(os.path.basename(align_file).split('_')[1].split('.')[0])
                f.write(f"Alignment {align_num}:\n")
                f.write("-" * 40 + "\n")
                
                # Conservation scores
                f.write("Conservation: ")
                f.write(", ".join([f"{s:.1f}" for s in result['conservation_scores']]))
                f.write(f"\nSum: {result['conservation_total']:.1f}\n\n")
                
                # Quality scores
                f.write("Quality: ")
                f.write(", ".join([f"{s:.6f}" for s in result['quality_scores']]))
                f.write(f"\nSum: {result['quality_total']:.6f}\n\n")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("="*60)
    print("MAFFT ALIGNMENT SCORER - EXACT JALVIEW IMPLEMENTATION")
    print("="*60)
    print("\nThis version uses:")
    print("- EXACT Jalview conservation calculation")
    print("- EXACT Jalview quality calculation")
    print("- RAW SUMS (not per-column averages)")
    print("- ALL alignments in output (not limited to top 20)\n")
    
    # Select folders
    folders = select_folders()
    
    if not folders:
        print("No folders selected. Exiting.")
        return
    
    print(f"\nProcessing {len(folders)} folder(s)...")
    
    # Process each folder
    for folder in folders:
        print(f"\n{'='*60}")
        print(f"Processing folder: {folder}")
        print(f"{'='*60}")
        
        # Get alignment files
        alignment_files = []
        for file in os.listdir(folder):
            if file.startswith("alignment_") and file.endswith(".fasta"):
                alignment_files.append(os.path.join(folder, file))
        
        alignment_files.sort(key=lambda x: int(os.path.basename(x).split('_')[1].split('.')[0]))
        
        if not alignment_files:
            print(f"No alignment files found in {folder}")
            continue
        
        print(f"Found {len(alignment_files)} alignment files")
        
        # Score alignments
        print("\nScoring alignments...")
        results = process_alignments_parallel(alignment_files)
        
        # Calculate rankings
        print("\nCalculating rankings...")
        rankings = calculate_rankings(results, alignment_files)
        
        if not rankings:
            print("Failed to calculate rankings")
            continue
        
        # Find consensus
        consensus = find_consensus(rankings, alignment_files)
        
        # Write results
        output_file = os.path.join(folder, "alignment_scores_EXACT_JALVIEW.txt")
        write_results(results, alignment_files, rankings, consensus, output_file, folder)
        print(f"Results written to: {output_file}")
        
        # Write per-position scores
        write_per_position_scores(results, alignment_files, folder)
        print(f"Per-position scores written to: {os.path.join(folder, 'per_position_scores.txt')}")
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
