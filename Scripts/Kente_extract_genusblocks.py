#!/usr/bin/env python3
import sys, os, glob, argparse, re
import numpy as np
from collections import defaultdict

"""
Kente_extract_genusblocks.py

Goal:
  The goal here is to convert graph aligned per-genus GAF files into a single TSV of "genus blocks" per contig.

Pipeline logic:
  1) Read all .gaf files (each file corresponds to mapping contigs to ONE genus graph).
  2) For each alignment line, parse the CIGAR (cg:Z:...) to mark which QUERY bases are truly covered.
     - '='/'X'/'M consume query -> counts as coverage
     - 'I' consumes query -> NO coverage for the graph (gap on query)
     - 'D' consumes reference only -> does not move query pointer
  3) Weight coverage by mapping quality (MAPQ) and identity (1 - dv).
  4) In fixed windows (e.g., 500bp), pick the genus with the highest total support.
  5) Merge adjacent windows with the same winner -> blocks.
"""
#This function infers the genus name from gaf file generated in the previous step
def clean_genus_from_filename(gaf_path):
    base = os.path.basename(gaf_path).rsplit(".", 1)[0]
    if "vs_" in base:
        return base.split("vs_")[-1].split("_")[0]
    return base.split("_")[-1]

def parse_dv(cols):
    for c in cols[12:]:
        if c.startswith("dv:f:"):
            try: return float(c.split(":")[-1])
            except: return 0.1
    return 0.1

def parse_cigar_coverage(cigar_str, start_bp, mapq, identity_factor, scoreboard, g_idx):
    """
    Parses minigraph CIGAR (cg:Z:...) and updates the scoreboard array.
    Operations:
      M, =, X: Alignment match/mismatch -> ADDS SCORE
      I: Insertion in Query (Query has it, Graph doesn't) -> NO SCORE, ADVANCE IDX
      D: Deletion in Query (Graph has it, Query doesn't) -> NO SCORE, NO ADVANCE
    """
    # HEre, we are finding patterns in our cigar string and giving scores
    # Minigraph uses =, X, I, D. Sometimes M, I think.
    ops = re.findall(r'(\d+)([=XIDM])', cigar_str)
    
  # query pointer begins at query-start from GAF 
    curr_pos = start_bp
    score_val = (mapq + 1) * identity_factor
    
    for length_str, op in ops:
        length = int(length_str)
        
        if op == '=' or op == 'X' or op == 'M':
            # These bases are covered by the graph
            # We place a cap at end_pos to avoid indexing a past contig lenght 
            end_pos = min(curr_pos + length, len(scoreboard))
            if end_pos > curr_pos:
                scoreboard[curr_pos:end_pos, g_idx] += score_val
            curr_pos += length
            
        elif op == 'I':
            # Insertion in Query: The query has sequence here that the Graph DOES NOT.
            
            curr_pos += length
            
        elif op == 'D':
            # Deletion in Query: The Graph has sequence, Query does not.
            # Pointer does not move on Query.
            pass

def parse_gaf_files(gaf_dir, min_mapq):
    # We need to know contig lengths first to build arrays
    contig_lengths = {}
    files = glob.glob(os.path.join(gaf_dir, "*.gaf"))
    
    # Pass 1: Get Lengths from the GAF file
    print(f"[INFO] Scanning {len(files)} files ...", file=sys.stderr)
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"): continue
                cols = line.split("\t")
                try: contig_lengths[cols[0]] = int(cols[1])
                except: continue

    # Initialize Scoreboards: [Length x Genera]
    # We will identify genera indices on the fly
    genus_set = set()
    for f in files:
        genus_set.add(clean_genus_from_filename(f))
    genera = sorted(list(genus_set))
    gen2idx = {g:i for i,g in enumerate(genera)}
    
    # Dict of numpy arrays: { contig: np.zeros((len, n_genera)) }
    # Allocate scoreboard per contig:
    #   rows = contig length (per-base)
    #   cols = number of genera we mapped against
    # float32 saves memory
    
    scoreboards = { c: np.zeros((l, len(genera)), dtype=np.float32) for c,l in contig_lengths.items() }
    
    # Pass 2: Filling the Scoreboards using the CIGAR
    print(f"[INFO] Parsing GAF for exact coverage...", file=sys.stderr)
    for f in files:
        genus = clean_genus_from_filename(f)
        g_idx = gen2idx[genus]
        
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                cols = line.split("\t")
              #Core fields needed from the GAF file. The GAF is very important, lol!   
                try:
                    qname = cols[0] #Query name (contig ID)
                    # qlen = int(cols[1])
                    qs = int(cols[2])
                    # qs (Query start)
                    mapq = int(cols[11])
                except: continue

                #Let us filter noisy alignments to get rid of false positives down the line
                if mapq < min_mapq: continue
                
                # Extract CIGAR
                cigar = None
                dv = 0.1
                for tag in cols[12:]:
                    if tag.startswith("cg:Z:"):
                        cigar = tag[5:]
                    elif tag.startswith("dv:f:"):
                        dv = float(tag[5:])
                
                if not cigar: continue
                
                identity_factor = max(0.1, 1.0 - dv)
                
                if qname in scoreboards:
                    parse_cigar_coverage(cigar, qs, mapq, identity_factor, scoreboards[qname], g_idx)

    return scoreboards, genera
#Winners here are the best fit for the contigs and which genus they aligned to the most. Our default window size is 500bp. We can always change it
def get_winners(scoreboards, genera, win_size=500):
    blocks = []
    
    for contig, board in scoreboards.items():
        c_len = board.shape[0]
        n_wins = (c_len + win_size - 1) // win_size
        
        # Calculate winner for each window based on sum of scores
        win_labels = []
        for w in range(n_wins):
            s = w * win_size
            e = min((w+1) * win_size, c_len)
            
            # Sum scores in this window for each genus
            slice_sums = np.sum(board[s:e, :], axis=0)
            max_val = np.max(slice_sums)

            # We deal with UNLABeled regions: No genus has any support in this window if your score is 0. 
            if max_val == 0:
                win_labels.append("UNLABELED")
            else:
                winners = np.where(slice_sums == max_val)[0]
                if len(winners) == 1:
                    win_labels.append(genera[winners[0]])
                else:
                    win_labels.append("UNLABELED") # Strict tie
        
        # Merge Labels
        if not win_labels: continue
        curr = win_labels[0]
        start_idx = 0
        for i in range(1, len(win_labels)):
            if win_labels[i] != curr:
                blocks.append((contig, start_idx*win_size, i*win_size, curr))
                curr = win_labels[i]
                start_idx = i
        blocks.append((contig, start_idx*win_size, c_len, curr))
        
    return blocks

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaf_dir", required=True, help="Directory containing per-genus .gaf files")
    ap.add_argument("--out", required=True, help="Output TSV of genus blocks")
    ap.add_argument("--win", type=int, default=500, help="Window size in bp (default: 500)")
    ap.add_argument("--min_mapq", type=int, default=1, help="Minimum MAPQ to keep an alignment (default: 1)")
    args = ap.parse_args()
    
    scoreboards, genera = parse_gaf_files(args.gaf_dir, args.min_mapq)
    
    print(f"[INFO] Computing best alignment matches along the graph...", file=sys.stderr)
    blocks = get_winners(scoreboards, genera, args.win)
    
    with open(args.out, "w") as out:
        out.write("contig_id\tstart_bp\tend_bp\tgenus\n")
        for b in blocks:
            out.write(f"{b[0]}\t{b[1]}\t{b[2]}\t{b[3]}\n")

if __name__ == "__main__":
    main()
