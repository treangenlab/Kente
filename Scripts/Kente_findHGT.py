#!/usr/bin/env python3
import argparse, csv
from collections import defaultdict, Counter, namedtuple

"""
Kente_findHGT.py
Goal:
  The goal here is to take the genus-block TSV from previous step (GAF->CIGAR->blocks) and call candidate HGT events.

Key definitions:
  - A = "host" genus for the contig (dominant genus by total labeled bp)
  - B = any non-host genus block that passes length threshold

Event types (strict-ish):
  - Sandwich (A-B-A): B is flanked by host blocks on both sides (strongest evidence)
  - Tip (A-B) or (B-A): B is at one end of contig with host on the other side
  - Mosaic: everything else (internal switches not forming strict sandwich/tip)

Important:
  UNLABELED blocks represent windows where no genus wins (no support, ties, etc).
  We allow skipping over small unlabeled gaps when searching for host flanks,
  but we cap the total skipped unlabeled bp with max_unlab_bp.
"""

Block = namedtuple("Block", "contig start end genus")

UNLAB = {"UNLABELED", "", None}

#Read the tsv file produced in previous step
def read_blocks(path):
    by_contig = defaultdict(list)
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        need = {"contig_id", "start_bp", "end_bp", "genus"}
        # Allow flexible naming if needed, but sticking to my previous input format
        missing = need - set(reader.fieldnames or [])
        if missing:
            # SHould be contig_id but incase the label is contig
            if "contig" in reader.fieldnames:
                pass # Acceptable variation
            else:
                raise SystemExit(f"Missing columns in {path}: {', '.join(sorted(missing))}")

        for row in reader:
            c = (row.get("contig_id") or row.get("contig") or "").strip()
            if not c:
                continue
            try:
                s = int(row.get("start_bp") or row.get("start"))
                e = int(row.get("end_bp") or row.get("end"))
                g = (row.get("genus") or "").strip()
            except:
                continue
            if e <= s:
                continue
            by_contig[c].append(Block(c, s, e, g))

    # Here, we ensure the blocks are in order for neighbour searching
    for c in by_contig:
        by_contig[c].sort(key=lambda b: b.start)
    return by_contig

def contig_length(blocks):
    return max(b.end for b in blocks) if blocks else 0

def infer_host(blocks):
    """
    Infers the host (Clade_A) for THIS SPECIFIC CONTIG.
    """
    counts = Counter()
    for b in blocks:
        if b.genus in UNLAB:
            continue
        counts[b.genus] += (b.end - b.start)
    return counts.most_common(1)[0][0] if counts else None

"""
Walk left/right from index i to find the nearest *labeled* neighbor block.
We allow skipping over small UNLABELED gaps
However, we cap the *total* amount of unlabeled bp we can skip on that side.
"""
def skip_unlabeled(blocks, i, direction, max_unlab_bp):
    j = i + direction
    skipped = 0
    while 0 <= j < len(blocks):
        if blocks[j].genus in UNLAB:
            skipped += (blocks[j].end - blocks[j].start)
            if skipped > max_unlab_bp:
                return None, skipped
            j += direction
            continue
        return j, skipped
    return None, skipped

def flank_bp(blocks, host, i, direction, max_unlab_bp):
    j, skipped = skip_unlabeled(blocks, i, direction, max_unlab_bp)
    if j is None:
        return 0, skipped, False
    if blocks[j].genus != host:
        return 0, skipped, False
    return (blocks[j].end - blocks[j].start), skipped, True

def classify_event(blocks, host, i, max_unlab_bp, min_flank_bp):
    left_bp, left_unlab, left_is_host = flank_bp(blocks, host, i, -1, max_unlab_bp)
    right_bp, right_unlab, right_is_host = flank_bp(blocks, host, i, +1, max_unlab_bp)

    at_left_edge = (skip_unlabeled(blocks, i, -1, max_unlab_bp)[0] is None)
    at_right_edge = (skip_unlabeled(blocks, i, +1, max_unlab_bp)[0] is None)

    if left_is_host and right_is_host and left_bp >= min_flank_bp and right_bp >= min_flank_bp:
        return "Sandwich (A-B-A)"
    if left_is_host and at_right_edge:
        return "Tip (A-B)"
    if right_is_host and at_left_edge:
        return "Tip (B-A)"
    if at_left_edge or at_right_edge:
        return "Fragmented"
    return "Complex/Mosaic"

def find_hgt(by_contig, min_insert_bp, max_unlab_bp, min_flank_bp):
    events = []

    for contig, blocks in by_contig.items():
        # 1. Infer Host (Clade_A) for this contig only
        host = infer_host(blocks)
        if not host:
            continue

        c_len = contig_length(blocks)

        for i, curr in enumerate(blocks):
            if curr.genus in UNLAB:
                continue
            if curr.genus == host:
                continue

            length = curr.end - curr.start
            if length < min_insert_bp:
                continue

            etype = classify_event(blocks, host, i, max_unlab_bp, min_flank_bp)

            conf = (length / c_len) if c_len > 0 else 0.0

            events.append({
                "contig_ID": contig,
                "contig_len": c_len,
                "Clade_A": host,           
                "Clade_B": curr.genus,   
                "start": curr.start,
                "end": curr.end,
                "length": length,
                "event_type": etype,
            
            })

    return events

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Blocks TSV from classify_hybrid")
    
    # Splitting the outputs high and low confidence HGT calls. I trust the A-B-A calls more. 
    ap.add_argument("--out_aba", required=True, help="Output for SANDWICH events (High Confidence)")
    ap.add_argument("--out_other", required=True, help="Output for Other events (Tips, Mosaic)")
    
    ap.add_argument("--min_insert", type=int, default=500, help="Min BP for foreign block")
    ap.add_argument("--max_unlab_bp", type=int, default=2000, help="Bridge unlabeled gaps")
    ap.add_argument("--min_flank_bp", type=int, default=500, help="Min host flank bp")
    args = ap.parse_args()

    data = read_blocks(args.input)
    candidates = find_hgt(data, args.min_insert, args.max_unlab_bp, args.min_flank_bp)

    print(f"Scanned {len(data)} contigs.")
    print(f"Total HGT Candidates Found: {len(candidates)}")

    # Split the results
    sandwich_evs = [ev for ev in candidates if "Sandwich" in ev["event_type"]]
    other_evs = [ev for ev in candidates if "Sandwich" not in ev["event_type"]]

    print(f"  -> Sandwich (A-B-A): {len(sandwich_evs)} written to {args.out_aba}")
    print(f"  -> Other (Tips/Mosaic): {len(other_evs)} written to {args.out_other}")

    cols = ["contig_ID", "contig_len", "Clade_A", "Clade_B", "start", "end", "length", "event_type"]

    # Write Sandwich events
    with open(args.out_aba, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(sandwich_evs)

    # Write Other possible events
    with open(args.out_other, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(other_evs)

if __name__ == "__main__":
    main()