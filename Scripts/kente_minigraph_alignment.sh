#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Kente minigraph mapping (batch, single script)
# For each contig file in CONTIGS_DIR:
#   - reads OUT_ROOT/<sample>/<sample>_genera.topN.txt
#   - maps contigs to each genus graph in GRAPHS_ROOT
#   - writes GAFs to OUT_ROOT/<sample>/gaf_by_genus/<genus>.gaf
#   - runs multiple genera in parallel (JOBS), each minigraph uses THREADS_PER_JOB threads
#   - skips output GAF if it already exists (This is incase the batch stops, resume-safe)
#
# Usage:
#   ./kente_minigraph_alignment.sh CONTIGS_DIR OUT_ROOT GRAPHS_ROOT [TOPN] [JOBS] [THREADS_PER_JOB]
#
# Example:
#   ./kente_minigraph_alignment.sh \
#     /path/your/metagenome_samples \
#     /path/to/your/output_folder \
#    /path/to/pangenome_graphs_db \
#     30 8 4
# ------------------------------------------------------------

CONTIGS_DIR="$1"
OUT_ROOT="$2"
GRAPHS_ROOT="$3"
TOPN="${4:-30}"
JOBS="${5:-8}"
THREADS_PER_JOB="${6:-4}"

shopt -s nullglob
FILES=( "$CONTIGS_DIR"/*.fna )
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No .fna files found in: $CONTIGS_DIR" >&2
  exit 1
fi

echo "[INFO] Found ${#FILES[@]} contig files"
echo "[INFO] OUT_ROOT=$OUT_ROOT"
echo "[INFO] GRAPHS_ROOT=$GRAPHS_ROOT"
echo "[INFO] TOPN=$TOPN  JOBS=$JOBS  THREADS_PER_JOB=$THREADS_PER_JOB"
echo

for CONTIGS in "${FILES[@]}"; do
  base="$(basename "$CONTIGS")"
  sample="${base%.fna}"

  OUTDIR="$OUT_ROOT/$sample"
  GENERA_TXT="$OUTDIR/${sample}_genera.top${TOPN}.txt"
  GAF_DIR="$OUTDIR/gaf_by_genus"

  echo "============================================================"
  echo "[INFO] Sample:  $sample"
  echo "[INFO] Contigs: $CONTIGS"
  echo "[INFO] Outdir:  $OUTDIR"

  if [[ ! -s "$GENERA_TXT" ]]; then
    echo "[WARN] Missing genera list (run preprocess first): $GENERA_TXT" >&2
    continue
  fi

  mkdir -p "$GAF_DIR"

  export CONTIGS GRAPHS_ROOT GAF_DIR THREADS_PER_JOB

  run_one_genus() {
    genus="$1"
    [[ -z "$genus" ]] && exit 0

    gfa="${GRAPHS_ROOT}/${genus}/${genus}.rgfa"
    if [[ ! -f "$gfa" ]]; then
      echo "[WARN] missing graph: $gfa" >&2
      exit 0
    fi

    out="${GAF_DIR}/${genus}.gaf"
    if [[ -s "$out" ]]; then
      echo "[SKIP] $sample $genus (exists)"
      exit 0
    fi

    echo "[RUN ] $sample $genus"
    minigraph -cx lr -t "$THREADS_PER_JOB" "$gfa" "$CONTIGS" > "$out"
    echo "[OK  ] $sample $genus"
  }
  export -f run_one_genus

  echo "[INFO] Mapping genera from: $GENERA_TXT"
  echo "[INFO] Writing GAFs to:     $GAF_DIR"
  echo "[INFO] Parallel: JOBS=$JOBS  THREADS_PER_JOB=$THREADS_PER_JOB"
  echo

  # Run genera in parallel
  grep -v '^[[:space:]]*$' "$GENERA_TXT" \
    | xargs -I{} -P "$JOBS" bash -c 'run_one_genus "$@"' _ {}

  echo "[DONE] Finished sample: $sample"
done

echo "============================================================"
echo "[DONE] All samples complete."
