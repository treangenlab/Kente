#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./kente_kraken2.sh IN_DIR OUT_ROOT Kraken2_DB [TOP_No OF GENERA TO FILTER TO] [THREADS]
#
# Example:
#   ./kente_kraken2.sh \
#     /path/your/metagenome_samples \
#     /path/to/your/output_folder \
#     /path/to/your/kraken2/database \
#     30 16

IN_DIR="$1"
OUT_ROOT="$2"
K2_DB="$3"
TOPN="${4:-30}"
THREADS="${5:-16}"

mkdir -p "$OUT_ROOT"

# finding all the contigs in your input directory with specific file extensions
mapfile -d '' FILES < <(
  find "$IN_DIR" -type f \( -name "*.fna" -o -name "*.fna.gz" \) -print0
)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No .fna or .fna.gz files found under: $IN_DIR" >&2
  exit 1
fi

echo "[INFO] Found ${#FILES[@]} contig files under $IN_DIR"
echo "[INFO] OUT_ROOT=$OUT_ROOT  TOPN=$TOPN  THREADS=$THREADS"

for CONTIGS in "${FILES[@]}"; do
  base="$(basename "$CONTIGS")"
  sample="$base"
  sample="${sample%.fna.gz}"
  sample="${sample%.fna}"

  OUTDIR="$OUT_ROOT/$sample"
  mkdir -p "$OUTDIR"

  KRAKEN_OUT="$OUTDIR/${sample}.kraken.out"
  KRAKEN_REPORT="$OUTDIR/${sample}.kraken.report.txt"
  GENERA_TXT="$OUTDIR/${sample}_genera.top${TOPN}.txt"

  echo "------------------------------------------------------------"
  echo "[INFO] Processing: $CONTIGS"
  echo "[INFO] Sample: $sample"
  echo "[INFO] Output dir: $OUTDIR"

  # Step1: Kraken2
  if [[ -s "$KRAKEN_OUT" && -s "$KRAKEN_REPORT" ]]; then
    echo "[SKIP] Kraken2 outputs already exist:"
    echo "       $KRAKEN_OUT"
    echo "       $KRAKEN_REPORT"
  else
    echo "[RUN ] Kraken2 -> $KRAKEN_OUT"
    if ! kraken2 --db "$K2_DB" --threads "$THREADS" --use-names --quick \
        --report "$KRAKEN_REPORT" \
        "$CONTIGS" > "$KRAKEN_OUT"; then
      echo "[WARN] Kraken2 failed for: $CONTIGS (skipping)" >&2
      continue
    fi
    echo "[OK  ] Kraken2 done"
  fi

  # Step2: Top-N genera
    
  if [[ -s "$GENERA_TXT" ]]; then
    echo "[SKIP] Genera list already exists: $GENERA_TXT"
  else
    echo "[RUN ] Extract top-${TOPN} genera -> $GENERA_TXT"

    # report columns:
    # 1:%  2:reads_clade  3:reads_taxon  4:rank  5:taxid  6..:name
    awk '$4=="G"{name=$6; for(i=7;i<=NF;i++) name=name" "$i; sub(/^[ \t]+/,"",name); print $2"\t"name}' \
      "$KRAKEN_REPORT" \
      | sort -nr -k1,1 \
      | awk -v top="$TOPN" 'NR<=top{print $2}' \
      > "$GENERA_TXT"

    echo "[OK  ] Wrote genera list: $GENERA_TXT"
  fi


done

echo "============================================================"
echo "[DONE] Batch preprocess complete."
