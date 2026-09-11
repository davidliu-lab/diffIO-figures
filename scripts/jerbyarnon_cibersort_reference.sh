#!/usr/bin/env bash
# Build the CIBERSORTx "Create Signature Matrix" reference file from the
# Jerby-Arnon single-cell data.
#
# Usage (from the project root):
#   bash scripts/jerbyarnon_cibersort_reference.sh
#
# Takes the two-file export from jerbyarnon_cibersort_matrix_export.R and
# produces the single file CIBERSORTx actually wants, applying three fixes:
#
#   1. Header. CIBERSORTx expects the header row to be cell-type labels, one
#      per cell, not barcodes with a separate label map.
#   2. Log space. so@assays$RNA@data is Seurat LogNormalize output -- natural
#      log, CP10K. Verified: sum(exp(x)-1) over a cell == 10000.0 exactly.
#      CIBERSORTx requires non-log linear space, so we apply expm1.
#   3. Labels. Renamed to LM22 conventions, because
#      bulk-analysis-other-immune-signatures.qmd reads the result using LM22
#      column names ("T cells CD8" etc).
#
# Cells labelled "?" are dropped.
#
# Streams with awk rather than loading in R -- the input is ~834 MB.

set -euo pipefail
cd "$(dirname "$0")/.."

EXPR=data/scRNAseq/so_jerbyarnon_expression.txt
META=data/scRNAseq/so_jerbyarnon_metadata.txt
OUTDIR=results/my_cache/bulk-analysis-other-immune-signatures
OUT="$OUTDIR/jerbyarnon_cibersort_reference.txt"

for f in "$EXPR" "$META"; do
  [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

awk -F'\t' -v OFS='\t' '
  # LM22-style names for the populations the chapter sums.
  BEGIN {
    rename["T.CD8"]      = "T cells CD8"
    rename["T.CD4"]      = "T cells CD4"
    rename["B.cell"]     = "B cells"
    rename["Macrophage"] = "Macrophages"
    rename["NK"]         = "NK cells"
    rename["T.cell"]     = "T cells other"
    rename["Mal"]        = "Malignant"
    rename["Endo."]      = "Endothelial"
    rename["CAF"]        = "CAF"
  }

  # Pass 1: barcode -> cell type.
  NR == FNR {
    if (FNR > 1) type[$1] = $2
    next
  }

  # Pass 2, header: decide which columns survive, and emit their labels.
  FNR == 1 {
    printf "GeneSymbol"
    for (i = 2; i <= NF; i++) {
      t = type[$i]
      if (t == "" || t == "?") continue          # unlabelled or ambiguous
      keep[++n] = i
      printf "%s%s", OFS, (t in rename ? rename[t] : t)
      kept[(t in rename ? rename[t] : t)]++
    }
    printf "\n"
    for (t in kept) print "  " t ": " kept[t] " cells" > "/dev/stderr"
    print "  kept " n " of " (NF-1) " cells" > "/dev/stderr"
    next
  }

  # Pass 2, data: expm1 back to linear CP10K.
  {
    printf "%s", $1
    for (j = 1; j <= n; j++) {
      x = $(keep[j]) + 0
      printf "%s%s", OFS, (x == 0 ? "0" : sprintf("%.6g", exp(x) - 1))
    }
    printf "\n"
  }
' "$META" "$EXPR" > "$OUT"

echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
