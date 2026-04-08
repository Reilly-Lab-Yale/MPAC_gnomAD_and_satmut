#!/usr/bin/env bash
set -euo pipefail

# Usage: ./count_ref_alt_lengths.sh /path/to/top/dir
ROOT_DIR="${1:-.}"

if ! command -v gzip >/dev/null 2>&1; then
  echo "Error: gzip not found in PATH." >&2
  exit 1
fi

grand_total=0
grand_ref_long=0
grand_alt_long=0
grand_both_long=0

# Find all *.csv.gz files under ROOT_DIR (recursively)
while IFS= read -r -d '' file; do
  # total = all rows (excluding header)
  # ref_long = rows with length(col4) > 10
  # alt_long = rows with length(col5) > 10
  # both_long = rows with both col4>10 AND col5>10
  read -r total ref_long alt_long both_long < <(
    gzip -cd -- "$file" | awk -F'\t' '
      NR > 1 {
        total++
        if (length($4) > 10) ref_long++
        if (length($5) > 10) alt_long++
        if (length($4) > 10 && length($5) > 10) both_long++
      }
      END {
        print (total+0), (ref_long+0), (alt_long+0), (both_long+0)
      }
    '
  )

  printf '%s\t%d\t%d\t%d\t%d\n' \
    "$file" "$total" "$ref_long" "$alt_long" "$both_long"

  grand_total=$((grand_total + total))
  grand_ref_long=$((grand_ref_long + ref_long))
  grand_alt_long=$((grand_alt_long + alt_long))
  grand_both_long=$((grand_both_long + both_long))
done < <(find "$ROOT_DIR" -type f -name '*.csv.gz' -print0)

# Overall summary across all files
printf 'TOTAL\t%d\t%d\t%d\t%d\n' \
  "$grand_total" "$grand_ref_long" "$grand_alt_long" "$grand_both_long"
