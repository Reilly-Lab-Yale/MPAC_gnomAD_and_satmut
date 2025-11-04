# Create temporary directory for test
OUTDIR=$(mktemp -d)
echo "Writing outputs to: $OUTDIR"

# Make a small gzipped TSV example
cat <<'EOF' | gzip > "$OUTDIR/input.tsv.gz"
col1	col2	col3
a1	b1	c1
a2	b2	c2
a3	b3	c3
a4	b4	c4
a5	b5	c5
EOF

# Count total lines (including header)
TOTAL=$(zcat "$OUTDIR/input.tsv.gz" | wc -l)
N=2

echo "Splitting $TOTAL total lines into $N parts..."

module load miniconda
conda activate speedracer

# Run the splitter
zcat "$OUTDIR/input.tsv.gz" | pypy3 break.py "$N" "$TOTAL" part_ "$OUTDIR"

# Show outputs
echo "=== Output files ==="
for f in "$OUTDIR"/part_*.tsv; do
    echo "--- $f ---"
    cat "$f"
    echo
done

# Cleanup note
echo "Test complete. Files are in $OUTDIR"