# Split stdin (e.g., from zcat) into N sequential output files, given total line count
# Usage:
#   zcat input.tsv.gz | pypy3 split_stream.py N total_lines [prefix] [outdir]
# Behavior:
#   - Reads text lines from stdin.
#   - Expects total line count (integer) as a second argument.
#   - Treats the first line as a header and writes it to all output files.
#   - Splits remaining lines sequentially into N files of roughly equal size.
#   - Single-pass, streaming; no memory buffering or temp files.

import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: pypy3 split_stream.py N total_lines [prefix] [outdir]", file=sys.stderr)
        sys.exit(1)

    try:
        n = int(sys.argv[1])
        total = int(sys.argv[2])
        if n <= 0 or total < 1:
            raise ValueError
    except ValueError:
        print("Error: N and total_lines must be positive integers.", file=sys.stderr)
        sys.exit(2)

    prefix = sys.argv[3] if len(sys.argv) >= 4 else "part_"
    outdir = sys.argv[4] if len(sys.argv) >= 5 else "."
    os.makedirs(outdir, exist_ok=True)

    base = (total - 1) // n  # subtract header line from total
    rem = (total - 1) % n
    sizes = [base + (1 if i < rem else 0) for i in range(n)]

    width = max(2, len(str(n - 1)))
    paths = [os.path.join(outdir, f"{prefix}{str(i).zfill(width)}.tsv") for i in range(n)]
    files = [open(p, "w") for p in paths]

    try:
        # Read and duplicate header line
        header = sys.stdin.readline()
        for f in files:
            f.write(header)

        # Stream the remaining lines to outputs
        i = 0
        remaining = sizes[i]
        for line in sys.stdin:
            files[i].write(line)
            remaining -= 1
            if remaining == 0 and i + 1 < n:
                i += 1
                remaining = sizes[i]
    finally:
        for f in files:
            try:
                f.flush()
                f.close()
            except Exception as e:
                print(f"Warning: error closing file {f.name}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
