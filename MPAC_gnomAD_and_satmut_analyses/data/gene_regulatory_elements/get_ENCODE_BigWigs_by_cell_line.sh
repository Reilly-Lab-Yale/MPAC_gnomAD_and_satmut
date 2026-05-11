#!/bin/sh

# Download DNase, H3K4me3, and H3K27ac BigWigs from ENCODE
# Cell lines: K562, HepG2, SK-N-SH

URL="https://www.encodeproject.org/files"

for acc in \
    ENCFF414OGC ENCFF806YEZ ENCFF849TDM \
    ENCFF546MZK ENCFF732PJK ENCFF795ONN \
    ENCFF280RMA ENCFF651WOM ENCFF262UEH
do
    echo "Downloading ${acc}..."
    wget -q "${URL}/${acc}/@@download/${acc}.bigWig"
done

echo "Done."
