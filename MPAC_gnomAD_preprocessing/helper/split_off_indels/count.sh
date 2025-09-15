cat slurm* | grep "^After filtering" | grep -v "Individual" > counts.txt
cat counts.txt | awk '{sum4 += $4; sum9 += $9} END {print "col4:", sum4, "col9:", sum9}' > summary.txt
