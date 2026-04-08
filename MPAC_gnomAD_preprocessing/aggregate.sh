echo "total"
awk '{sum += $2} END {print sum}' count.txt
echo "ref_long"
awk '{sum += $3} END {print sum}' count.txt
echo "alt_long"
awk '{sum += $4} END {print sum}' count.txt
echo "both_long"
awk '{sum += $5} END {print sum}' count.txt
