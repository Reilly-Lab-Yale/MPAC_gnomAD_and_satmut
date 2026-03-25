Some small scripts for querying gnomad v3.1.2

`wrapper.sh` produces `output.txt`

`grep -c '^      1 \[' output.txt` reports 115,776 gnomAD variants do not have allele frequency designations.

**Definitions:**
| Term    | Meaning |
|---------|---------|
| common  | MAF ≥ 0.01 |
| rare    | not common |
| long    | len(alt) ≥ 10 or len(ref) ≥ 10 |
| short   | not long |

(These match the main pipeline definitions.)

---

### Counts from `output.txt`

| Type  | Frequency | Length | Count         |
|-------|-----------|--------|---------------|
| indel | common    | long   | 794,381       |
| indel | common    | short  | 3,917,119     |
| indel | rare      | long   | 29,627,155    |
| indel | rare      | short  | 71,557,092    |
| snv   | common    | short  | 13,410,979    |
| snv   | rare      | short  | 608,421,662   |

All SNVs are "short" by definition.

---

### "SNVs are the plurality of human common variation"

|                          | Count      |
|--------------------------|------------|
| SNV, common              | 13,410,979 |
| indel, common, long      | 794,381    |
| indel, common, short     | 3,917,119  |
| **All common variants**  | **18,122,479** |

> 13,410,979 / 18,122,479 = **74%**

---

### "Short indels are the vast majority of common indels"

|                          | Count     |
|--------------------------|-----------|
| indel, common, short     | 3,917,119 |
| indel, common, long      | 794,381   |
| **All common indels**    | **4,711,500** |

> 3,917,119 / 4,711,500 = **83%**
