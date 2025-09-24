Some small scripts for querying gnomad v3.1.2

wrapper.sh produces output.txt

`grep -c '^      1 \[' output.txt`

reports 115776 gnomad variants do not have allele frequency designations. 

Definitions:
- common : MAF >=0.01. 
- rare : not common. (these are the same as the main pipeline, but repeated here for clairity).
- long : len(alt)>=10 or len(ref)>=10.
- short : not long.

Looking at the tail of output.txt we see

 794381 indel   common  long
3917119 indel   common  short
29627155 indel  rare    long
71557092 indel  rare    short
13410979 snv    common  short
608421662 snv   rare    short

(all SNVs are "short", of course.)

"SNVs are the plurality of human common variation"

(13410979 snv common) / (794381 indel common long + 3917119 indel common short + 13410979 snv common short)
 = 13410979 / (794381 + 3917119 + 13410979)
 = 13410979 / 18122479
 = 74%


"indels < 10bp in length are the vast majority of common indels"

(3917119 indel   common  short)/(3917119 indel   common  short +  794381 indel   common  long)
=3917119/(3917119+794381)
=3917119/4711500
=83%
