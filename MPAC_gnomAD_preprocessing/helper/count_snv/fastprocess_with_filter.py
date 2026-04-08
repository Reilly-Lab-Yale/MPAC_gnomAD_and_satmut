# Like fastprocess.py but also tracks FILTER=PASS vs non-PASS
# Run with pypy3 and pipe file input

import sys

def main():
    for line in sys.stdin.buffer:
        line = line.decode()
        line = line.rstrip('\n').split('\t')

        info = line[-1].split(';')
        AF_str = [el for el in info if el.startswith("AF=")]
        if len(AF_str) != 1:
            continue

        AF_str = AF_str[0]
        AF = float(AF_str[len("AF="):])
        MAF = min(AF, 1 - AF)

        if MAF >= 0.01:
            varfreq = "common"
        else:
            varfreq = "rare"

        ref = line[3]
        alt = line[4]

        if alt == "." or ref == "." or len(alt) != 1 or len(ref) != 1:
            vartype = "indel"
        else:
            vartype = "snv"

        if len(alt) >= 10 or len(ref) >= 10:
            varlen = "long"
        else:
            varlen = "short"

        filt = line[6]
        if filt == "PASS":
            filtcat = "PASS"
        else:
            filtcat = "nonPASS"

        spit = [vartype, varfreq, varlen, filtcat]
        print('\t'.join(spit))

if __name__ == "__main__":
    main()
