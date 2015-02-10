#!/usr/bin/env python

import sys, os
from analyze_common import *
use_message = '''
'''

def intron_distribution():
    chr_dic = None
    trans_dic, trans_ids = extract_transcripts("genes.gtf", verbose = False)
    junctions_dic = extract_junctions(trans_dic)

    max, interval = 1 << 20, 1 << 10
    assert interval < max
    counts = [0 for i in range(max / interval)]
    total = 0
    for key, value in junctions_dic.items():
        chr, left, right = key.split("-")
        left, right = int(left), int(right)
        size = right - left + 1

        if size >= max:
            size = max - 1
            
        counts[size / interval] += 1
        total += 1

    so_far = 0
    for i in range(len(counts)):
        so_far += counts[i]

        # print "%d - %d : %.4f" % (i * interval, (i+1) * interval - 1, float(so_far) / total * 100.0)
        print "%d\t%d" % ((i+1) * interval, counts[i])

        
    
if __name__ == "__main__":
    intron_distribution()
