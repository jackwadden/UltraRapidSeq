#!/usr/bin/python

import sys
import os

fn = "pileup.txt"
coverage_thresh = 5

if not os.path.isfile(fn):
    print("File not found...")
    sys.exit()


with open(fn) as fp:

    hotspot_count = 0
    hotspot_read_count = 0
    in_hotspot = False
    max_coverage = 0
    hotspot_chr = ""
    hotspot_start = 0
    hotspot_end = 0
    
    for line in fp:

        # parse each line (chr, loc, base, coverage, codes, quality)
        fields = line.split()

        coverage = int(fields[3])

        if not in_hotspot :
            if coverage > coverage_thresh :
                in_hotspot = True
                hotspot_chr = str(fields[0])
                hotspot_start = str(fields[1])
                hotspot_count = hotspot_count + 1
                max_coverage = coverage
                #print(line)
                
        else :
            #print(line)

            if coverage > max_coverage:
                max_coverage = coverage                
            
            if coverage < coverage_thresh and in_hotspot:
                #print(max_coverage)
                hotspot_read_count = hotspot_read_count + max_coverage
                hotspot_end = str(fields[1])
                print(hotspot_chr + "\t" + hotspot_start + "\t" + hotspot_end)
                max_coverage = 0
                in_hotspot = False


    #print("Num hotspots: ", hotspot_count)
    #print("Hotspot reads: ", hotspot_read_count)

            
