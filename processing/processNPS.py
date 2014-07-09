#!/usr/bin/python
import sys
from math import floor

def writeNPSBed(npsfn, npsbedfn):
    '''Write the NPS file in proper bed format.  Add the center as thick'''
    try:
        npsfh = open(npsfn)
    except:
        print "Error opening NPS file:  ", npsfn
        return
    try:
        npsbedfh = open(npsbedfn, 'w')
    except:
        print "Error opening NPS bed (outfile):  ", npsbedfn

    # burn the header
    npsfh.readline()
    
    print "Writing full .bed file for NPS data to %s..." % (npsbedfn)
    for line in npsfh:
        line = line.strip().split('\t')
        chrom = line[0]
        start = int(line[1])
        stop = int(line[2])
        name = line[3]
        score = line[4]
        strand = '.'
        thickStart = int(start + floor((stop - start) / 2))
        thickEnd = thickStart + 1
        newline = '\t'.join([chrom, str(start), str(stop), name, score, strand, str(thickStart), str(thickEnd)])
        npsbedfh.write(newline + '\n')
        
    npsfh.close()
    npsbedfh.close()



if __name__ == '__main__':

    npsfn = sys.argv[1]
    npsbedfn = sys.argv[2]


    writeNPSBed(npsfn, npsbedfn)
    
