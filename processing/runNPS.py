#!/usr/bin/python
import sys
import subprocess
from math import floor

def generateNPSCalls(sample_metadata_dict):
    '''
    '''
    runNPS(sample_metadata_dict)
    writeNPSBed(sample_metadata_dict)
    

def runNPS(sample_metadata_dict):
    '''
    '''
    nps_bin = 'bin/SeqTag.py'
    nps_par_fn = 'processing/NPS.par'

    print "Calculating NPS peaks for sample %s" % (sample_metadata_dict['sample_name'])

    bed_fn = sample_metadata_dict['bam_fn'].replace('.bam', '.bed')
    nps_out_fn = bed_fn.replace('.bed', '.nps')

    nps_log_fn = bed_fn.replace('.bed', '.nps.log')
    nps_log_fh = open(nps_log_fn, 'w')

    cmd_args = [nps_bin, nps_par_fn, bed_fn, nps_out_fn]

    proc = subprocess.Popen(cmd_args, stderr=nps_log_fh, stdout=nps_log_fh)
    proc.wait()

    sample_metadata_dict['nps_out_fn'] = nps_out_fn

def writeNPSBed(sample_metadata_dict):
    '''
    Write the NPS file in proper bed format.  Add the center as thick
    '''
    nps_fn = sample_metadata_dict['nps_out_fn']
    nps_bed_fn = nps_fn.replace('.nps', '.processed.nps')

    nps_fh = open(nps_fn)
    nps_bed_fh = open(nps_bed_fn, 'w')

    # burn the header
    nps_fh.readline()

    print "Writing full BED-format file for NPS data to %s..." % (nps_bed_fn)
    for line in nps_fh:
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
        nps_bed_fh.write(newline + '\n')

    nps_fh.close()
    nps_bed_fh.close()

    sample_metadata_dict['nps_bed_fn'] = nps_bed_fn

