#!/usr/bin/python
import sys
sys.path.append('./lib/python2.7/site-packages/')
import os
import subprocess
from tools import *

def filterSAMByMapq(sam_fn, min_mapq_score=37):
    '''
    Removes all reads of quality < min_mapq_score from the SAM file
    '''
    # removes reads with less than the specified mapq_score
    print "Removing reads with a mapq_score less than %s..." % (min_mapq_score)
    tmp_sam_fn = sam_fn.replace('.sam', '.tmp.sam')
    tmp_sam_fh = open(tmp_sam_fn, 'w')
    with open(sam_fn) as sam_fh:
        good_count = 0
        total_count = 0
        for line in sam_fh:
            if line.startswith('@'):
                tmp_sam_fh.write(line)
                continue
            total_count += 1
            mapq_score = int(line.split('\t')[4])
            if mapq_score >= min_mapq_score:
                good_count += 1
                tmp_sam_fh.write(line)
        print "...retained %s of %s reads (%s)." % (good_count, total_count, 100.*good_count/total_count)
    tmp_sam_fh.close()

    os.rename(tmp_sam_fn, sam_fn)

def convertSAMtoBAM(sam_fn):
    bam_fn = sam_fn.replace('.sam', '.bam')
    print 'Converting %s to .bam format...' % (sam_fn)
    cmd_args = ['bin/samtools', 'view', '-bS', sam_fn, '-o', bam_fn]
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    stdout, stderr = proc.communicate()
    print stdout, stderr
    print '\t...done.'

    return bam_fn

def sortBED(bed_fn):
    '''
    use GNU sort to sort the bed file
    '''
    print "Sorting the bed file %s with GNU sort..." % (bed_fn)
    cmd_args = ['sort', '-k1,1', '-k2,2n', bed_fn, '-o', bed_fn]
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sout, serr = proc.communicate()
    if serr:
        print serr

def sortBAM(bam_fn):
    print 'Sorting %s...' % (bam_fn)
    # this split below is required b/c samtools sort wants to always add .bam to the file
    cmd_args = ['bin/samtools', 'sort', bam_fn, bam_fn.split('.bam')[0]]
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    stdout, stderr = proc.communicate()
    print stdout, stderr
    print '\t...done.'

def makeBAI(bam_fn):
    bai_fn = bam_fn.replace('.bam', '.bai')
    print 'Making .bai file for %s...' % (bam_fn)
    cmd_args = ['bin/samtools', 'index', bam_fn, bai_fn]
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    stdout, stderr = proc.communicate()
    print stdout, stderr
    print '\t...done.'

def removeDuplicateReads(bam_fn):
    print "Removing exact mapped duplicates..."
    rmdup_fn = bam_fn.replace('.bam', '.rmdup')
    cmd_args = ['bin/samtools', 'rmdup', '-s', bam_fn, rmdup_fn]
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    stdout, stderr = proc.communicate()
    print stdout, stderr
    print '\t...done.'

    os.rename(bam_fn, bam_fn.replace('.bam', '.withduplicates.bam') )
    os.rename(rmdup_fn, bam_fn)

def makeBAMfromSAM(sam_fn, min_mapq_score=37, remove_duplicate_reads_flag=False):
    '''
    Makes a BAMfile from a SAM file with the follow steps:
    1) filter on the given min_mapq_score
    2) Convert to BAM format
    3) Sort the BAM file
    4) Make a BAI file
    5) Remove exact duplicate reads
    6) Remake the BAI file
    '''
    # how samtools sort works...

    filterSAMByMapq(sam_fn, min_mapq_score)
    bam_fn = convertSAMtoBAM(sam_fn)
    sortBAM(bam_fn)
    makeBAI(bam_fn)
    if remove_duplicate_reads_flag:
        removeDuplicateReads(bam_fn)
    makeBAI(bam_fn)

    return bam_fn

def calculateMappedLibrarySize(bam_fn):
    '''
    Uses SAM Tools to calculate the mapped library size of the given
    BAM file.
    '''
    print "Calculating the mapped library size in the BAM file %s" % (bam_fn)

    cmd_args = ['bin/samtools', 'idxstats', bam_fn]
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # releveant output is from stderr
    stderr, stdout = proc.communicate()
    size = 0
    for line in stderr.split('\n'):
        line = line.split()
        if line:
            size += int(line[2])

    return size

def makeWigFromBED(bed_fn, chromosome_size_fn):
    '''
    make the density coverage (WIG) file

    returns the wig_fn
    '''
    print "Making the density coverage file with genomeCoverageBed..."

    wig_fn = bed_fn.replace('.bed', '.wig')
    wig_fh = open(wig_fn, 'w')

    cmd_args = ['bin/genomeCoverageBed', '-i', bed_fn, '-g', chromosome_size_fn, '-dz']
    #cmd_args = ['bin/genomeCoverageBed', '-i', bed_fn, '-g', chromosome_size_fn, '-d']
    proc = subprocess.Popen(cmd_args, stdout=wig_fh, stderr=wig_fh)
    proc.wait()
    #densities = proc.communicate()[0].split('\n')
    wig_fh.close()

    return wig_fn

def scaleWigFile(wig_fn, scaled_library_size, mapped_library_size):
    '''
    scales density coverage (WIG) file to the scaled library size.

    returns the name of the scaled WIG file.
    '''
    print "Writing scaled density file..."
    wig_fh = open(wig_fn, 'r')
    scaled_wig_fn = wig_fn.replace('.wig', '.scaled.wig')
    scaled_wig_fh = open(scaled_wig_fn, 'w')
    chrom = ''
    #lastline = ''

    for line in wig_fh:
        line = line.split('\t')
        try:
            current_chrom = line[0]
            # wig files are 1-based
            pos = int(line[1]) + 1
            count = float(line[2]) * scaled_library_size / mapped_library_size
        except Exception as error:
            # in the event of a malformed BED file line during the wig 
            # conversion, this exception will be thrown
            print "Error from BAM to BED: %s, %s" % ( current_chrom, wig_fn )
            #import pdb
            #pdb.set_trace()
            return

        # this if statement starts the new section for each chromosome
        if not current_chrom == chrom:
            print "\t...%s..." % (current_chrom)
            chrom = current_chrom
            scaled_wig_fh.write('variableStep chrom=%s\n' % (chrom))
        # and this writes the coordinates/value for each line
        out_line = '\t'.join([str(pos), str(count)])
        scaled_wig_fh.write(out_line + '\n')
        #lastline = line

    scaled_wig_fh.close()
    wig_fh.close()

    return scaled_wig_fn

def makeBigWigFromWIG(wig_fn, chromosome_size_fn):
    '''
    Create BigWig
    '''
    bw_fn = wig_fn.replace('.wig', '.bw')
    print "Creating BigWig file %s..." % (bw_fn)
    cmd_args = ['bin/wigToBigWig', wig_fn, chromosome_size_fn, bw_fn]
    proc = subprocess.Popen(cmd_args)
    proc.wait()

    return bw_fn

def makeBedGraphFromBigWig(bw_fn):
    '''
    Create BedGraph
    '''
    bg_fn = bw_fn.replace('.bw', '.bedGraph')
    print "Creating BedGraph file %s..." % (bg_fn)
    cmd_args = ['bin/bigWigToBedGraph', bw_fn, bg_fn]
    #print cmd_args
    proc = subprocess.Popen(cmd_args)
    proc.wait()

    return bg_fn
