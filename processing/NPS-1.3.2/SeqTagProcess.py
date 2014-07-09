#SeqTagProcess.py
# zy@jimmy.harvard.edu, shin@jimmy.harvard.edu

import sys, os.path, copy, time

# add the path to SeqTag
#sys.path.append(r'/Users/jin/Research/lib/NPS-1.0.1')
#sys.path.append(r'/home/shin/Research/lib/SeqTag')

from numpy import *
from GenomeInfo import chrlength
from Bed import *
from ReadParam import *
from WaveletDenoise import *
#from PeakFinder import *

def mainloop(sortedbedfile, par = {}):
    """mainloop(sortedbedfile, parameter_hash)
    The parameters should be read from parameter file
    for example:
    par['SPENAME'] = 'hg18'
    par['EXTENSION'] = 75
    par['SHIFT'] = 37
    par['PVALUE'] = 1e-5
    par['LOG] = 3
    par['MIN_WIDTH'] = 80
    par['MAX_WIDTH] = 250
    par['MIN_HEIGHT'] = 17
    par['PEAK_INFLECTION_RATIO'] = 1.2
    par['BIAS_RATIO'] = 4
    par['OUTFILE'] = filename
    par['TAG_THR'] = tag threshold. If tag pile is less than this threshold over an 1kb region, this is ignored for the further processing.
    and some other parameters needed
    """
    
    chrold = ''
    cluster = []
    plus = []
    wigfile = open(par['OUTFILE'], 'w')
    #outfilename=os.path.join(os.path.split(sortedbedfile)[0],os.path.split(sortedbedfile)[1].split('.bed')[0]+'.wig')
    #wigfile = open(outfilename, 'w')
    #wigfile = open(par['OUTFILE'], 'w')
    #wigfile=open(par['INFILE'].split('.')[0]+'_wd.wig','w')
    tag_thr = int(par['TAG_THR'])
    
    for line in open(sortedbedfile).xreadlines():
        if not line[:3] == 'chr':
            continue
        line = line.strip().split()
        if not chrlength[par['SPENAME']].has_key(line[0]):
            continue
        
        if chrold == '':
            chrold = line[0]
            print >>sys.stderr, 'Reading for', chrold, '......', time.asctime()
            if chrlength[par['SPENAME']][line[0]] % 10 == 0:
                cluster = [0] * (chrlength[par['SPENAME']][line[0]] / 10)       # Tag number (10 bp space)
                plus = [0] * (chrlength[par['SPENAME']][line[0]] / 10)          # Tag number in plus strand (10 bp space)
            else:
                cluster = [0] * (chrlength[par['SPENAME']][line[0]] / 10 + 1)
                plus = [0] * (chrlength[par['SPENAME']][line[0]] / 10 + 1)
            
        elif line[0] != chrold:
            print >>sys.stderr, 'Changing for', chrold, '......', time.asctime()
            package = {}
            package = change2postion(chrold, cluster, tag_thr) # package['position_file'], package['new_file']
#            change2postion_saveas_wig(wigfile,chrold, cluster)
            ##################################
            #    test
            ##################################
            #print >>sys.stderr, 'writing for', chrold, '......', time.asctime()
            #filetmp = open(chrold + '.new', 'w')
            #for c in package['new_file']:
            #    print >>filetmp, c
            #filetmp.close()
            #filepos = open(chrold + '.pos', 'w')
            #for p in package['position_file']:
            #    print >>filepos, p
            #filepos.close()
            #peakregions = []
            #for m in open('/home/liulab/yzhang/ChIP_seq/num_pos/nucleosome_cluster_75/Nucleosome_chr1.cluster_10dec_wavelet.bed').xreadlines():
            #    if m[:3] == 'chr':
            #        peakregions.append(m.strip())
            ##################################
            
            ##################################
            #    to be finished
            ##################################
            if par['WANT_DENOISE'] == 'yes':
                print >>sys.stderr, 'Denoising for', chrold, '......', time.asctime()
                denoised = denoiseChIPSeq(package['NEW_FILE'], package['POSITION_FILE'], par)
            else:
                denoised = package['NEW_FILE']
            
            print>>sys.stderr, 'Saving as a wig file', chrold, '......', time.asctime()
            change2wig(wigfile,package['POSITION_FILE'],denoised)

            package = {}
            peakregions = []
            peakregions_filtered  = []
            ##################################
            chrold = line[0]
            print >>sys.stderr, 'Reading for', chrold, '......', time.asctime()
            
            if chrlength[par['SPENAME']][line[0]] % 10 == 0:
                cluster = [0] * (chrlength[par['SPENAME']][line[0]] / 10)       # Tag number (10 bp space)
                plus = [0] * (chrlength[par['SPENAME']][line[0]] / 10)          # Tag number in plus strand (10 bp space)
            else:
                cluster = [0] * (chrlength[par['SPENAME']][line[0]] / 10 + 1)
                plus = [0] * (chrlength[par['SPENAME']][line[0]] / 10 + 1)
        
        try:
            if line[5] == '+':           # Tag in plus strand
                b = int(line[1]) + int(par['SHIFT'])
                e = b + int(par['EXTENSION'])
                if (max(b, 1) - 1) / 10 == 0:
                    beginpos = max(b, 1)
                else:
                    beginpos = (max(b, 1) / 10 + 1) * 10 + 1
                for k in xrange(beginpos, min(e, chrlength[par['SPENAME']][line[0]]), 10):
                    cluster[(k - 1) / 10] += 1
                    plus[(k - 1) / 10] += 1
                                
            elif line[5] == '-':        # Tag in minus strand
                e = int(line[2]) - int(par['SHIFT'])
                b = e - int(par['EXTENSION'])
                if (max(b, 1) - 1) / 10 == 0:
                    beginpos = max(b, 1)
                else:
                    beginpos = (max(b, 1) / 10 + 1) * 10 + 1
                for k in xrange(beginpos, min(e, chrlength[par['SPENAME']][line[0]]), 10):
                    cluster[(k - 1) / 10] += 1
            else:
                continue
        except:
            print >> sys.stderr, 'Tag position file error: ', sys.exc_info()[0], sys.exc_info()[1]
            sys.exit()
    
    print >>sys.stderr, 'changing for', chrold, '......', time.asctime()
    package = change2postion(chrold, cluster, tag_thr) # package['position_file'], package['new_file']
#    change2postion_saveas_wig(wigfile,chrold, cluster)
    ##################################
    #    test
    ##################################
    #print >>sys.stderr, 'writing for', chrold, '......', time.asctime()
    #filetmp = open(chrold + '.new', 'w')
    #for c in package['new_file']:
    #    print >>filetmp, c
    #filetmp.close()
    #filepos = open(chrold + '.pos', 'w')
    #for p in package['position_file']:
    #    print >>filepos, p
    #filepos.close()
    ##################################
    
    ##################################
    #    to be finished
    ##################################
    if par['WANT_DENOISE'] ==  'yes':
        print >>sys.stderr, 'Denoising for', chrold, '......', time.asctime()
        denoised = denoiseChIPSeq(package['NEW_FILE'], package['POSITION_FILE'], par)
    else:
        denoised = package['NEW_FILE']

    print>>sys.stderr, 'Saving as a wig file', chrold, '......', time.asctime()
    change2wig(wigfile,package['POSITION_FILE'],denoised)

    package = {}
    peakregions = []
    peakregions_filtered  = []
    ##################################
    cluster = []
    plus = []
    wigfile.close()

def change2postion(chrname,cluster, tag_thr=2):
    """
    change2position(chrname, cluster)
    change cluster format to new & position format
    input: cluster = [0, 1, 1, ...] (cluster format, 10bp resolution)
    output: package
            package['position_file'] = ['chrY	2709761	2710961	1	121', 'chrY	2712321	2732421	122	2132'] (position format)
            package['new_file'] = [0, 1, 1, ...] (new format, 10bp resolution)
    """
    package = {}
    package['POSITION_FILE'] = []
    package['NEW_FILE'] = []
    
    cutoff = 1000 / 10
    oceanbegin = 0        # ocean: tag num <= 2
    oceanflag = 1
    
    num = []
    for k in xrange(len(cluster)):
        num.append(cluster[k])
        
    for k in xrange(len(num) - 1):
        if num[k] > tag_thr:
            if oceanflag == 1:
                oceanflag = 0
                if (k - oceanbegin) >= cutoff:
                    oceanflag = 0
                    for m in xrange(oceanbegin, k):
                        num[m] = -1
                
        elif num[k] <= tag_thr and oceanflag == 0:
            oceanbegin = k
            oceanflag = 1
    if oceanflag == 1:
        for m in xrange(oceanbegin, len(num)):
            num[m] = -1

    linenum = 0
    islandflag = 0
    islandbegin = 0
    islandline = 0
    for k in xrange(len(num) - 1):
        if islandflag == 0 and num[k] > -1:
            islandflag = 1
            linenum += 1
            islandbegin = k + 1
            islandline = linenum
            package['NEW_FILE'].append(num[k])
        elif islandflag == 1 and num[k] > -1:
            package['NEW_FILE'].append(num[k])
            linenum += 1
        elif islandflag == 1 and num[k] == -1:
            package['POSITION_FILE'].append('\t'.join([chrname, str(islandbegin * 10 - 9) , str(k * 10 - 9), str(islandline), str(linenum)]))
            islandflag = 0

    if islandflag == 1:
        package['NEW_FILE'].append(num[len(num) - 1])
        linenum += 1
        package['POSITION_FILE'].append('\t'.join([chrname, str(islandbegin * 10 - 9), str(len(num) * 10 - 9), str(islandline), str(linenum)]))
    
    num = []
    return package

def change2postion_saveas_wig(wigfile,chrname, cluster, tag_thr=2):
    """
    change2position(chrname, cluster)
    change cluster format to new & position format
    input: cluster = [0, 1, 1, ...] (cluster format, 10bp resolution)
    output: package
            package['position_file'] = ['chrY    2709761    2710961    1    121', 'chrY    2712321    2732421    122    2132'] (position format)
            package['new_file'] = [0, 1, 1, ...] (new format, 10bp resolution)
    """
    package = {}
   
    cutoff = 1000 / 10
    oceanbegin = 0        # ocean: tag num <= 2
    oceanflag = 1
            
    num = []
    for k in xrange(len(cluster)):
        num.append(cluster[k])
    
    #put a header for each chromosome
    print >>wigfile,"track type=wiggle_0\nvariableStep chrom=%s span=%d" %(chrname,10)
    
    for k in xrange(len(num) - 1):
        if num[k] > tag_thr:
            if oceanflag == 1:
                oceanflag = 0
                if (k - oceanbegin) >= cutoff:
                    oceanflag = 0
                    for m in xrange(oceanbegin, k):
                        num[m] = -1
                
        elif num[k] <= tag_thr and oceanflag == 0:
            oceanbegin = k
            oceanflag = 1
    if oceanflag == 1:
        for m in xrange(oceanbegin, len(num)):
            num[m] = -1

    linenum = 0
    islandflag = 0
    islandbegin = 0
    islandline = 0
    for k in xrange(len(num) - 1):
        if islandflag == 0 and num[k] > -1:
            islandflag = 1
            linenum += 1
            islandbegin = k + 1
            islandline = linenum
            print >>wigfile, "%d\t%d" %(islandbegin*10-9,num[k])
            
        elif islandflag == 1 and num[k] > -1:
            linenum += 1
            print >>wigfile, "%d\t%d" %(k*10+1,num[k])
        elif islandflag == 1 and num[k] == -1:
            islandflag = 0

    if islandflag == 1:
        linenum += 1
        print >>wigfile, "%d\t%d" %(len(num)*10-9,num[len(num)-1])
    
    num = []

def change2wig(wigfile,position,new):
               
    # get the chromosome name from the POSITION_FILE.
    chrname=position[0].split()[0]
    print >>wigfile,"track type=wiggle_0\nvariableStep chrom=%s span=%d" %(chrname,1)
    
    itr= 0
    for pos in position:
        pos = pos.split()
        islandbeg=int(pos[1])
        islandend=int(pos[2])
        linebeg=int(pos[3])
        lineend=int(pos[4])
        
        i = 0
        segment=new[linebeg-1:lineend]
        for loc in xrange(islandbeg,islandend+1,10):
            print >> wigfile, "%d\t%d" %(loc,int(round(segment[i])))
            i+=1
            itr+=1
    
    # check
    if itr!=len(new):
        print >>sys.stderr, 'Mismatch between POSITION_FILE and NEW_FILE!'
        
    
def  ratio_filter(peakregions = [], cluster = [], plus = []):
    """ratio_filter(peakregions, cluster, plus)
    return peakregion_filtered
    Note: plus/minus ratio cutoff: 4 fold
    """
    peakregion_filtered = []
    for lineold in peakregions:
        line = lineold.strip().split()
        plusnum = 0
        minusnum = 0
        for k in xrange(int(line[1]) / 10, int(line[2]) / 10):
            plusnum += plus[k]
            minusnum += cluster[k] - plus[k]
        
        if minusnum == 0 or plusnum == 0:
            continue
        ratio = float(plusnum) / float(minusnum)
        if 0.25 < ratio < 4:
            peakregion_filtered.append('\t'.join([lineold, str(plusnum), str(minusnum)]))
    return peakregion_filtered
            
def test1():
    par = {}
    par['SPENAME'] = 'hg18'
    par['EXTENSION'] = 75
    par['SHIFT'] = 37
    par['OUTFILE'] = 'tmp'
    mainloop('HM_hg18.bed.sorted', par = par)
    
def usage():
    print >> sys.stderr, "USAGE: python SeqTagProcess.py parameter_file"
    #print >> sys.stderr, "NOTE: OUTFILE in the par file will be ignored in this program"
    
if __name__ == '__main__':
    
    if len(sys.argv) < 2: usage()
    else:
        f = ReadParam(sys.argv[1])
        PARAMS = f.readParam()
        
        # if 'WANT_SORT' is yes, then sort
        infile = PARAMS['INFILE']
        if PARAMS['WANT_SORT'].lower() == 'yes':
            infile_sort = sortBed(infile)
        else: infile_sort = infile
        
        mainloop(infile_sort, par = PARAMS)
        
        print "Used parameters for this analysis are ..."
        printParam(PARAMS)
