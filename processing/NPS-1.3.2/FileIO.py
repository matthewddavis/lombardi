# jssong
# 12/19/2006

import math, StringUtil
import sys, copy

class FileIO:
    ''' Read nuclesome data'''

    def __init__(self, filename, filename_pos):
        self.filename = filename
        self.filename_pos = filename_pos
       
    ###################################################
    # READ DATA
    #  Read a table of wavelet denoised signals and find contiguous regions.
    #  new = [enrichment1, enrichment2, ...]
    #  position = [[start in chromosome, end in chromosome, start index, end index], ...]


    def readData(self):
        ''' sampleList = [ sampleName1, sampleName2, ...]
               spacing = tiled spacing
        '''
        # read the position file
        f= open(self.filename_pos, 'r')
        
        se = []
        FIRSTLINE = True
        for line in f.xreadlines():
            l = line.strip().split()
            if FIRSTLINE == True:
                if l[0] == "Spacing":
                    spacing = int(l[1])
                    
                else:
                    spacing = 1 
                    se.append([int(l[1]), int(l[2]), int(l[3]), int(l[4])])
                    chr = l[0]
                FIRSTLINE = False
                continue
            chr = l[0]
            se.append([int(l[1]), int(l[2]), int(l[3]), int(l[4])])
        
        f.close()
        
        # read the new file

        f = open(self.filename, 'r')
        
        enrichment = []
        for line in f.xreadlines():
            l = line.strip()
            enrichment.append(float(l))

        seqs = []    
        for i in xrange(0, len(se)):
            subseq = [[se[i][0]-1 + j*spacing, enrichment[se[i][2]-1 + j]] for j in xrange(0, se[i][3]-se[i][2]+1)]
            seqs.append(subseq)
            
        f.close()
        
        return chr, seqs        
       


if __name__ == "__main__":
    print >> sys.stderr, "Use in another python script"