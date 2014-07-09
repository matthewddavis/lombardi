#!/home/liulab/yzhang/python/bin/python
# zy@jimmy.harvard.edu

import sys, os.path, time
from numpy import *

class Bed(object):
    """ Bed file object
    1) if the bed file is not very large:
    bed1 = Bed()
    bed1.readbed(bedfilename)
    bed1.bedsort()
    bed1.release_space()
    2) if the bed file is very large:
    bed1 = Bed()
    bed1.bedsort_by_chr(bedfilename)
    bed1.release_space()
    """
    def __init__(self):
        self.chroms = []
        self.chr = []
        self.begin = []
        self.end = []
        self.content = []
        self.bedfilename = None
        
    def readbed(self, bedfilename):
        """readbed(bedfilename)
        read bed file
        """
        self.bedfilename = bedfilename
        for lineold in open(bedfilename).xreadlines():
            if lineold[:3] == 'chr':
                line = lineold.strip().split()
                self.chr.append(line[0])
                if not (line[0] in self.chroms):
                    self.chroms.append(line[0])
                self.begin.append(line[1])
                self.end.append(line[2])
                self.content.append(lineold.strip())
                
        self.chr = char.array(self.chr)
        self.begin = array([int(i) for i in self.begin], int32)
        self.end = array([int(i) for i in self.end], int32)
        self.content = char.array(self.content)
        
    def bedsort(self, newfile = ''):
        """bedsort(newfilename = '')
        The default new file name is self.bedfilename + '.sorted'
        """
        self.chroms.sort(cmp)
        
        if newfile:
            ofn = newfile
        else:
            l = self.bedfilename.split('.')
            if len(l) == 2:
                ofn = '.'.join([l[0] + '_sorted', l[1]])
            else:
                ofn = '.'.join(['.'.join(l[0:-1]) + '_sorted', l[-1]])
        
        outputfile = open(ofn, 'w')
        
        for k in xrange(len(self.chroms)):
            segment = where(self.chr == self.chroms[k])[0]
            chrbegin = self.begin[segment]
            chrcontent = self.content[segment]
            
            if not segment.shape[0]:
                continue
            index = chrbegin.argsort()
            for i in index:
                print >>outputfile, chrcontent[i]
        
        return ofn
    
    def release_space(self):
        """release_space()
        release the memory
        """
        self.chroms = []
        self.chr = []
        self.begin = []
        self.end = []
        self.content = []
        self.bedfilename = None
    
    def bedsort_by_chr(self, bedfilename, newfile = ''):
        """bedsort_by_chr(bedfilename)
        If the bed file is larger than 0.5G (smaller than 10G), to use chr by chr sort
        """
        self.bedfilename = bedfilename
        if newfile:
            ofn = newfile
        else:
            l = self.bedfilename.split('.')
            if len(l) == 2:
                ofn = '.'.join([l[0] + '_sorted', l[1]])
            else:
                ofn = '.'.join(['.'.join(l[0:-1]) + '_sorted', l[-1]])
        
        outputfile = open(ofn, 'w')
            
        chrname = 'chr1'   # initial from chr1
        self.chroms = []
        self.chr = []
        self.begin = []
        self.end = []
        self.content = []
        
        print >>sys.stderr, 'Beginning with chr1 ....', time.asctime()
        
        filelist = {}
        tagnumbers = {}
        
        for lineold in open(bedfilename).xreadlines():
            if lineold[:3] == 'chr':
                lineold = lineold.strip()
                line = lineold.split()
                if chrname == line[0]:
                    if not tagnumbers.has_key(line[0]):
                        tagnumbers[line[0]] = 1
                    else:
                        tagnumbers[line[0]] += 1
                    self.begin.append(line[1])
                    self.content.append(lineold)
                else:
                    if not filelist.has_key(line[0]):
                        if os.path.dirname(bedfilename):
                            filetmp = os.path.dirname(bedfilename) + os.path.sep + '.' + line[0] + '.bed'
                        else:
                            filetmp = '.' + line[0] + '.bed'
                        if os.access(filetmp, os.F_OK):       # to delete the tmp file, if it exists
                            os.unlink(filetmp)
                        filelist[line[0]]  = open(filetmp, 'a')
                        tagnumbers[line[0]] = 1
                        print >>filelist[line[0]], lineold
                    else:
                        tagnumbers[line[0]] += 1
                        print >>filelist[line[0]], lineold
        
        for chrname in filelist.keys():
            filelist[chrname].close()
                    
        if len(self.begin) > 0:
            print >>sys.stderr, 'Sorting chr1 ....', time.asctime()
            self.begin = array([int(i) for i in self.begin], int32)
            self.content = char.array(self.content)
            index = self.begin.argsort()
            for i in index:
                print >>outputfile, self.content[i]
        
        #backitems = [ [v[1],v[0]] for v in chrhash.items()]         # to order the chrnames by tag numbers (reversely)
        #backitems.sort(reverse = True)
        #sortedchr = [ backitems[i][1] for i in range(0,len(backitems))]
        
        sortedchr = filelist.keys()
        sortedchr.sort(cmp)
        
        for k in xrange(len(sortedchr)):
            chrname = sortedchr[k]
            if os.path.dirname(bedfilename):
                filetmp = os.path.dirname(bedfilename) + os.path.sep + '.' + chrname + '.bed'
            else:
                filetmp ='.' + chrname + '.bed'
            self.begin = []
            self.content = []
            print >>sys.stderr, 'Handling', chrname, '....', time.asctime()
            for lineold in open(filetmp).xreadlines():
                if lineold[:3] == 'chr':
                    lineold = lineold.strip()
                    line = lineold.split()
                    if chrname == line[0]:
                        self.begin.append(line[1])
                        self.content.append(lineold)
            if len(self.begin) > 0:
                print >>sys.stderr, 'Sorting', chrname, '....', time.asctime()
                self.begin = array([int(i) for i in self.begin], int32)
                self.content = char.array(self.content)
                index = self.begin.argsort()
                for i in index:
                    print >>outputfile, self.content[i]
            os.unlink(filetmp)
        outputfile.close()
        return ofn, tagnumbers

def sortBed(bedfilename):
    """sortBed(bedfilename)
    If the input file is too large (>10G), the program cann't sort the bed file. Users can try other ways to sort the bed file.
    If the size is larger than 0.5G, the program use chr by chr sort instead of sorting all chr regions together.
    """
    filesize = os.path.getsize(bedfilename)
    if filesize > 10000000000:
        print >>sys.stderr, 'The file is too large to be sorted!'
        sys.exit(0)
    elif filesize < 500000000:
        print >>sys.stderr, 'Sorting bed file ....'
        print >>sys.stderr, 'Start time:', time.asctime()
        bed = Bed()
        bed.readbed(bedfilename)
        ofn = bed.bedsort()
        bed.release_space()
        print >>sys.stderr, 'End time:', time.asctime() 
    else:
        tagnumbers = {}
        bed = Bed()
        ofn, tagnumbers = bed.bedsort_by_chr(bedfilename)
        bed.release_space()
        print >>sys.stderr, 'End time:', time.asctime()
        #print tagnumbers
    
    return ofn
    
if __name__ == '__main__':
    
    if len(sys.argv) < 2:
        print "Usage: %prog bedfile"
        sys.exit(0)
        
    sortBed(sys.argv[1])
    