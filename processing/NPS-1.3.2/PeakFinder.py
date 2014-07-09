#!/home/liulab/shin/python/bin/python
# 
# PeakFinder.py detects peaks from ChIP-sequence data using LoG, which originally senses the edges given a signal.

import sys,  math, operator, warnings #numpy, scipy
#from scipy.stats import *
#from rpy import *
import Prob
from Prob import poisson_cdf as pcdf
from Prob import poisson_cdf_inv as pcdfi

# a debuggin code to check the p value calc is all right
DEBUG=False

class PeakFinder:

    def __init__(self, new, pos, par):
 
        if len(new) == 0 or len(pos) == 0:
            raise ValueError("Empty NEW or POSITION")
        
        # new and position
        self.new = new
        self.pos = pos
        
        # general attributes
        self.genomeLen = 2.4e9        # accessible genome length
        self.chr = pos[0].strip().split()[0]    # chromosome
        self.pc = 0                    # peak count
        
        # parameters for peak finding
        # interval: ChIP sequence resolution, default 10 bps
        if par['INTERVAL'] == '':
            par['INTERVAL'] = 10
        self.interval = int(par['INTERVAL'])
        
        # genome legnth
        if par['GENOME_LEN'] == '':
            par['GENOME_LEN'] = 2.4e9
            
        self.genomeLen = float(par['GENOME_LEN'])
        
        # p value: the cut-off p value for peak finding
        if par['PVALUE'] == '':
            par['PVALUE'] = 1e-5
        self.pvalCutOff = int(-100*math.log10(float(par['PVALUE'])))/10.0
       
        # extension: the tag length
        if par['EXTENSION'] == '':
            par['EXTENSION'] = 75
        self.tagLen = int(par['EXTENSION'])
       
        # tag number: the number of tags
        if par['TAG_NUM'] == '':
            par['TAG_NUM'] = 186e6
        self.tagNum = float(par['TAG_NUM'])
        
        # minimum width: the minimum peak width in peak finding
        if par['MIN_WIDTH'] == '':
            par['MIN_WIDTH'] = 80
        self.minWid = int(par['MIN_WIDTH'])
        
        # maximum width: the maximum peak width in peak finding 
        if par['MAX_WIDTH'] == '':
            par['MAX_WIDTH'] = 250
        self.maxWid = int(par['MAX_WIDTH'])
        
        # maximum height: the maximum peak height in peak finding
        if par['MAX_HEIGHT'] == '':
            par['MAX_HEIGHT'] = 10000
        self.maxHeight = int(par['MAX_HEIGHT'])
        
        # minimum height: the minimum peak height in peak finding
        # if it is not given, it is automatically calculated based on the other parameters
        if par['MIN_HEIGHT'] == '':
            par['MIN_HEIGHT'] = self.calcLowerCutOff(pVal = 1e-4)
        self.minHeight = int(par['MIN_HEIGHT'])
        print self.calcLowerCutOff(pVal=1e-4), 'is the lower cutoff'
        print self.minHeight, 'is the min height'

        # LoG: Laplacian of Gaussian, the width of the Gaussian envelope
        if par['LOG'] == '':
            par['LOG'] = -1 
        self.LoG = int(par['LOG'])
        
        # slope: The first derivative of Gaussian, the width of the Gaussian envelope
        if par['SLOPE'] == '':
            par['SLOPE'] = -1 
        self.slope = int(par['SLOPE'])
        
        # peak to inflection ratio: the ratio of the peak height and inflection point height
        if par['PEAK_INFLECTION_RATIO'] == '':
            par['PEAK_INFLECTION_RATIO'] = 1.2 
        self.peak2inflection = float(par['PEAK_INFLECTION_RATIO'])
        
        # median window size: the size of window in median filtering
        if par['MED_WSIZE'] == '':
            par['MED_WSIZE'] = 5
        self.medWsize = int(par['MED_WSIZE'])
              
        

    ###################################################
    # Laplacian of Gaussian filter; sigma = STD
    #
    def LOG(self,sigma):
        mask = []
        sigma2 = sigma*sigma*1.0
        for x in range(-sigma*3,sigma*3+1):
            mask.append((x*x/sigma2 - 1 ) *math.exp(-x*x/(2.0*sigma2)))
        return mask
    
    ###################################################
    # Simple Laplican operator
    #
    
    def Laplace(self):
        LapOp = [1, -2, 1]      
 #       LapOp = [5.0/84.0, 0.0, -3.0/84.0, -4.0/84.0, -3.0/84.0, 0.0, 5.0/84.0]
        
        return LapOp
    
    ###################################################
    # Simple Gradient operator
    #
    
    def Gradient(self):
        GradOp = [0.5, 0.0, -0.5]
#        GradOp = [3.0/28.0, 2.0/28.0, 1.0/28.0, 0.0, -1.0/28.0, -2.0/28.0, -3.0/28.0]
        
        return GradOp

    ###################################################
    # Slope of Gaussian
    #
    def GaussSlope(self,sigma):
        mask = []
        sigma2 = sigma*sigma*1.0
        for x in range(-sigma*3,sigma*3+1):
            mask.append(2.0*x*math.exp(-x*x/(2.0*sigma2)))
        #norm = sum(mask)
        #mask = [ x/norm for x in mask]
        return mask
   
    ###################################################
    # Median
    # in case a even number of data points, the average of two n/2 th and n/2 + 1th points is returned.
    #
    
    def median(self, data):
        
        temp = data[:]
        temp.sort()
        dataLen = len(data)
        if dataLen % 2 == 0: # even number of data points
            med = (temp[dataLen/2 -1] + temp[dataLen/2])/2.0
        else:
            med = temp[dataLen/2]
        
        return med
   
    ###################################################
    # medFilter
    # median filtering based on a sliding window
    # 
    def medFilter(self, seq, wsize):
                
        if wsize != 0:
            seqLen = len(seq)           
            padded = self.padData(seq, wsize/2)
           
            newseq = []
            for i in xrange(0, seqLen):
                # median filtering
                s = self.median(padded[i:i+wsize]) 
            
                #save the filtered value in a new array
                newseq.append(s)
            
            return newseq
        else: # if wsize is 0, do not median filtering
            return seq    
   
    
    ###################################################
    # CONVOLUTION
    # padded = padded version of array
    # filter = symmetric filter, e.g. LOG
    #
    def convolve(self, arraylen, padded, filter):
        temp = []
        L = len(filter)
        for x in xrange(arraylen):
            newvalue = 0.0 
            for y in xrange(L):
               newvalue += (padded[x+y]*filter[y])
                               
            temp.append(newvalue)
        return temp

    ###################################################
    # Pad the edges of the data so that we can apply LOG
    #
    def padData(self, array, padlen):
        newdata=[]
        PAD = [0 for x in range(padlen)]
        newdata.extend(PAD)
        newdata.extend(array)
        newdata.extend(PAD)
        return newdata
    
    ###################################################
    # Get the peak to trough ratio of a peak
    #
    def getRatio(self, peak, trough):
        if trough == 0:
            peak2trough = (peak + 0.1)/ (trough + 0.1)
        elif trough < 0:
            peak2trough = (peak - trough + 0.1)/0.1
        else:
            peak2trough = (peak/trough)
        
        return peak2trough
    
    ###################################################
    # Calculate the pvalue of each peak (return log pvalue)
    #
    
    def CalcPval(self, ratios, peakLoc):
        
        # calculate the null probability
        left = peakLoc[0]
        right = peakLoc[1]
        width = right-left+1
        l = self.tagLen/self.interval
        p = 1.0*width/(self.genomeLen/self.interval)
        LAMBDA = self.tagNum*p
        
        # linear interpolation for calculating the effect of tags on the peak
        start = max(0, left-l/2)
        end = min(len(ratios)-2, right+l/2)
        
        # weight
        weightL = range(start-left+l/2,l)
        for i in range(0, len(weightL)):
            weightL[i] = 1.0*weightL[i]/l
        weightR = range(l-1, right+l/2-end-1, -1)
        for i in range(0, len(weightR)):
            weightR[i] = 1.0*weightR[i]/l
        
        weight=[]
        # when peak width is smaller than tag extension, a warning is raised
        if width<=l:
            warnings.warn("Peak width, %d bp is recommended to be larger than tag extension (EXTENSION in *.par)" %(self.interval*width))
            weight.extend(weightL[:(end-start+1)/2])
            weight.extend(weightR[-1*(end-start+1)/2:])
            
            diff = end-start+1 - len(weight)
            if diff < 0:
                diff = -1*diff
                lp = diff/2
                rp = diff - lp
                r = [0]*lp + ratios[start:(end+1)] + [0]*rp 
                w = weight
            elif diff == 0:
                r = ratios[start:(end+1)]
                w = weight
            else:
                lp = diff/2
                rp = diff - lp
                r = ratios[start:(end+1)]
                w = [0]*lp + weight + [0]*rp
            
            tagCount = sum( map( lambda x, y: x*y, w, r ) )
            tagCount = 1.0*tagCount/l
            #tagCount = 0
            #for i in range(start):
            #    tagCount += weight[i-start]*ratios[i]
            #tagCount = 1.0*tagCount/l

        else:
            weightM = [1]*(end-start+1-len(weightL)-len(weightR))
            weight.extend(weightL)
            weight.extend(weightM)
            weight.extend(weightR)
            
            tagCount =sum( map(lambda x, y: x*y, weight, ratios[start:end+1] ) )
            tagCount = 1.0*tagCount/l
            #tagCount = 0
            #for i in range(start, end+1):
            #    tagCount += weight[i-start]*ratios[i]
            #tagCount = 1.0*tagCount/l
                        
        # count the tags within the bucket (W+l/2+l/2)
        #tagCount = 0
        #for i in range(start, end+1):
        #    tagCount += weight[i-start]*ratios[i]
        #tagCount = 1.0*tagCount/l      
        
         # calculate the P value using poison distribution: Log p value is used to cover very small p values
        #pval = poisson.sf(tagCount, LAMBDA)               # p value using scipy
#        pVal=1-r.ppois(tagCount,LAMBDA)
#        if pVal<1e-5:
#            if pVal==0.0: pVal=-1*r.ppois(tagCount,LAMBDA,log_p=r.T)    # using Taylor series: lp=r.ppois(log_p=r.T), 1-p=1-exp(lp)~= 1-(1+lp)=-lp
#            if pVal==0.0: pVal=3250.0
#            else:
#                pVal=int(-100*math.log10(pVal))/10.0
#        else:
#            pVal=int(-100*math.log10(pVal))/10.0

        # calculate P values using Prob.py. I put round(tagCount, 6) because the mismatch between R and Python in handling floating point numbers
        pVal = pcdf(round(tagCount,6), LAMBDA, lower=False)
#        if pVal == 0.0: pVal = int(-100*math.log10(3250.0))/10.0
#        else: pVal = int(-100*math.log10(pVal))/10.0
        if pVal == 0.0: 
            pVal = 3250.0
        else:
            pVal = int(-100*math.log10(pVal))/10.0
#        if abs(LAMBDA - 0.587626583333) <= 1e-5 and abs(tagCount - 6.0) <= 1e-5:
#            print l,tagCount
        return pVal, LAMBDA, tagCount        # this p value increases as the real p value decreases.
    
    ###################################################
    #
    # Calculate the lower cut off threshold based on a simple probabilistic modeling given the number of tags, the length of a tag, and the length of sequentiable genome
    
    def calcLowerCutOff(self, pVal=1e-4):
        
        LAMBDA = 1.0 *self.tagNum * self.tagLen / self.genomeLen 
        #return poisson.isf(1-pVal, LAMBDA)
        #return r.qpois(1-pVal,LAMBDA)
       
        # calculate the lower cut-off using Prob.py
        return pcdfi(1-pVal, LAMBDA, maximum=1000)
    
    ###################################################
    # Find the peak candidates based on LoG
    #
    
    def strictCandidates(self, subseq, ratios, logseq,  slopeseq, minFoldChange, maxFoldChange, minWid, maxWid, minPeak2Inflection, pvalCutoff):
        candidates = []
        #backtrack = int(50/self.interval)
        backtrack = 1
        LENGTH = len(subseq)
        FOUND = False       # True if negative 2nd derivative
        L = 3               # To compute the peak and trough heights, average the top and bottom L probes
        overL = 1.0/float(L)

        
        #********  initiation **********
        if logseq[0] < 0:
            start = 0
            FOUND = True
        
        x = 0
        while x < LENGTH:
            if FOUND and logseq[x] >=0:
                FOUND = False
                end = x
                boundary1 = start + 1 
                boundary2 = end -1  
                
                # search for a point closer to 0 in logseq
                if abs(logseq[start]) < abs(logseq[boundary1]):
                    boundary1 = start
                
                if abs(logseq[end]) < abs(logseq[boundary2]):
                    boundary2 = end

                # if two boundaries are too close, continue
                if boundary1 >= boundary2 or boundary2 - boundary1 <= 2:
                    continue
               
                #****** To compute the peak and inflection heights, average the top and bottom L probes
                temp = []
                for x in range(boundary1, boundary2+1):   # find peak and inflection
                    temp.append(ratios[x])
                temp.sort()
                
                
                peak = sum(temp[-L:])* overL
                inflectionL = sum(ratios[boundary1-L+1:boundary1+1])*overL
                inflectionR = sum(ratios[boundary2:boundary2+L])*overL
            
                inflectionRatio = self.getRatio(inflectionL, inflectionR)
                inflectionMin = min(inflectionL,inflectionR)
                inflectionMax = max(inflectionL,inflectionR)
                inflection = (inflectionL +inflectionR)/2       # take the mean of two inflections

                #****** Make sure that the peak is tall enough   ********
                peak2inflection = self.getRatio(peak, inflection)
                            
                if peak < minFoldChange or peak > maxFoldChange or peak2inflection < minPeak2Inflection: continue


               
                #****** Find left edge ***************
                LEDGE = boundary1
                REDGE = boundary2
                #****** Find right edge ***************
                
                LEFTPOS = subseq[LEDGE][0]
                RIGHTPOS = subseq[REDGE][0]
                tmp = []
                for x in xrange(LEDGE, REDGE+1):
                    tmp.append(ratios[x])
                PEAKPOS = subseq[tmp.index(max(tmp))+LEDGE][0]
                del tmp
                WIDTH = RIGHTPOS - LEFTPOS +1
                if ( WIDTH >= minWid) and  ( WIDTH <= maxWid):             
                    pval, LAMBDA, tc = self.CalcPval(ratios, [LEDGE, REDGE]) 
                    if pval >= pvalCutoff:
                        if DEBUG:
                            candidates.append("%s\t%d\t%d\t%d\t%e\t%f\t%f" %(self.chr, LEFTPOS, RIGHTPOS, self.pc, pval, LAMBDA, tc))# PEAKPOS, peak, inflectionL, inflectionR, pval])                            
                        else:
                            candidates.append("%s\t%d\t%d\t%d\t%e" %(self.chr, LEFTPOS, RIGHTPOS, self.pc, pval))# PEAKPOS, peak, inflectionL, inflectionR, pval])
                        self.pc += 1
                x = REDGE
            elif not FOUND and logseq[x] < 0:
                start = x-1
                FOUND = True
            x +=1 
            
        if FOUND and slopeseq[LENGTH-3] < 0:
            end = LENGTH-1
            boundary1 = start +1 
            boundary2 = end
                        
            # search for a point closer to 0 in logseq
            if abs(logseq[start]) < abs(logseq[boundary1]):
                boundary1 = start
                        
            if (subseq[boundary2][0] - subseq[boundary1][0] +1 >= minWid) :#and (boundary2 - boundary1 +1 <= maxWid):
                temp = []
                for x in range(boundary1, boundary2+1):
                    temp.append(ratios[x])
                temp.sort()

                L = 3
                peak = sum(temp[-L:])/float(L)
                
                inflectionL = sum(ratios[boundary1-L+1:boundary1+1])*overL
                inflectionR = sum(ratios[boundary2:boundary2+L])*overL
                
                inflection = (inflectionL +inflectionR)/2
                inflectionRatio = self.getRatio(inflectionL, inflectionR)
                inflectionMin = min(inflectionL,inflectionR)
                inflectionMax = max(inflectionL,inflectionR)
                LEDGE = boundary1
                REDGE = boundary2
                
                LEFTPOS = subseq[LEDGE][0]
                RIGHTPOS = subseq[REDGE][0]
                tmp = []
                for x in xrange(LEDGE, REDGE+1):
                    tmp.append(ratios[x])
                PEAKPOS = subseq[tmp.index(max(tmp))+LEDGE][0]
                del tmp
                WIDTH = RIGHTPOS - LEFTPOS +1
                if ( WIDTH  >= minWid) and ( WIDTH <= maxWid):
                    peak2inflection = self.getRatio(peak, inflection)        # need to think about. Is this OK?
                    if peak >= minFoldChange and peak <= maxFoldChange and peak2inflection >= minPeak2Inflection:
                        pval, LAMBDA, tc = self.CalcPval(ratios, [LEDGE, REDGE])           

                        if pval >= pvalCutoff:
                            if DEBUG:
                                candidates.append("%s\t%d\t%d\t%d\t%e\t%f\t%f" %(self.chr, LEFTPOS, RIGHTPOS, self.pc, pval, LAMBDA, tc))# PEAKPOS, peak, inflectionL, inflectionR, pval])
                            else:
                                candidates.append("%s\t%d\t%d\t%d\t%e" %(self.chr, LEFTPOS, RIGHTPOS, self.pc, pval))# PEAKPOS, peak, inflectionL, inflectionR, pval])
                            self.pc += 1
        return candidates
    
                                         
       
    ###################################################
    # slopeLevel = level of Gaussian in finding slopes
    # LoG = 2, 3, 4, 5 or 6
    # minWid = minimum width of peaks
    # cutoff = minimum absolute height of peaks
    # maxHeight = maximum absolute height of peaks
    # minHeight = min height of peak compared to inflection

    def findPeaks(self):
        
        peakList = []
        
        for pln in self.pos:
            p = pln.strip().split()
            subseq = [[int(p[1]) + j*self.interval, self.new[int(p[3])-1 + j]] for j in xrange(0, int(p[4])-int(p[3])+1)]
            
            if self.slope != -1:    # if SLOPE was not given, no Gaussian smoothing
                slope = self.GaussSlope(self.slope)
            else:
                slope = self.Gradient()
                
            # get the padding length for slope
            L2S = (len(slope)-1)/2
                
            if self.LoG != -1:
                log1 = self.LOG(1)
                log2 = self.LOG(2)
                log3 = self.LOG(3)
                log4 = self.LOG(4)
                log5 = self.LOG(5)
                log6= self.LOG(6)
                if self.LoG > 6:
                    logCustom = self.LOG(self.LoG)
                    LCustom = (len(logCustom)-1)/2
                L1 = (len(log1)-1)/2
                L2 = (len(log2)-1)/2
                L3 = (len(log3)-1)/2
                L4 = (len(log4)-1)/2
                L5 = (len(log5)-1)/2
                L6 = (len(log6)-1)/2
            else:
                laplace = self.Laplace()
                # get the padding length for Laplacian
                LL = (len(laplace)-1)/2
            
            arraylen = len(subseq)
            ratios = [0]*arraylen 
            Slope = [0]*arraylen 
            
            if self.LoG != -1:
                LoG1 = [0]*arraylen
                LoG2= [0]*arraylen 
                LoG3 = [0]*arraylen
                LoG4 = [0]*arraylen
                LoG5 = [0]*arraylen
                LoG6 = [0]*arraylen
                if self.LoG > 6:
                    LoGCustom = [0]*arraylen
            else:
                Laplace = [0]*arraylen
                
         
            for x in range(len(subseq)):
                ratios[x]  = subseq[x][1]
                        
            padded = self.padData(ratios,L2S)
            Slope = self.convolve(arraylen,padded,slope)
            
            if self.LoG == 6:
                ratios = self.medFilter(ratios, self.medWsize)
                padded6 = self.padData(ratios,L6)
                LoG6 = self.convolve(arraylen,padded6,log6)
                candidates = self.strictCandidates(subseq, ratios, LoG6,  Slope, \
                                                   self.minHeight,self.maxHeight, self.minWid, self.maxWid, self.peak2inflection,self.pvalCutOff)

            elif self.LoG == 5:
                ratios = self.medFilter(ratios, self.medWsize)
                padded5 = self.padData(ratios,L5)               
                LoG5 = self.convolve(arraylen,padded5,log5)
                candidates = self.strictCandidates(subseq, ratios, LoG5,  Slope,\
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)
                

            elif self.LoG == 4:
                ratios = self.medFilter(ratios, self.medWsize)
                padded4 = self.padData(ratios,L4)
                LoG4 = self.convolve(arraylen,padded4,log4)
                candidates = self.strictCandidates(subseq, ratios, LoG4,  Slope, \
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)

            elif self.LoG == 3:
                ratios = self.medFilter(ratios, self.medWsize)
                padded3 = self.padData(ratios,L3)
                LoG3 = self.convolve(arraylen,padded3,log3)
                candidates = self.strictCandidates(subseq, ratios, LoG3,  Slope,\
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)

            elif self.LoG == 2:
                ratios = self.medFilter(ratios, self.medWsize)
                padded2 = self.padData(ratios,L2)
                LoG2 = self.convolve(arraylen,padded2,log2)
                candidates = self.strictCandidates(subseq, ratios, LoG2,  Slope,\
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)

            elif self.LoG == 1:
                ratios = self.medFilter(ratios, self.medWsize)
                padded1 = self.padData(ratios,L1)
                LoG1 = self.convolve(arraylen,padded1,log1)
                candidates = self.strictCandidates(subseq, ratios, LoG1,  Slope,\
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)

               
            elif self.LoG > 6:
                ratios = self.medFilter(ratios, self.medWsize)
                padded = self.padData(ratios,LCustom)
                LoGCustom = self.convolve(arraylen,padded,logCustom)
                
                candidates = self.strictCandidates(subseq, ratios, LoGCustom,  Slope,\
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)
 
            # in case we do not use Guassian smoothing
            elif self.LoG == -1:
                ratios = self.medFilter(ratios, medWsize)
                padded = self.padData(ratios,LL)
                Laplace = self.convolve(arraylen, padded, laplace)
                candidates = self.strictCandidates(subseq, ratios, Laplace,  Slope,\
                                                   self.minHeight,self.maxHeight, self.minWid,self.maxWid,self.peak2inflection,self.pvalCutOff)
 

            peakList.extend(candidates)
 
        return peakList 
     


if __name__ == "__main__":
    print >> sys.stderr, "Use this by importing it in another python file"
