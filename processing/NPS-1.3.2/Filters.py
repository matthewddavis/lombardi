#!/usr/bin/python
# 
# Filters.py
#
# A collection of filters
#
# moving median filter
# moving average filter
#
# shin@jimmy.harvard.edu
#

import sys,  math, operator

              
###################################################
# Laplacian of Gaussian filter; sigma = STD
#
def LOG(sigma):
    mask = []
    sigma2 = sigma*sigma*1.0
    for x in range(-sigma*3,sigma*3+1):
        mask.append((x*x/sigma2 - 1 ) *math.exp(-x*x/(2.0*sigma2)))
    return mask
    
###################################################
# Simple Laplican operator
#
    
def Laplace():
    LapOp = [1, -2, 1]      
#       LapOp = [5.0/84.0, 0.0, -3.0/84.0, -4.0/84.0, -3.0/84.0, 0.0, 5.0/84.0]
        
    return LapOp
    
###################################################
# Simple Gradient operator
#
    
def Gradient():
    GradOp = [0.5, 0.0, -0.5]
#    GradOp = [3.0/28.0, 2.0/28.0, 1.0/28.0, 0.0, -1.0/28.0, -2.0/28.0, -3.0/28.0]
        
    return GradOp

###################################################
# Slope of Gaussian
#
def GaussSlope(sigma):
    mask = []
    sigma2 = sigma*sigma*1.0
    for x in range(-sigma*3,sigma*3+1):
        mask.append(2.0*x*math.exp(-x*x/(2.0*sigma2)))
    #norm = sum(mask)
    #mask = [ x/norm for x in mask]
    return mask

###################################################
# Gaussian smoothing
#

def Gauss(sigma):
    envelope=[]
    sigma2=sigma*sigma*1.0
    for x in range(-sigma*3,sigma*3+1):
        envelope.append(math.exp(-x*x/(2.0*sigma2))/((2.0*math.pi)**0.5 *sigma))
    norm=sum(envelope)
    envelope=[x/norm for x in envelope]
    return envelope

###################################################
# CONVOLUTION
# padded = padded version of array
# filter = symmetric filter, e.g. LOG
#
def convolve(arraylen, padded, filter):
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
def padData(array, padlen):
    newdata=[]
    PAD = [0 for x in range(padlen)]
    newdata.extend(PAD)
    newdata.extend(array)
    newdata.extend(PAD)
    return newdata

###################################################
# Median
# in case a even number of data points, the average of two n/2 th and n/2 + 1th points is returned.
#
    
def median(data):
        
    temp = data[:]
    temp.sort()
    dataLen = len(data)
    if dataLen % 2 == 0: # even number of data points
        med = (temp[dataLen/2 -1] + temp[dataLen/2])/2.0
    else:
        med = temp[dataLen/2]
     
    return med

###################################################
# Mean
# return the mean of a given sequence of numbers
#
def mean(data):
    
    dataLen = len(data)
    
    return 1.0*sum(data)/dataLen
    
###################################################
# medFilter
# median filtering based on a sliding window
# 
def medFilter(seq, wsize):
                
    if wsize != 0:
        seqLen = len(seq)           
        padded = padData(seq, wsize/2)
       
        newseq = []
        for i in xrange(0, seqLen):
            # median filtering
            s = median(padded[i:i+wsize]) 
            
            #save the filtered value in a new array
            newseq.append(s)
            
        return newseq
    else: # if wsize is 0, do not do median filtering
        return seq    
    
###################################################
# meanFilter
# moving average (mean) filter
#
   
def meanFilter(seq, wsize):
    
    if wsize != 0:
        seqLen = len(seq)           
        padded = padData(seq, wsize/2)
       
        newseq = []
        for i in xrange(0, seqLen):
            # mean filtering
            s = mean(padded[i:i+wsize]) 
            
            #save the filtered value in a new array
            newseq.append(s)
            
        return newseq
    else: # if wsize is 0, do not do mean filtering
        return seq    
    
###################################################
# gaussianFilter
# Gaussian smoothing
#

def gaussianFilter(seq, sigma):
    
    gauss=Gauss(sigma)
    padded=padData(seq,len(gauss)/2)
    newseq=convolve(len(seq),padded,gauss)
    
    return newseq


if __name__ == "__main__":
    print >> sys.stderr, "Use this by importing it in another python file"