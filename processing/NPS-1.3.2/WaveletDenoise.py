#!/home/liulab/shin/python/bin/python
# 
# WaveletDenoise.py
#
# WaveletDenoise.py contains functions for wavelet denoising of ChIP-seq data. 
#
# Programmed by Hyunjin Shin 
#

import pywt
import sys
from numpy import *
from math import *
from WDen import *
      
#################################
#
# wden(x, tptr, sorh, scal, n, wname) does wavelet denoising. 
#
# x, input signal to be denoised
# tptr, threshold selection rule. See thselect.
# sorh, threshold type. see wthresh
# scal = 'one', for no threshold rescaling
#      = 'sln', for rescaling using a single estimation of level noise based 
#               on the first detail coefficients
#      = 'mln', for rescaling done usig level dependent estimation
# wname, wavelet name
#

def denoiseChIPSeq(new, pos, par):
    
    
    # Decomposition level, default = 2
    if par['DECOMP_LEVEL'] == '':
        par['DECOMP_LEVEL'] = 2
    level = int(par['DECOMP_LEVEL'])
    
    # wavelet, default = coif4
    if par['WAVELET'] == '':
        par['WAVELET'] = 'coif4'
    wname = par['WAVELET'].lower()
    
    # thresholding estimation method, default = heursure
    if par['THRESHOLD_EST'] == '':
        par['THRESHOLD_EST'] = 'heursure'
    tptr = par['THRESHOLD_EST'].lower()
        
    # threshold type, default = soft
    if par['THRESHOLD_TYPE'] == '':
        par['THRESHOLD_TYPE'] = 'soft'
    sorh = par['THRESHOLD_TYPE'].lower()
    
    # rescaling, default = 'mln'
    if par['SCALE'] == '':
        par['SCALE'] = 'mln'
    scal = par['SCALE'].lower()
    
    # perform wavelet denoising on the chromosome
    denoised = []    
    for pln in pos:
        p = pln.strip().split()
        subseq = [new[int(p[3])-1 + j] for j in xrange(0, int(p[4])-int(p[3])+1)]
            
        dsubseq = wden(subseq, tptr, sorh, scal, level, wname)
        denoised.extend(dsubseq)
        
    return denoised
    
            
    