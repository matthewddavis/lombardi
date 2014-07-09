# ------------------------------------
# python modules
# ------------------------------------
from math import exp
from math import log
# ------------------------------------
# constants
# ------------------------------------
LSTEP = 200
EXPTHRES = exp(LSTEP)
EXPSTEP  = exp(-LSTEP)
# ------------------------------------
# import modules - 2
# ------------------------------------

import math

# ------------------------------------
# functions
# ------------------------------------

def poisson_cdf (n, lam,lower=True):
    """Poisson CDF evaluater.

    This is a more stable CDF function. It can tolerate large lambda
    value. While the lambda is larger than 700, the function will be a
    little slower.

    Parameters:
    n     : your observation
    lam   : lambda of poisson distribution
    lower : if lower is False, calculate the upper tail CDF
    """
    k = int(n)
    if lam <= 0.0:
        raise Exception("Lambda must > 0")

    if lower:
        if lam > 700:
            return __poisson_cdf_large_lambda (k, lam)
        else:
            return __poisson_cdf(k,lam)
    else:
        if lam > 700:
            return __poisson_cdf_Q_large_lambda (k, lam)
        else:
            return __poisson_cdf_Q(k,lam)

def __poisson_cdf (k,a):
    """Poisson CDF For small lambda. If a > 745, this will return
    incorrect result.

    """
    if k < 0:
        return 0                        # special cases
    next = exp( -a )
    cdf = next
    for i in xrange(1,k+1):
        last = next
        next = last * a / i
        cdf = cdf + next
    if cdf > 1:
        return 1
    else:
        return cdf
    
def __poisson_cdf_large_lambda ( k,a ):
    """Slower poisson cdf for large lambda.
    
    """
    if k < 0:
        return 0                        # special cases
    num_parts = int(a/LSTEP)
    last_part = a % LSTEP
    lastexp = exp(-last_part)
    next = EXPSTEP
    num_parts -= 1
    cdf = next
    for i in xrange(1,k+1):
        last = next
        next = last * a / i
        cdf = cdf + next
        if next > EXPTHRES or cdf > EXPTHRES:
           if num_parts>=1:
               cdf *= EXPSTEP
               next *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1

    for i in xrange(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

def __poisson_cdf_Q (k,a):
    """internal Poisson CDF evaluater for upper tail with small
    lambda.

    """
    if k < 0:
        return 1                        # special cases
    next = exp( -a )

    for i in xrange(1,k+1):
        last = next
        next = last * a / i

    cdf = 0
    i = k+1
    while next >0:
        last = next
        next = last * a / i
        cdf += next
        i+=1
    return cdf

def __poisson_cdf_Q_large_lambda (k,a):
    """Slower internal Poisson CDF evaluater for upper tail with large
    lambda.
    
    """
    if k < 0:
        return 1                        # special cases
    num_parts = int(a/LSTEP)
    last_part = a % LSTEP
    lastexp = exp(-last_part)
    next = EXPSTEP
    num_parts -= 1

    for i in xrange(1,k+1):
        last = next
        next = last * a / i
        if next > EXPTHRES:
           if num_parts>=1:
               next *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1
    cdf = 0
    i = k+1
    while next >0:
        last = next
        next = last * a / i
        cdf += next
        i+=1
        if next > EXPTHRES or cdf > EXPTHRES:
           if num_parts>=1:
               cdf *= EXPSTEP
               next *= EXPSTEP
               num_parts -= 1
           else:
               cdf *= lastexp
               lastexp = 1

    for i in xrange(num_parts):
        cdf *= EXPSTEP
    cdf *= lastexp
    return cdf

def poisson_cdf_inv ( cdf, lam, maximum=1000):
    """inverse poisson distribution.

    cdf : the CDF
    lam : the lambda of poisson distribution

    note: maxmimum return value is 1000
    and lambda must be smaller than 740.
    """
    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception ("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0
    sum2 = 0
    newval = exp( -lam )
    sum2 = newval
#     if cdf <= sum2:
#         return i

    for i in xrange(1,maximum+1):
        sumold = sum2
#         if i == 0:
#             newval = exp( -a )
#             if newval==0:
#                 newval = 4.9406564584124654e-324
#             sum2 = newval
#         else:
        last = newval
        newval = last * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    
    return maximum

def poisson_cdf_Q_inv ( cdf, lam, maximum=1000):
    """inverse poisson distribution.

    cdf : the CDF
    lam : the lambda of poisson distribution

    note: maxmimum return value is 1000
    and lambda must be smaller than 740.
    """
    assert lam < 740
    if cdf < 0 or cdf > 1:
        raise Exception ("CDF must >= 0 and <= 1")
    elif cdf == 0:
        return 0
    sum2 = 0
    newval = exp( -lam )
    sum2 = newval

    for i in xrange(1,maximum+1):
        sumold = sum2
        last = newval
        newval = last * lam / i
        sum2 = sum2 + newval
        if sumold <= cdf and cdf <= sum2:
            return i
    
    return maximum

def poisson_pdf ( k, a ):
    """Poisson PDF.

    PDF(K,A) is the probability that the number of events observed in
    a unit time period will be K, given the expected number of events
    in a unit time.
    """
    if a <= 0:
        return 0
    return exp(-a) * pow (a, k) / factorial (k)

