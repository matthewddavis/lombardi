def getChromSizes(species, load_heterochromatin=False):
    '''
    Loads the chromosome lengths for each chromosome from a local file (fetched from UCSC).
    Returns them in in a dict keyed by chromosome.
    '''

    if species == 'sacCer3':
        chrom_size_fn='../reference/sacCer3.chrom.sizes'

    try:
        chrom_size_fh = open(chrom_size_fn, 'r')
    except:
        print "Error opening the chromosome size file %s." % chrom_size_fn

    chrom_size_dict = {}
    for line in chrom_size_fh:
        line = line.strip('\n').split('\t')
        chrom = line[0]
        length = int(line[1])
        chrom_size_dict[chrom] = length

    return chrom_size_dict


def loadChromScoreVectorDict(bedGraph_fn, species='sacCer3'):
    '''
    Accepts a bedGraph filename and loads all chromosomes into float vectors
    of the scores across the bedGraph regions.
    '''
    from numpy import zeros

    try:
        bedGraph_fh = open(bedGraph_fn, 'r')
    except:
        print "Error opening bedGraph file %s." % (bedGraph_fn)

    print "Loading bedGraph scores from %s..." % (bedGraph_fn)

    chrom_size_dict = getChromSizes(species)
    chrom_score_vector_dict = {}
    for chrom in chrom_size_dict:
        chrom_score_vector_dict[chrom] = zeros(chrom_size_dict[chrom], dtype=float)

    for line in bedGraph_fh:
        if line.startswith('track'):
            continue
        line = line.strip('\n').split('\t')
        chrom = line[0]
        start = int(line[1]) - 1  # bed format is 1-indexed
        stop = int(line[2]) - 1
        score = float(line[3])
        if chrom in chrom_score_vector_dict:
            chrom_score_vector_dict[chrom][start:stop] = score

    print "\t...done."
    return chrom_score_vector_dict

def getColors(num_colors, lightness=.5, saturation=1.0):
    '''
    stolen/adapted from http://stackoverflow.com/q/470690
    '''
    import colorsys
    from numpy import arange

    if num_colors == 0:
        num_colors = 1

    colors=[]
    for i in arange(0., 360., 360. / num_colors):
        hue = i/360.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))

    return colors



def smooth(x,window_len=11,window='hanning'):
    """taken from http://www.scipy.org/Cookbook/SignalSmooth

    smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    import numpy

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    #return y
    return y[(window_len/2-1):-(window_len/2)]
