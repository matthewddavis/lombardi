#!/home/liulab/shin/python/bin/python
#
# ReadParam.py
#
# The class ReadParam reads the parameters needed for peak finding of ChIP sequence data and parses them.
#
# Programmed by Hyunjin Shin (based on JSSONG's ReadInput.py for Nucleosome ChIP-chip data analysis)
#


import sys

class ReadParam:

    """Read the parameters specified in the parameter file """
    
    #################################
    #
    # constructor
    #
    
    def __init__(self, filename):
        self.filename = filename
        self.PARAMS ={}
            
    #################################
    #
    # readParam: read the parameters and parse them
    #
    
    def readParam(self):
        f = open(self.filename)
        tag = ""
        for line in f.xreadlines():
            line = line.strip()
            
            # ignore lines with # predeceding
            
            if line == "" or line[0] == "#": continue
            else:
                l = line.split("=")
                l = [x.strip() for x in l]
                if len(l) > 1:
                    tag = l[0]
                    self.PARAMS[tag] = l[1]
                else:
                    self.PARAMS[tag] = ''
        f.close() 
   
        for t in self.PARAMS:
            self.PARAMS[t] = self.PARAMS[t].strip()

        return self.PARAMS

###################################
#
# usage: show how to use this class
#
 
def usage():   
    print >> sys.stderr, "USAGE: python ReadParam.py PARAMETER_FILE or python ReadParam.py PARAMETER_FILE OUT_FILE"   

###################################
#
# printParam: print the parameter values to standard I/O
#
    
def printParam(param):
    
    if param == {}:
        raise ValueError('Empty parameter set')

    names = param.keys()
    names.sort()
    
    for n in names:
        if param[n] != -1 or param[n] != '':
            print n, "=", param[n]

####################################
#
# main: read the parameters and print them on the standard out
#     
            
if __name__ == "__main__":
    
    if len(sys.argv) < 2: usage()
    else:
        f = ReadParam(sys.argv[1])
        PARAMS = f.readParam()
        printParam(PARAMS)
    