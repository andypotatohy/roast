#!/usr/bin/python

# given input train log file, extract loss function on validation and training. Feed into python script for plotting

import sys, getopt, os
import matplotlib.pyplot as plt
from subprocess import PIPE, Popen
import numpy as np
import pandas as pd

def makeFiles(argv):
    file = str(argv[0])
    bins = 20
    if len(argv) > 1:
    	bins = int(argv[1])
    bashCommand_getDSC = "grep 'DCS' " +  file + " | awk '{print $2}' "
    p = Popen(bashCommand_getDSC, stdout=PIPE, shell=True)
    output = p.communicate()
    DSC = output[0].split()
    DSC = DSC[:-1]	
    for i in range(0, len(DSC)):
    	DSC[i] = float(DSC[i])
    df = pd.DataFrame(DSC)
    df.plot.hist(bins=bins, edgecolor='black')
    plt.show()
    
if __name__ == "__main__":
	makeFiles(sys.argv[1:])
    
    
    
