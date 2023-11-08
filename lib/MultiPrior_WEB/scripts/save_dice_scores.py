#!/usr/bin/python

# given input train log file, extract loss function on validation and training. Feed into python script for plotting

import sys, getopt, os
import matplotlib.pyplot as plt
from subprocess import PIPE, Popen
import numpy as np



def makeFiles(argv):
    file = ''
    w = 1
    save = False
    try:
		opts, args = getopt.getopt(argv,"hf:",["file="])
    except getopt.GetoptError:
		print('plotLossFunctionKeras.py -f <train log full file address: /home/...> -m <moving average>')
		sys.exit(2)
    for opt, arg in opts:
		if opt == '-h':
			print('plotLossFunctionKeras.py -f <train log full file address : /home/...> -m <moving average>')
			sys.exit()
		elif opt in ("-f","--file"):
			myFile = str(arg)

    bashCommand_getDSC = "grep 'Overall' " + myFile +  "  | awk '{print $3}'  | sed 's/]//' "

    OUTPUT_PATH = '/'.join(myFile.split('/')[:-1])

    p = Popen(bashCommand_getDSC, stdout=PIPE, shell=True)
    output = p.communicate()
    DSC = output[0].split()

    for i in range(0, len(DSC)):
    	DSC[i] = str(DSC[i])

    f = open(OUTPUT_PATH + '/Dice_scores.txt', 'w+')
    for dsc in DSC:
	f.write(dsc)
	f.write('\n')
    f.close()

    
if __name__ == "__main__":
	makeFiles(sys.argv[1:])
    
    
    
