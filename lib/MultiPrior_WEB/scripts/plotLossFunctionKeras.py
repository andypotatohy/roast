#!/usr/bin/python

# given input train log file, extract loss function on validation and training. Feed into python script for plotting

import sys, getopt, os
import matplotlib.pyplot as plt
from subprocess import PIPE, Popen
import numpy as np



def movingAverageConv(a, window_size=1) :
    if not a : return a
    window = np.ones(int(window_size))
    result = np.convolve(a, window, 'full')[ : len(a)] # Convolve full returns array of shape ( M + N - 1 ).
    slotsWithIncompleteConvolution = min(len(a), window_size-1)
    result[slotsWithIncompleteConvolution:] = result[slotsWithIncompleteConvolution:]/float(window_size)
    if slotsWithIncompleteConvolution > 1 :
        divisorArr = np.asarray(range(1, slotsWithIncompleteConvolution+1, 1), dtype=float)
        result[ : slotsWithIncompleteConvolution] = result[ : slotsWithIncompleteConvolution] / divisorArr
    return result

def makeFiles(argv):
    file = ''
    w = 1
    v = 1
    save = False
    try:
		opts, args = getopt.getopt(argv,"hsf:m:v:",["file=", "movingAverage=", "movingAverageVal="])
    except getopt.GetoptError:
		print('plotLossFunctionKeras.py -f <train log full file address: /home/...> -m <moving average> -v <moving average val> -s <save>')
		sys.exit(2)
    for opt, arg in opts:
		if opt == '-h':
			print('plotLossFunctionKeras.py -f <train log full file address: /home/...> -m <moving average> -v <moving average val> -s <save>')
			sys.exit()
		elif opt in ("-f","--file"):
			file = str(arg)
		elif opt in ("-m","--movingAverage"):
			w = int(arg)   
		elif opt in ("-v","--movingAverageVal"):
			v = int(arg)   
		elif opt in ("-s"):
			save = True   


    bashCommand_getTrain = "grep 'Train cost and metrics' " +  file + " | awk '{print $5}' | grep -o '[0-9].*' | sed 's/,//' "
    bashCommand_getDICE1Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $6}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//' "
    bashCommand_getDICE2Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $7}' | grep -o '[0-9].*' | sed 's/]//'| sed 's/,//' "
    bashCommand_getDICE3Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $8}' | grep -o '[0-9].*' | sed 's/]//'| sed 's/,//' "
    bashCommand_getDICE4Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $9}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//' "
    bashCommand_getDICE5Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $10}' | grep -o '[0-9].*' | sed 's/]//'| sed 's/,//' "
    bashCommand_getDICE6Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $11}' | grep -o '[0-9].*' | sed 's/]//'| sed 's/,//' "
    bashCommand_getDICE7Train = "grep 'Train cost and metrics' " + file +  " | awk '{print $12}' | grep -o '[0-9].*' | sed 's/]//'| sed 's/,//' "

    bashCommand_getVal = "grep 'Validation cost and accuracy' " +  file + " | awk '{print $5}' | grep -o '[0-9].*' | sed 's/,//' "
    bashCommand_getDICE1Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $6}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//'"
    bashCommand_getDICE2Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $7}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//' "
    bashCommand_getDICE3Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $8}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//'"
    bashCommand_getDICE4Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $9}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//' "    
    bashCommand_getDICE5Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $10}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//'"
    bashCommand_getDICE6Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $11}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//' "
    bashCommand_getDICE7Val = "grep 'Validation cost and accuracy' " + file +  " | awk '{print $12}' | grep -o '[0-9].*' | sed 's/]//' | sed 's/,//' "

    bashCommand_getDSC = "grep 'Overall' " + file +  "  | awk '{print $3}'  | sed 's/]//' "


    p = Popen(bashCommand_getTrain, stdout=PIPE, shell=True)
    output = p.communicate()
    train = output[0].split()

    p = Popen(bashCommand_getDICE1Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice1 = output[0].split()
    	
    p = Popen(bashCommand_getDICE2Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice2 = output[0].split()

    p = Popen(bashCommand_getDICE3Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice3 = output[0].split()
    	
    p = Popen(bashCommand_getDICE4Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice4 = output[0].split()

    p = Popen(bashCommand_getDICE5Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice5 = output[0].split()
    	
    p = Popen(bashCommand_getDICE6Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice6 = output[0].split()

    p = Popen(bashCommand_getDICE7Train, stdout=PIPE, shell=True)
    output = p.communicate()
    Tdice7 = output[0].split()


    p = Popen(bashCommand_getVal, stdout=PIPE, shell=True)
    output = p.communicate()
    Val = output[0].split()

    p = Popen(bashCommand_getDICE1Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice1 = output[0].split()

    p = Popen(bashCommand_getDICE2Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice2 = output[0].split()

    p = Popen(bashCommand_getDICE3Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice3 = output[0].split()

    p = Popen(bashCommand_getDICE4Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice4 = output[0].split()    

    p = Popen(bashCommand_getDICE5Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice5 = output[0].split()

    p = Popen(bashCommand_getDICE6Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice6 = output[0].split()

    p = Popen(bashCommand_getDICE7Val, stdout=PIPE, shell=True)
    output = p.communicate()
    Valdice7 = output[0].split()

    p = Popen(bashCommand_getDSC, stdout=PIPE, shell=True)
    output = p.communicate()
    DSC = output[0].split()



    for i in range(0,len(train)-1):
        train[i] = float(train[i])
    train = train[:-1]

    for i in range(0,len(Tdice1)-1):
        Tdice1[i] = float(Tdice1[i])
    Tdice1 = Tdice1[:-1]
    
    for i in range(0, len(Tdice2)-1):
        Tdice2[i] = float(Tdice2[i])
    Tdice2 = Tdice2[:-1]

    for i in range(0,len(Tdice3)-1):
        Tdice3[i] = float(Tdice3[i])
    Tdice3 = Tdice3[:-1]
    
    for i in range(0, len(Tdice4)-1):
        Tdice4[i] = float(Tdice4[i])
    Tdice4 = Tdice4[:-1]

    for i in range(0,len(Tdice5)-1):
        Tdice5[i] = float(Tdice5[i])
    Tdice5 = Tdice5[:-1]
    
    for i in range(0, len(Tdice6)-1):
        Tdice6[i] = float(Tdice6[i])
    Tdice6 = Tdice6[:-1]

    for i in range(0, len(Tdice7)-1):
        Tdice7[i] = float(Tdice7[i])
    Tdice7 = Tdice7[:-1]


    for i in range(0, len(Val)-1):
    	Val[i] = float(Val[i])
    Val = Val[:-1]

    for i in range(0, len(Valdice1)):
    	Valdice1[i] = float(Valdice1[i])
    Valdice1 = Valdice1[:-1]
    for i in range(0, len(Valdice2)):
    	Valdice2[i] = float(Valdice2[i])
    Valdice2 = Valdice2[:-1]
    for i in range(0, len(Valdice3)):
    	Valdice3[i] = float(Valdice3[i])
    Valdice3 = Valdice3[:-1]
    for i in range(0, len(Valdice4)):
    	Valdice4[i] = float(Valdice4[i])
    Valdice4 = Valdice4[:-1]
    for i in range(0, len(Valdice5)):
    	Valdice5[i] = float(Valdice5[i])
    Valdice5 = Valdice5[:-1]
    for i in range(0, len(Valdice6)):
    	Valdice6[i] = float(Valdice6[i])
    Valdice6 = Valdice6[:-1]
    for i in range(0, len(Valdice7)):
    	Valdice7[i] = float(Valdice7[i])
    Valdice7 = Valdice7[:-1]


    for i in range(0, len(DSC)):
    	DSC[i] = float(DSC[i])
    #DSC = DSC[:-1]


    train = movingAverageConv(train, window_size = w)
    Tdice1 = movingAverageConv(Tdice1, window_size = w)
    Tdice2 = movingAverageConv(Tdice2, window_size = w)
    Tdice3 = movingAverageConv(Tdice3, window_size = w)
    Tdice4 = movingAverageConv(Tdice4, window_size = w)
    Tdice5 = movingAverageConv(Tdice5, window_size = w)
    Tdice6 = movingAverageConv(Tdice6, window_size = w)
    Tdice7 = movingAverageConv(Tdice7, window_size = w)

    Val = movingAverageConv(Val, window_size = v)
    Valdice1 = movingAverageConv(Valdice1, window_size = v)
    Valdice2 = movingAverageConv(Valdice2, window_size = v)
    Valdice3 = movingAverageConv(Valdice3, window_size = v)
    Valdice4 = movingAverageConv(Valdice4, window_size = v)
    Valdice5 = movingAverageConv(Valdice5, window_size = v)
    Valdice6 = movingAverageConv(Valdice6, window_size = v)

    #plt.clf()

    
    
    if len(Valdice1) == 0:
            plt.figure(figsize=(4,8))
	    ax1 = plt.subplot(311)
	    plt.plot(range(len(train)),train,'k-')
	    plt.title('Training Data - moving average {}'.format(w),)
	    plt.axis('tight')
	    plt.legend(('Training Loss',))

	    ax3 = plt.subplot(312)
	    plt.plot(range(len(Tdice1)),Tdice1,'k-')
	    plt.plot(range(len(Tdice2)),Tdice2,'k--')
	    plt.plot(range(len(Tdice3)),Tdice3,'g-')
	    plt.plot(range(len(Tdice4)),Tdice4,color='yellow')
	    plt.plot(range(len(Tdice5)),Tdice5,color='orange')
	    plt.plot(range(len(Tdice6)),Tdice6,'r-')
	    plt.plot(range(len(Tdice7)),Tdice7,'b-')
	    plt.legend(('Air','Cav', 'GM','WM','CSF','Bone','Skin'), loc='lower center')


	    ax5 = plt.subplot(313)#, sharey=ax3)
	    #plt.ylim([np.min(DSC)-0.01,np.max(DSC)+0.01])
	    plt.plot(range(2,len(DSC)+2),DSC,'b-o')
	    plt.legend(('Full Segmentation Dice',), loc='lower right')
	    plt.xlabel('Epochs')



    else:
	    plt.figure(figsize=(12,12))
	    ax1 = plt.subplot(321)
	    plt.plot(range(len(train)),train,'k-')
	    plt.title('Training Data - moving average {}'.format(w),)
	    plt.axis('tight')
	    plt.legend(('Training Loss',))

	    ax2 = plt.subplot(322, sharey = ax1)


	    plt.plot(np.linspace(0,len(Val),len(train),endpoint=True), train,'b--')
	    plt.plot(range(0,len(Val)),Val,'k-')
	    plt.legend(('Training Loss','Validation Loss',))
	    plt.title('Validation Data- moving average {}'.format(v),)

	    ax3 = plt.subplot(323)
	    plt.plot(range(len(Tdice1)),Tdice1,'k-')
	    plt.plot(range(len(Tdice2)),Tdice2,'k--')
	    plt.plot(range(len(Tdice3)),Tdice3,'g-')
	    plt.plot(range(len(Tdice4)),Tdice4,color='yellow')
	    plt.plot(range(len(Tdice5)),Tdice5,color='orange')
	    plt.plot(range(len(Tdice6)),Tdice6,'r-')
	    plt.plot(range(len(Tdice7)),Tdice7,'b-')
	    plt.legend(('Air','Cav', 'GM','WM','CSF','Bone','Skin'), loc='lower center')

	    ax4 = plt.subplot(324, sharey = ax3)
	    plt.plot(range(0,len(Valdice1)),Valdice1,'k-')
	    plt.plot(range(0,len(Valdice2)),Valdice2,'k--')
	    plt.plot(range(0,len(Valdice3)),Valdice3,'g-')
	    plt.plot(range(0,len(Valdice4)),Valdice4,color='yellow')
	    plt.plot(range(0,len(Valdice5)),Valdice5,color='orange')
	    plt.plot(range(0,len(Valdice6)),Valdice6,'r-')
	    plt.plot(range(0,len(Valdice7)),Valdice7,'b-')
	    plt.legend(('Air','Cav', 'GM','WM','CSF','Bone','Skin'), loc='lower center')

	    ax5 = plt.subplot(326)#, sharey=ax4)
	    plt.ylim([0.6,1])
	    plt.plot(range(len(DSC)),DSC,'b-o')
	    plt.legend(('Full Segmentation Dice',), loc='lower right')
	    plt.xlabel('Epochs')
    
    if save:
        out = '/'.join(file.split('/')[:-1])
	print('saved in {}'.format(out))
    	plt.savefig(out+'/Training_session.png')
	plt.clf()

    else:
	#print('bye')
        plt.show()
    
    
if __name__ == "__main__":
	makeFiles(sys.argv[1:])
    
    
    
