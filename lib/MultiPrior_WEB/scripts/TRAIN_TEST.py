#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 17:36:27 2018

@author: lukas
"""

if __name__ == "__main__":

    import warnings
    warnings.filterwarnings("ignore")
    import os
    import sys
    sys.path.append(os.getcwd())
    sys.path.append(os.getcwd() + '/scripts')
        
    from scripts.lib import train_test_model

    workingDir = os.getcwd()
    
    if not os.path.exists(workingDir):    
        os.mkdir(workingDir)
    os.chdir(workingDir)
    
    if len(sys.argv) < 2:
        print('Please include a training configuration file:')
        print('>>> python TRAIN_TEST.py </path_to/config.py> \n')
        sys.exit()
        
    if len(sys.argv) == 3:
        print('Using included model')
        
    configFile = sys.argv[1:][0]
    
    train_test_model(configFile, workingDir)


