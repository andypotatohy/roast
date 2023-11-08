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
     
    if len(sys.argv) < 2:
        print('Please include a model configuration file:')
        print('>>> python SEGMENT.py </path_to/segmentation_config.py> \n')
        sys.exit()

    workingDir = os.getcwd()

    from scripts.lib import segment
    configFile = sys.argv[1:][0]#'/home/hirsch/Documents/projects/brainSegmentation/multiScale_TPM/configFile.py'
    segment(configFile, workingDir)


