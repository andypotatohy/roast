# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 13:22:35 2018

@author: hirsch
"""

if __name__ == "__main__":

    import warnings
    import os
    import sys
    warnings.filterwarnings("ignore")
    sys.path.append(os.getcwd())
    sys.path.append(os.getcwd() + '/scripts')
       
    print(sys.path)

    if len(sys.argv) < 2:
        print('Please include a model configuration file:')
        print('>>> python MAKE_MODEL.py </path_to/model_config.py> \n')
        sys.exit()

    print(sys.argv[1:][0])
    from scripts.lib import make_model
    configFile = sys.argv[1:][0]
    make_model(configFile)

