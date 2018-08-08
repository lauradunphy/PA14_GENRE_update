
# coding: utf-8

# In[ ]:
# LJD_Functions.py
# Purpose: To contain all cobra functions that I use regularly such that I can import them into whatever main script I am using 
# Functions contained in this script:
    # changeMedia_PA_LJD - changes the media conditions of a cobrapy model

# import packages
from copy import *

# changeMedia_PA_LJD adapted from the Matlab function with the same name
    # Changes the media condition of the model by altering the bounds of the appropriate exchange reactions
    # Inputs:
        # model: a cobrapy model
        # media: the media condition to change the environment to
            # 1 - LB media
            # 2 - SCFM
            # 3 - Minimal media
        # limEX: a list of exchanges ('EX_cpd#####(e)') that are added as limited exchanges to minimal media
    # Outputs:
        # modelOutput: a cobrapy model in the new media condition (with adjusted exchanges)
def changeMedia_PA_LJD(model, media, limEX=[]):
    modelOutput = deepcopy(model)
    
    # LB Components:
    # Metabolites that are in excess in LB
    LB_open_exchanges = ['EX_cpd00001(e)', 
                     'EX_cpd00009(e)',
                     'EX_cpd00011(e)',
                     'EX_cpd00021(e)',
                     'EX_cpd00034(e)',
                     'EX_cpd00048(e)',
                     'EX_cpd00058(e)',
                     'EX_cpd00205(e)',
                     'EX_cpd00254(e)',
                     'EX_cpd00971(e)',
                     'EX_cpd01012(e)',
                     'EX_cpd00067(e)']
    
    # Metabolites that are more limited but available in LB
    LB_limited_exchanges = ['EX_cpd00023(e)',
                            'EX_cpd00027(e)',
                            'EX_cpd00033(e)',
                            'EX_cpd00035(e)',
                            'EX_cpd00039(e)',
                            'EX_cpd00041(e)',
                            'EX_cpd00051(e)',    
                            'EX_cpd00054(e)',    
                            'EX_cpd00060(e)',    
                            'EX_cpd00065(e)',    
                            'EX_cpd00066(e)',    
                            'EX_cpd00069(e)',    
                            'EX_cpd00084(e)',    
                            'EX_cpd00107(e)',    
                            'EX_cpd00119(e)',    
                            'EX_cpd00129(e)',    
                            'EX_cpd00156(e)',    
                            'EX_cpd00161(e)',    
                            'EX_cpd00305(e)',   
                            'EX_cpd00322(e)',    
                            'EX_cpd00092(e)',    
                            'EX_cpd00307(e)',
                            'EX_cpd03091(e)']
    
    # SCFM Components:
    SCFM_open_exchanges = ['EX_cpd00001(e)',
                           'EX_cpd00009(e)',
                           'EX_cpd00011(e)',
                           'EX_cpd00021(e)',
                           'EX_cpd00023(e)',
                           'EX_cpd00027(e)',
                           'EX_cpd00033(e)',
                           'EX_cpd00035(e)',    
                           'EX_cpd00039(e)',   
                           'EX_cpd00041(e)',    
                           'EX_cpd00048(e)',    
                           'EX_cpd00051(e)',   
                           'EX_cpd00054(e)',   
                           'EX_cpd00060(e)',    
                           'EX_cpd00064(e)',   
                           'EX_cpd00065(e)',    
                           'EX_cpd00066(e)',    
                           'EX_cpd00067(e)',    
                           'EX_cpd00069(e)',   
                           'EX_cpd00084(e)',    
                           'EX_cpd00107(e)',    
                           'EX_cpd00119(e)',    
                           'EX_cpd00129(e)',    
                           'EX_cpd00156(e)',    
                           'EX_cpd00221(e)',   
                           'EX_cpd00161(e)',   
                           'EX_cpd00205(e)',   
                           'EX_cpd00209(e)',    
                           'EX_cpd00254(e)',   
                           'EX_cpd00322(e)',    
                           'EX_cpd00971(e)',    
                           'EX_cpd00013(e)']
    
    # Minimal Media Components: 
    minimalMedia_open_exchanges = ['EX_cpd00001(e)', 
                                   'EX_cpd00009(e)', 
                                   'EX_cpd00011(e)',
                                   'EX_cpd00021(e)',
                                   'EX_cpd00030(e)',
                                   'EX_cpd00034(e)',
                                   'EX_cpd00048(e)',
                                   'EX_cpd00058(e)',
                                   'EX_cpd00067(e)',
                                   'EX_cpd00149(e)',
                                   'EX_cpd00205(e)',
                                   'EX_cpd00254(e)',
                                   'EX_cpd00528(e)',
                                   'EX_cpd00971(e)', 
                                   'EX_cpd00013(e)', 
                                   'EX_cpd01012(e)', 
                                   'EX_cpd10516(e)',
                                   'EX_cpd00244(e)']
   
    # Set the new media conditions
    for i in range(0,len(modelOutput.exchanges)):
        # Set all EX rxns lb = 0, up = 1000
        if modelOutput.exchanges[i].name[0:2] == 'EX':
            modelOutput.exchanges[i].upper_bound = 1000
            modelOutput.exchanges[i].lower_bound = 0
        # Set aerobic exchange (O2) to 20 mmol/gDW/hr 
        if modelOutput.exchanges[i].id[0:11] == 'EX_cpd00007':
            modelOutput.exchanges[i].lower_bound = -20
    
        # Media == 1 Change bounds to LB
        if media == 1:
            # Set open lb exchanges to -1000
            for j in range(0,len(LB_open_exchanges)):
                if modelOutput.exchanges[i].id[0:11] == LB_open_exchanges[j][0:-3]:
                    modelOutput.exchanges[i].lower_bound = -1000
            # Set limited lb exchanges to -10
            for k in range(0,len(LB_limited_exchanges)):
                if modelOutput.exchanges[i].id[0:11] == LB_limited_exchanges[k][0:-3]:
                    modelOutput.exchanges[i].lower_bound = -10
                    #modelOutput.exchanges[i].upper_bound = 10

        
        # Media == 2 Change bounds to SCFM
        elif media == 2:
            for j in range(0,len(SCFM_open_exchanges)):
                if modelOutput.exchanges[i].id[0:11] == SCFM_open_exchanges[j][0:-3]:
                    modelOutput.exchanges[i].lower_bound = -10 #According to the matlab script...
        
        # Media == 3 Change bounds to minimal + exchanges?
        elif media == 3:
            for j in range(0,len(minimalMedia_open_exchanges)):
                if modelOutput.exchanges[i].id[0:11] == minimalMedia_open_exchanges[j][0:-3]:
                    modelOutput.exchanges[i].lower_bound = -1000
            if len(limEX) > 0:
                for k in range(0,len(limEX)):
                    if modelOutput.exchanges[i].id[0:11] == limEX[k][0:-3]:
                        modelOutput.exchanges[i].lower_bound = -10
        else:
            print('unrecognized media condition. Please enter 1 for LB; 2 for SCFM; 3 for minimal media')
   
    return(modelOutput) 

