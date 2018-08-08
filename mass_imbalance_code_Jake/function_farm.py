
# coding: utf-8

# In[ ]:


import cobra
import os
from os.path import join
import pandas as pd
print("This statement is printed whenever this file is imported\n")
print('Imported modules below are necessary for proper functioning\n')
print("To get names of functions enter:\ndir(function_farm)")
print('\nTo get help with available functionality of the functions in this list enter:\n\nhelp(FUNCTION_NAME) ')
print('\nTo import all functions enter:\n\nfrom function_farm import * ')
print('_______________________________________________________________________')
import inspect 
print('\nImported inspect module')
import pydoc
print('\nImported pydoc module')
import copy
print('\nImported copy module')
import csv
print('\nImport csv module')

from cobra.flux_analysis import sample
from cobra.medium import minimal_medium
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from cobra.util.util import format_long_string
from cobra.core import get_solution
print('Imported many flux analysis modules and gene deletion modules')

from cobra import Model, Reaction, Metabolite
print('\nImport Model, Reaction, Metabolite modules\n')

try:
    print('Importing Workbook from openpyxl (necessary to complete ) and all pandas functionality\n')
    from openpyxl import Workbook
    from pandas import *
except:
    print('Something went wrong with importing pandas or openpyxl\nIn cmd prompt try this:\n\npip install openpyxl')
    raise


# In[1]:


def checkFreeMass(raw_model, cytosol='c'):   
    '''
            Checks for free mass

            Made by Matt Jenior

            Enter model as first argument and 'c' as the second argument
    '''
    model = copy.deepcopy(raw_model) # changed the with statement -- JJ 7-30-2018
    for index in model.boundary:
        model.reactions.get_by_id(index.id).lower_bound = 0
    demand_metabolites = [x.reactants[0].id for x in model.demands if len(x.reactants) > 0] + [x.products[0].id for x in model.demands if len(x.products) > 0]
    free = []
    for index in model.metabolites: 
        if index.id in demand_metabolites:
            continue
        elif index.compartment != cytosol:
            continue
        else:
            demand = model.add_boundary(index, type='demand')
            model.objective = demand
            obj_val = model.slim_optimize(error_value=0.)
            if obj_val > 1e-8:
                free.append([index.id, obj_val])
            model.remove_reactions([demand])
    if len(free) > 0:
        print(str(len(free)) + ' metabolites are generated for free')
    return(free)


# In[2]:


def changeMedia_PA_LJD(model, mediaCondition, limEx = None):
    '''
        Function Purpose:

            Alters the in silico media. Current media offerings are LB, M9 Minimal media, and SCFM.

        Example Input Style:

            changeMedia_PA_LJD(model, 2, ['EX_cpd00001_LPAREN_e_RPAREN_','EX_cpd00009_LPAREN_e_RPAREN_', 'EX_cpd00011_LPAREN_e_RPAREN_'] )

        INPUTS (in order):

            model - COBRA model structure.

            mediaCondition - A scalar either 1, 2 or 3 to indicate the desired 
            in silico media condition. 
            1 - returns in silico LB.
            2 - returns in silico CFM. cystic fibrosis media
            3 - returns in silico glucose minimal media.

            limEx -some exchange reaction list that is used in Minimal media condition 3 only.

        OUTPUT:

            modelout - COBRA model structure with modified exchange reaction bounds

            Matthew Oberhardt, 1-7-2010
            Jennifer Bartell, 3-27-2013
            Anna Blazier, 9-18-2012 
            Laura Dunphy 9-22-2016
            Jacob Jakielaszek 7-12-2018
    '''
    try: 
        modelout = copy.deepcopy(model) #deep copy so that changes in modelout don't impact input model
        for index in modelout.exchanges: #Turn off all exchanges and allow secretion only
            modelout.reactions.get_by_id(index.id).lower_bound = 0          
            modelout.reactions.get_by_id(index.id).upper_bound = 1000   
        modelout.reactions.get_by_id('EX_cpd00007_LPAREN_e_RPAREN_').lower_bound = -20 #Oxygen bounds set to -20 like in Burk Paper
        try:  
            if mediaCondition == 1: #LB Media
                print('LB Media Selected\n')            
                # %Nutrients such as ions that are freely exchanged:  
                openexchanges = [
                'EX_cpd00001_LPAREN_e_RPAREN_',    #H2O 1379 set
                'EX_cpd00009_LPAREN_e_RPAREN_',    #%Phosphate set
                'EX_cpd00011_LPAREN_e_RPAREN_',   #%CO2 set
                # % 'EX_cpd10515_LPAREN_e_RPAREN_'    %Fe2+ % BC
                'EX_cpd00021_LPAREN_e_RPAREN_',    #%Fe2+ % PA - LJD 1358 set
                'EX_cpd00034_LPAREN_e_RPAREN_',    #%Zn2+ 1460 set
                'EX_cpd00048_LPAREN_e_RPAREN_',    #%Sulfate 1445 set
                'EX_cpd00058_LPAREN_e_RPAREN_',    #%Cu2+ 1331 set
                'EX_cpd00205_LPAREN_e_RPAREN_',    #%K+ 1386 set
                'EX_cpd00254_LPAREN_e_RPAREN_',    #%Mg 1422 set
                'EX_cpd00971_LPAREN_e_RPAREN_',    #%Na+ 1426 set
                'EX_cpd01012_LPAREN_e_RPAREN_',    #%Cd2+ 1324 set
                'EX_cpd00067_LPAREN_e_RPAREN_'    #%H+ 1378 set
                ]

                LBexchanges = [
                'EX_cpd00023_LPAREN_e_RPAREN_',    #%L-Glutamate 1396 set
                'EX_cpd00027_LPAREN_e_RPAREN_',    #%D-Glucose 1348 set
                'EX_cpd00033_LPAREN_e_RPAREN_',    #%Glycine 1368 set
                'EX_cpd00035_LPAREN_e_RPAREN_',    #%L-Alanine 1388 set
                'EX_cpd00039_LPAREN_e_RPAREN_',    #%L-Lysine 1403 set
                'EX_cpd00041_LPAREN_e_RPAREN_',    #%L-Aspartate 1393 set
                'EX_cpd00051_LPAREN_e_RPAREN_',    #%L-Arginine 1391 set
                'EX_cpd00054_LPAREN_e_RPAREN_',    #%L-Serine 1410 set
                'EX_cpd00060_LPAREN_e_RPAREN_',    #%L-Methionine 1405 set
                'EX_cpd00065_LPAREN_e_RPAREN_',    #%Tryptophan 1413 set (L-Tryoptophan)
                'EX_cpd00066_LPAREN_e_RPAREN_',    #L-Phenylalanine 1408 set
                'EX_cpd00069_LPAREN_e_RPAREN_',    #%L-Tyrosine 1414 set
                'EX_cpd00084_LPAREN_e_RPAREN_',    #%L-Cysteine 1395 set
                'EX_cpd00107_LPAREN_e_RPAREN_',    #%L-Leucine 1402 set
                'EX_cpd00119_LPAREN_e_RPAREN_',    #%L-Histidine 1398 set
                'EX_cpd00129_LPAREN_e_RPAREN_',    #%L-Proline 1409 set
                'EX_cpd00156_LPAREN_e_RPAREN_',    #%L-Valine 1415 set
                'EX_cpd00161_LPAREN_e_RPAREN_',    #%L-Threonine 1412 set
                'EX_cpd00305_LPAREN_e_RPAREN_',    #%Thiamin 1448 set
                'EX_cpd00322_LPAREN_e_RPAREN_',    #%L-Isoleucine 1400 set *
                'EX_cpd00092_LPAREN_e_RPAREN_',    #%Uracil 1452 set
                'EX_cpd00307_LPAREN_e_RPAREN_',    #%Cytosine 1336 set
                'EX_cpd03091_LPAREN_e_RPAREN_'    #%5'-Deoxyadenosine 1387 set (L-5"-Deoxy..)
                #% Glycerol exchange 861 (10% supplement)
                ]
                # %changes the lower bound of openexchanges to -1000, the lower bound of
                # %LBexchanges to -10 and the upper bound of LBexchanges to 10.  Also
                # %changes the upper bound of the glucose exchange reaction to 0.
                
                for index in openexchanges: # lower bounds set to -1000
                    modelout.reactions.get_by_id(index).lower_bound = -1000
#                     print(index, sep=' ', end='', flush=True)
                for index in LBexchanges:  # LB Media bounds -10 and upper bound 10 
                    modelout.reactions.get_by_id(index).lower_bound = -10      
                    modelout.reactions.get_by_id(index).upper_bound = 10
                modelout.id = 'LB'
                return modelout
            elif mediaCondition == 2:
                print('CFM Media Selected\n')
                # %Nutrients in CFM
                openexchanges = [
                'EX_cpd00001_LPAREN_e_RPAREN_',    #%H2O 872
                'EX_cpd00009_LPAREN_e_RPAREN_',    #%Phosphate 994
                'EX_cpd00011_LPAREN_e_RPAREN_',    #%CO2 774
                'EX_cpd00021_LPAREN_e_RPAREN_', #%Fe2+ LJD 6/13/17
                # % 'EX_cpd10515_LPAREN_e_RPAREN_'    %Fe2+ 813
                'EX_cpd00023_LPAREN_e_RPAREN_',    #%L-Glutamate  855
                'EX_cpd00027_LPAREN_e_RPAREN_',    #%D-Glucose 849
                'EX_cpd00033_LPAREN_e_RPAREN_',    #%Glycine 856
                'EX_cpd00035_LPAREN_e_RPAREN_',    #%L-Alanine 
                'EX_cpd00039_LPAREN_e_RPAREN_',    #%L-Lysine
                'EX_cpd00041_LPAREN_e_RPAREN_',    #%L-Aspartate
                'EX_cpd00048_LPAREN_e_RPAREN_',    #%Sulfate
                'EX_cpd00051_LPAREN_e_RPAREN_',    #%L-Arginine
                'EX_cpd00054_LPAREN_e_RPAREN_',    #%L-Serine
                'EX_cpd00060_LPAREN_e_RPAREN_',    #%L-Methionine
                'EX_cpd00064_LPAREN_e_RPAREN_',    #%Ornithine
                'EX_cpd00065_LPAREN_e_RPAREN_',    #%Tryptophan
                'EX_cpd00066_LPAREN_e_RPAREN_',    #%L-Phenylalanine
                'EX_cpd00067_LPAREN_e_RPAREN_',    #%H+
                'EX_cpd00069_LPAREN_e_RPAREN_',    #%L-Tyrosine
                'EX_cpd00084_LPAREN_e_RPAREN_',    #%L-Cysteine
                'EX_cpd00107_LPAREN_e_RPAREN_',    #%L-Leucine
                'EX_cpd00119_LPAREN_e_RPAREN_',    #%L-Histidine
                'EX_cpd00129_LPAREN_e_RPAREN_',    #%L-Proline
                'EX_cpd00156_LPAREN_e_RPAREN_',    #%L-Valine
                'EX_cpd00221_LPAREN_e_RPAREN_',    #%D-Lactate
                'EX_cpd00161_LPAREN_e_RPAREN_',    #%L-Threonine
                'EX_cpd00205_LPAREN_e_RPAREN_',    #%K+
                'EX_cpd00209_LPAREN_e_RPAREN_',    #%Nitrate
                'EX_cpd00254_LPAREN_e_RPAREN_',    #%Mg
                'EX_cpd00322_LPAREN_e_RPAREN_',    #%L-Isoleucine
                'EX_cpd00971_LPAREN_e_RPAREN_',    #%Na+
                'EX_cpd00013_LPAREN_e_RPAREN_'    #%NH3
                ]
                #UB CHECK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<JJ 2018, What is this for? DNA?
                ub_check = [
                'EX_cpd00305_LPAREN_e_RPAREN_'    #%Thiamin
                'EX_cpd00092_LPAREN_e_RPAREN_',    #%Uracil
                'EX_cpd00307_LPAREN_e_RPAREN_',    #%Cytosine
                'EX_cpd03091_LPAREN_e_RPAREN_'    #%5'-Deoxyadenosine
                ]   
                #%changes the lower bound for all available nutrients to -10.
                for index in openexchanges:     # lower bounds set to -10
                    modelout.reactions.get_by_id(index).lower_bound = -10
                modelout.id = 'CFM'
                return modelout
            elif mediaCondition == 3:
                print('M9 Minimal Media Selected\n')
                # Minimal M9 Media 
                #%Nutrients such as ions that are freely exchanged:    
                openexchanges = [     
                 'EX_cpd00001_LPAREN_e_RPAREN_', #%H2O
                 'EX_cpd00009_LPAREN_e_RPAREN_', #%Phosphate
                 'EX_cpd00011_LPAREN_e_RPAREN_', #%CO2
                #%  'EX_cpd10515_LPAREN_e_RPAREN_' %Fe2+
                 'EX_cpd00021_LPAREN_e_RPAREN_', #%Fe2+ LJD 5/26/17
#                  'EX_cpd00063_LPAREN_e_RPAREN_', #  Ca 2+ missing from model 7-13-18 JJ
#                  'EX_cpd00131_LPAREN_e_RPAREN_', #  Ca 2+ missing from model 7-13-18 JJ
#                  'EX_cpd00063_LPAREN_e_RPAREN_', #  Ca 2+ missing from model 7-13-18 JJ                    
                 'EX_cpd00030_LPAREN_e_RPAREN_', #%Mn2+
                 'EX_cpd00034_LPAREN_e_RPAREN_', #%Zn2+
                 'EX_cpd00048_LPAREN_e_RPAREN_', #%Sulfate
                 'EX_cpd00058_LPAREN_e_RPAREN_', #%Cu2+
                 'EX_cpd00067_LPAREN_e_RPAREN_', #%H+
                 'EX_cpd00149_LPAREN_e_RPAREN_', #%Co2+
                 'EX_cpd00205_LPAREN_e_RPAREN_', #%K+
                 'EX_cpd00254_LPAREN_e_RPAREN_', #%Mg
                #%  'EX_cJB00102_LPAREN_e_RPAREN_' #%Nitrogen
                 'EX_cpd00528_LPAREN_e_RPAREN_', #%nitrogen LJD 5/26/17
                 'EX_cpd00971_LPAREN_e_RPAREN_', #%Na+
                 'EX_cpd00013_LPAREN_e_RPAREN_', #%NH3
                 'EX_cpd01012_LPAREN_e_RPAREN_', #%Cd2+
                 'EX_cpd10516_LPAREN_e_RPAREN_', #%fe3
                 'EX_cpd00244_LPAREN_e_RPAREN_' #%Ni2+ ADDED BY PHIL
                ]
                limitedexchanges = copy.deepcopy(limEx)
                # %changes the lower bound of openexchanges to -1000 and the lower bound of
                # %the limitedexchanges to -10.  Also changes the upper bound of
                # %limitedexchanges to 0.
                for index in openexchanges:  # lower bounds set to -1000
                    modelout.reactions.get_by_id(index).lower_bound = -1000
                if isinstance(limitedexchanges, list):
                    try:
                        for index in limitedexchanges:  # M9 Media bounds -10 and upper bound 0 
                            modelout.reactions.get_by_id(index).lower_bound = -10
                        modelout.id = 'M9'
                        return modelout
                    except:
                        print('You entered a list, but one or more of the items in the list was not a reaction')
                        raise
                else:
                    print("You selected three (minimal media), but didn't input a limEx change correctly or at all")
                    raise
        except: 
            print('Something is wrong, check your inputs')
            print('Input 1, 2, or 3 in the function, could be related to other input issues') 
            riase
    except: 
        print("Something is terribly wrong with your input")
        print("enter:\nhelp(changeMedia_PA_LJD)\nfor additional information on how to utilize function")
        raise


# In[3]:


def minimal_media_conditions(model, limEx): 
    '''
        Function Purpose:

            Alters the in silico media to M9 ONLY and takes list of carbon sources rather than one carbon source. 

        Example Input Style:

            changeMedia_PA_LJD(model,['EX_cpd00001_LPAREN_e_RPAREN_','EX_cpd00009_LPAREN_e_RPAREN_', 'EX_cpd00011_LPAREN_e_RPAREN_'] )

        INPUTS (in order):

            model - COBRA model structure.

            limEx - Carbon exchange LIST.

        OUTPUT:

            modelout - COBRA model structures in a list with singly modified carbon source for each reaction in list

        Psuedo Example Input:

            minimal_media_conditions(model, [Glucose, Succinate, lysine])

        Psuedo Example Output:

            modelout = ....
            [ Model(glucose carbon source), Model(Succinate carbon source), Model(Lysine carbon source) ]


            Matthew Oberhardt, 1-7-2010
            Jennifer Bartell, 3-27-2013
            Anna Blazier, 9-18-2012 
            Laura Dunphy 9-22-2016
            Jacob Jakielaszek 7-12-2018

    '''
    try: 
        modelout = copy.deepcopy(model) #deep copy so that changes in modelout don't impact input model
        modelout_list = []              #returned list of models with various carbon sources turned on. 
        for index in modelout.exchanges:                                   #Turn off all exchanges and allow secretion only
            modelout.reactions.get_by_id(index.id).lower_bound = 0          
            modelout.reactions.get_by_id(index.id).upper_bound = 1000
        modelout.reactions.get_by_id('EX_cpd00007_LPAREN_e_RPAREN_').lower_bound = -20        #Oxygen bounds set to -20 like in Burk Paper
        try:
            print('M9 Minimal Media Selected\n')
            # Minimal M9 Media 
            #%Nutrients such as ions that are freely exchanged:    
            openexchanges = [     
             'EX_cpd00001_LPAREN_e_RPAREN_', #%H2O
             'EX_cpd00009_LPAREN_e_RPAREN_', #%Phosphate
             'EX_cpd00011_LPAREN_e_RPAREN_', #%CO2
            #%  'EX_cpd10515_LPAREN_e_RPAREN_' %Fe2+
             'EX_cpd00021_LPAREN_e_RPAREN_', #%Fe2+ LJD 5/26/17
#                  'EX_cpd00063_LPAREN_e_RPAREN_', #  Ca 2+ missing from model 7-13-18 JJ
#                  'EX_cpd00131_LPAREN_e_RPAREN_', #  Ca 2+ missing from model 7-13-18 JJ
#                  'EX_cpd00063_LPAREN_e_RPAREN_', #  Ca 2+ missing from model 7-13-18 JJ                    
             'EX_cpd00030_LPAREN_e_RPAREN_', #%Mn2+
             'EX_cpd00034_LPAREN_e_RPAREN_', #%Zn2+
             'EX_cpd00048_LPAREN_e_RPAREN_', #%Sulfate
             'EX_cpd00058_LPAREN_e_RPAREN_', #%Cu2+
             'EX_cpd00067_LPAREN_e_RPAREN_', #%H+
             'EX_cpd00149_LPAREN_e_RPAREN_', #%Co2+
             'EX_cpd00205_LPAREN_e_RPAREN_', #%K+
             'EX_cpd00254_LPAREN_e_RPAREN_', #%Mg
            #%  'EX_cJB00102_LPAREN_e_RPAREN_' #%Nitrogen
             'EX_cpd00528_LPAREN_e_RPAREN_', #%nitrogen LJD 5/26/17
             'EX_cpd00971_LPAREN_e_RPAREN_', #%Na+
             'EX_cpd00013_LPAREN_e_RPAREN_', #%NH3
             'EX_cpd01012_LPAREN_e_RPAREN_', #%Cd2+
             'EX_cpd10516_LPAREN_e_RPAREN_', #%fe3
             'EX_cpd00244_LPAREN_e_RPAREN_' #%Ni2+ ADDED BY PHIL
            ]
            limitedexchanges = limEx;
            # %changes the lower bound of openexchanges to -1000 and the lower bound of
            # %the limitedexchanges to -10.  Also changes the upper bound of
            # %limitedexchanges to 0.
            for index in openexchanges:                                   # lower bounds set to -1000
                modelout.reactions.get_by_id(index).lower_bound = -1000      
            if isinstance(limitedexchanges, list):
                try:
                    print('Function Minimal Media Conditions: M9 Media\nCreating Model(s)...')
                    for index in limitedexchanges:                                     # M9 Media bounds -10 and upper bound 0 
                        temp = copy.deepcopy(modelout) #deep copy so that modelout is not affected
#                         print(index, sep=' ', end='', flush=True)
                        temp.reactions.get_by_id(index).lower_bound = -10     #limited uptake of the carbon source
                        temp.id =  temp.reactions.get_by_id(index).name
                        temp.optimize()
                        modelout_list.append(temp)
                        del temp
                    return modelout_list
                except:
                    print('You entered a list, but one or more of the items in the list was not a reaction')
                    raise
            else:
                print("You selected three (minimal media), but didn't input a limEx change correctly or at all")
                raise
        except: 
            print('Something is wrong, check your inputs')
            print('Input 1, 2, or 3 in the function, could be related to other input issues')    
            raise
    except: 
        print("Something is terribly wrong with your input")
        print("enter:\nhelp(changeMedia_PA_LJD)\nfor additional information on how to utilize function")


# In[4]:


def minimal_media_varied_carbon(model):
    '''
        INPUT:

            Enter a valid model that has been loaded. The model must have the id's found in this function. 

        Output:

            Two excel files with growth ratios for every genes on carbon source.
            Returns a list with the models and carbon sources and their growths. 

            Please implement flux and genes as well...

        FROM MATLAB: 
            % LJD, 12/12/17
            % Goal: To perform gene essentiality simulations on the PA14 model on
            % minimal media with carbon sources from Biolog Phenotypic Microarray Plates PM1 and PM2a. 
            % Output: CSV file (titled '')of growth rate ratios for every gene on each
            % carbon source. This data is equivalent to S4 Data. 

    '''
    print("Assembling known carbon sources into lists...\n")

    # Carbon Source Compound IDs (model formality)
    limEX = [
        'EX_cpd00029_LPAREN_e_RPAREN_', # Acetic Acid
        'EX_cpd00137_LPAREN_e_RPAREN_', # Citric Acid
        'EX_cpd00080_LPAREN_e_RPAREN_', # D,L-alpha-Glycerol Phosphate
        'EX_cpd00117_LPAREN_e_RPAREN_', # D-Alanine
        'EX_cpd00082_LPAREN_e_RPAREN_', # D-Fructose
        'EX_cpd00182_LPAREN_e_RPAREN_', # Adenosine
        'EX_cpd00106_LPAREN_e_RPAREN_', # Fumaric Acid
        'EX_cpd00051_LPAREN_e_RPAREN_', # L-Arginine
        'EX_cpd00132_LPAREN_e_RPAREN_', # L-Asparagine
        'EX_cpd00041_LPAREN_e_RPAREN_', # L-Aspartic Acid
        'EX_cpd00023_LPAREN_e_RPAREN_', # L-Glutamic Acid
        'EX_cpd00053_LPAREN_e_RPAREN_', # L-Glutamine
        'EX_cpd00119_LPAREN_e_RPAREN_', # L-Histidine
        'EX_cpd00100_LPAREN_e_RPAREN_', # Glycerol
        'EX_cpd00033_LPAREN_e_RPAREN_', # Glycine
        'EX_cpd00380_LPAREN_e_RPAREN_', # Itaconic Acid
        'EX_cpd00064_LPAREN_e_RPAREN_', # L-Ornithine
        'EX_cpd00066_LPAREN_e_RPAREN_', # L-Phenylalanine
        'EX_cpd00129_LPAREN_e_RPAREN_', # L-Proline
        'EX_cpd00054_LPAREN_e_RPAREN_', # L-Serine
        'EX_cpd00308_LPAREN_e_RPAREN_', # Malonic Acid
        'EX_cpd00035_LPAREN_e_RPAREN_', # L-Alanine
        'EX_cpd00027_LPAREN_e_RPAREN_', # alpha-D-Glucose
        'EX_cpd00118_LPAREN_e_RPAREN_', # Putrescine
        'EX_cpd00020_LPAREN_e_RPAREN_', # Pyruvic Acid
        'EX_cpd00036_LPAREN_e_RPAREN_', # Succinic Acid
        'EX_cpd00024_LPAREN_e_RPAREN_', # alpha-Keto-Glutaric Acid
        'EX_cpd00281_LPAREN_e_RPAREN_', # gamma-Amino Butyric Acid
        'EX_cpd00322_LPAREN_e_RPAREN_', # L-Isoleucine
        'EX_cpd00130_LPAREN_e_RPAREN_', # L-Malic Acid
        'EX_cpd00159_LPAREN_e_RPAREN_', # L-Lactic Acid
        'EX_cpd00107_LPAREN_e_RPAREN_', # L-Leucine
        'EX_cpd00314_LPAREN_e_RPAREN_', # D-Mannitol
        'EX_cpd00162_LPAREN_e_RPAREN_', # 2-Aminoethanol
        'EX_cpd00072_LPAREN_e_RPAREN_', # D-Fructose-6-Phosphate
        'EX_cpd00280_LPAREN_e_RPAREN_', # D-Galacturonic Acid
        'EX_cpd00089_LPAREN_e_RPAREN_', # D-Glucose-1-Phosphate
        'EX_cpd00079_LPAREN_e_RPAREN_', # D-Glucose-6-Phosphate
        'EX_cpd00609_LPAREN_e_RPAREN_', # D-Saccharic Acid
        'EX_cpd00164_LPAREN_e_RPAREN_', # D-Glucuronic Acid
        'EX_cpd00588_LPAREN_e_RPAREN_', # D-Sorbitol
        'EX_cpd00154_LPAREN_e_RPAREN_', # D-Xylose
        'EX_cpd00139_LPAREN_e_RPAREN_', # Glycolic Acid
        'EX_cpd11589_LPAREN_e_RPAREN_', # Glycyl-L-Aspartic Acid
        'EX_cpd00060_LPAREN_e_RPAREN_', # L-Methionine
        'EX_cpd01242_LPAREN_e_RPAREN_', # 2-Deoxy-D-Ribose
        'EX_cpd00184_LPAREN_e_RPAREN_', # Thymidine
        'EX_cpd00156_LPAREN_e_RPAREN_', # L-Valine
        'EX_cpd00179_LPAREN_e_RPAREN_', # Maltose
        'EX_cpd01262_LPAREN_e_RPAREN_', # Maltotriose
        'EX_cpd00121_LPAREN_e_RPAREN_', # m-Inositol
        'EX_cpd00652_LPAREN_e_RPAREN_', # Mucic Acid
        'EX_cpd00224_LPAREN_e_RPAREN_', # L-Arabinose
        'EX_cpd00142_LPAREN_e_RPAREN_', # Acetoacetic Acid
        'EX_cpd00386_LPAREN_e_RPAREN_', # D-Malic Acid
        'EX_cpd00105_LPAREN_e_RPAREN_', # D-Ribose
        'EX_cpd00550_LPAREN_e_RPAREN_', # D-Serine
        'EX_cpd00246_LPAREN_e_RPAREN_', # Inosine
        'EX_cpd00161_LPAREN_e_RPAREN_', # L-Threonine
        'EX_cpd11592_LPAREN_e_RPAREN_', # Glycyl-L-Glutamic Acid
        'EX_cpd11588_LPAREN_e_RPAREN_', # Glycyl-L-Proline 
        'EX_cpd00851_LPAREN_e_RPAREN_', # Trans-3-Hydroxy-L- Proline
        'EX_cpd00211_LPAREN_e_RPAREN_', # Butyric Acid
        'EX_cpd11585_LPAREN_e_RPAREN_', # L-Alanyl-Glycine
        'EX_cpd00039_LPAREN_e_RPAREN_', # L-Lysine
        'EX_cpd00489_LPAREN_e_RPAREN_', # p-Hydroxy Phenyl Acetic Acid
        'EX_cpd00141_LPAREN_e_RPAREN_', # Propionic Acid
        'EX_cpd00249_LPAREN_e_RPAREN_', # Uridine
        'EX_cpd00797_LPAREN_e_RPAREN_', # beta-Hydroxy Butyric Acid
        'rJB00280'                      # D-Gluconic Acid
    ]  
    
    # Carbon Source Compound names (Consistent with Phenotypic Microarrays)
    limEX_Names = [
        'Acetic Acid',
        'Citric Acid',
        'D,L-alpha-Glycerol Phosphate',
        'D-Alanine',
        'D-Fructose',
        'Adenosine',
        'Fumaric Acid',
        'L-Arginine',
        'L-Asparagine',
        'L-Aspartic Acid',
        'L-Glutamic Acid',
        'L-Glutamine',
        'L-Histidine',
        'Glycerol',
        'Glycine',
        'Itaconic Acid',
        'L-Ornithine',
        'L-Phenylalanine',
        'L-Proline',
        'L-Serine',
        'Malonic Acid',
        'L-Alanine',
        'alpha-D-Glucose',
        'Putrescine',
        'Pyruvic Acid',
        'Succinic Acid',
        'alpha-Keto-Glutaric Acid',
        'gamma-Amino Butyric Acid',
        'L-Isoleucine',
        'L-Malic Acid',
        'L-Lactic Acid',
        'L-Leucine',
        'D-Mannitol',
        '2-Aminoethanol',
        'D-Fructose-6-Phosphate',
        'D-Galacturonic Acid',
        'D-Glucose-1-Phosphate',
        'D-Glucose-6-Phosphate',
        'D-Saccharic Acid',
        'D-Glucuronic Acid',
        'D-Sorbitol',
        'D-Xylose',
        'Glycolic Acid',
        'Glycyl-L-Aspartic Acid',
        'L-Methionine',
        '2-Deoxy-D-Ribose',
        'Thymidine',
        'L-Valine',
        'Maltose',
        'Maltotriose',
        'm-Inositol',
        'Mucic Acid',
        'L-Arabinose',
        'Acetoacetic Acid',
        'D-Malic Acid',
        'D-Ribose',
        'D-Serine',
        'Inosine',
        'L-Threonine',
        'Glycyl-L-Glutamic Acid',
        'Glycyl-L-Proline',
        'Hydroxy-L-Proline',
        'Butyric Acid',
        'L-Alanyl-Glycine',
        'L-Lysine',
        'p-Hydroxy Phenyl Acetic Acid',
        'Propionic Acid',
        'Uridine',
        'beta-Hydroxy Butyric Acid',
        'D-Gluconic Acid'
    ] 
    
    # changeMedia_PA_LJD takes input carbon source limEX and model and outputs list of...
    # models with only 1 carbon source from the list and minimal media conditions
    # also optimizes the data
    print('Optimizing the model with the known carbon sources...\n')
    try: 
        model_MM = minimal_media_conditions(model, limEX)
        print('Success!')
    except: 
        print('minimal_media_varied_carbon function failed\n')
        raise
    # Two list with growth or no growth based off model objective values
    grow_no_grow = []
    obj = []
    for i in range(0 , len(model_MM) ):
        obj.append(model_MM[i].objective.value)
        if model_MM[i].objective.value < 0.0001:
            grow_no_grow.append([model_MM[i].id, 0, 'None', model_MM[i].objective.value])
        else:
            grow_no_grow.append([model_MM[i].id, 1, 'Growth',model_MM[i].objective.value])
    
    print('\nCreating temporary exchange reactions from various compounds in list...\n')
    metList = [ # Carbon source Compound IDs
        'cpd00136_e', # 4-Hydroxy Benzoic Acid
        'cpd00266_e', # D,L-Carnitine
        'cpd00138_e', # D-Mannose
        'cpd00666_e', # D-Tartaric Acid
        'cpd00047_e', # Formic Acid
        'cpd00666_e', # L-Tartaric Acid
        'cpd00794_e', # D-Trehalose
        'cpd00666_e', # m-Tartaric Acid
        'cpd01502_e',  # Citraconic Acid
        'cpd00477_c'  # N_Acetyl Glutamic Acid    #Added by jake was missing from list made _c to work
    ]
    
    metNames = [ # Carbon source Compound names
        '4-Hydroxy Benzoic Acid',
        'D,L-Carnitine',
        'D-Mannose',
        'D-Tartaric Acid',
        'Formic Acid',
        'L-Tartaric Acid',
        'D-Trehalose',
        'm-Tartaric Acid',
        'Citraconic Acid',
        'N_Acetyl Glutamic Acid' #Added by jake was missing from list
    ]

    grRatios_EX = [] # Objective Values of Models
    grow_no_grow_unknown = [] # List of objective values of carbon sources that support growth or didn't
    #Adds the carbon sources in metList to the modelTemporary as exchange reactions
    model_temp = [] # tempoeary models list
    for i in range(0, len(metList)):
        model_temp.append( copy.deepcopy(model) ) # create temporary model 
        reaction = Reaction(metList[i]) # create reaction called the metList name 
        reaction.name = metList[i]  # create reaction called the metList name 
        reaction.subsystem = 'Exchange' #specify as exchange rxn 
        reaction.lower_bound = -1000 # This is the default
        reaction.upper_bound = 1000  # This is the default
        metaboliteExchange = model_temp[i].metabolites.get_by_id(metList[i]) # get the actual metabolite name and properties
        reaction.add_metabolites({
            metaboliteExchange: -1.0 #add the metabolite above to the left side of the exchange rxn eq.
        })
        model_temp[i].add_reactions([reaction]) #add this new exchange reaction to the model
        model_temp[i] = changeMedia_PA_LJD(model_temp[i], 3, [metList[i]] ) #access the single carbon source and add to the function
        model_temp[i].id = metNames[i] # call the temporary model the name of the metabolite (carbon source)
        model_temp[i].optimize() #optimize model
        grRatios_EX.append(model_temp[i].objective.value)  # List of Objective values with each change
        if model_temp[i].objective.value < 0.0001: #Append List
            grow_no_grow_unknown.append([metNames[i], 0, 'None', grRatios_EX[i]])
        else:
            grow_no_grow_unknown.append([metNames[i], 1, 'Growth', grRatios_EX[i]])
    
    try:
        all_models = model_MM + model_temp
    except:
        print('Couldnt add the two models list')
        raise
    try:
        big_growth_list = grow_no_grow + grow_no_grow_unknown
        relevant_results = [ [x[0] for x in big_growth_list], [x[3] for x in big_growth_list]]
    except:
        print('Something wrong with your growth lists')
        raise

    growth_results = input('Please enter a name for the excel file containing Growth/No Growth and Carbon Sources:\n\n')
    growth_results = growth_results + '.xlsx'
    df = DataFrame({'Names': [x[0] for x in big_growth_list], 
                    'Binary Grow': [x[1] for x in big_growth_list],
                    'Text Grow': [x[2] for x in big_growth_list],
                    'Obj_Values': [x[3] for x in big_growth_list]
                   })
    try:
        df.to_excel(growth_results, sheet_name='sheet1', index=False)
    except:
        print("Something wrong with your name that you chose, probably already taken")
        raise
    print('\nPrinted excel files with results. Look for your file name called: ', growth_results)
    
    try:
        return [relevant_results, big_growth_list, all_models]
    except:
        print('Couldnt return the list of growths and the models... weird')
        raise


# In[5]:


def results_genes_flux(model, rxn_name, flux_gene_off_on):
    '''
        results_genes_flux utility:

            essential genes of model outputted nicely AND flux through the reactions being changed outputted

        Input: 

            (1) model: model structure
            (2) rxn_name: string reaction name in list e.g. ['rxn00123','rxn00456',....]
            (3) flux_gene_off_on is list containing two binary values to turn off the reaction functionalities either
                [1,0] [0,0] or [0,1] first element is flux second is gene deletions is on

        Output:

            essential genes in pandas data arrary
    '''
    #two models taken from the original
    try:
        if flux_gene_off_on[0] == 0:
            del_results = []
            print('no gene del results returned')
        else:    
            #Essential Genes Before - Biomass Results Before - Gene Knockouts 
            print('Single Gene Deletion of Model Entered')
            model_genes = copy.deepcopy(model)
            del_results = single_gene_deletion(model_genes)
        if flux_gene_off_on[1] == 0:
            fluxes = []
            print('No fluxes returned')
        else:
            model_flux = copy.deepcopy(model)
            #flux through reaction
            fluxes = []
            loop_reactions = []
            print('Rxn flux loop entered')
            
            for i in range(0, len(rxn_name)):
                loop_reactions.append( model_flux.reactions.get_by_id(rxn_name[i]) )
                if i % 10 == 0:
                    print('Flux reaction number',i)
            fluxes = (flux_variability_analysis(model_flux, reaction_list=loop_reactions, loopless=True) )  #range of fluxes through the reaction
    except:
        print('You must enter a cobra model structure')
        raise

    return [del_results,fluxes]
    #     #the mcomms file virulence linkes genes string names previously acquired into a nice list    
    #     try:
    #         file_name = 'ncomms14631-s3.xlsx'
    #         xl_workbook = pandas.ExcelFile(file_name)  # Load the excel workbook
    #         df = xl_workbook.parse("Virulence-linked genes")  # Parse the sheet into a dataframe
    #         aList = df['names'].tolist()  # Cast the desired column into a python list
    #         [str(i) for i in aList]
    #         print(aList)
    #     except:
    #         print('not necessary but you didnt have the correct file (probably) to output from the prior model')


# In[6]:


def compared(Obj_before, Obj_after):
    '''
    Function: 
    
        compared
    
    Purpose: 
    
        compares objective values in list...
        checks to see if separated by certain threshold 10 % of 1% of Obj_before
        checks to see the total difference between the two obj values (.001) or greater
    
    Input: 
    
        Obj_before and Obj_after. Two list of the same length containing obj values for a given carbon source
        obj_before .g. [[Glucose, 1.12423], [Succinate, 1.232512], [Lysine, .00001]]
        obj_before e.g. [[Glucose, 1.1240], [Succinate, .75000], [Lysine, .00001]]
        first element should be 
        
    Output: 
    
        List containing if the rxn changes had any significant impact and if they were *good or not
        *good defined as a fixing the carbon source 
    '''
    
    try: #just a test
        names = Obj_before[0]
        obj_1 = Obj_before[1]
        obj_2 = Obj_after[1]
    except: 
        print('Your inputs are not list or are way off')
        raise    
    

    #List of the media that gives false positives/negativites

    bad_carbons = [
        'D-Malic Acid',
        'D-Ribose',
        'D-Serine',
        'Inosine',
        'L-Threonine',
        'N-Acetyl-L-Glutamic Acid',
        'Glycyl-L-Glutatmic Acid',
        'Glycyl-L-Proline',
        'Hydroxy-L-Proline',
        'Butric Acid',
        'L-Alanyl-Glycine',
        'L-Lysine',
        'D,L-Carnitine',
        'D-Gluconic Acid',
        'p-Hydroxy Phenyl Acetic Acid',
        'Propionic Acid',
        'Uridine',
        'B-Hydroxy Butryic Acid'
    ]

    bad_carbon_ids = [
        'EX_cpd00386_LPAREN_e_RPAREN_', 
        'EX_cpd00105_LPAREN_e_RPAREN_',
        'EX_cpd00550_LPAREN_e_RPAREN_',
        'EX_cpd00246_LPAREN_e_RPAREN_', 
        'EX_cpd00161_LPAREN_e_RPAREN_',
        'cpd00477_e',
        'EX_cpd11592_LPAREN_e_RPAREN_',
        'EX_cpd11588_LPAREN_e_RPAREN_',
        'EX_cpd00851_LPAREN_e_RPAREN_',
        'EX_cpd00211_LPAREN_e_RPAREN_',
        'EX_cpd11585_LPAREN_e_RPAREN_',
        'EX_cpd00039_LPAREN_e_RPAREN_',
        'cpd00266_e',
        'rJB00280',
        'EX_cpd00489_LPAREN_e_RPAREN_',
        'EX_cpd00141_LPAREN_e_RPAREN_',
        'EX_cpd00249_LPAREN_e_RPAREN_',
        'EX_cpd00797_LPAREN_e_RPAREN_'
    ]

    exp_left_comp_right = []    #growth you tabulted before adapted from figure
    results = [] #Out put [Carbon source, sig. mag. change, sig. percent change,  percent change #, mag. change #,]   ]
    meaningful_changes = []
    
    for i in range(0,len(bad_carbons) ):
        if i < 5:
            exp_left_comp_right.append([0, 1]) #0 is no growth 1 is growth
        else:
            exp_left_comp_right.append([1, 0])
    
    if (len(obj_1) == len(obj_2)):
        for i in range(0, len(names)):
            
            #big picture growth change results
            if obj_1[i] < 0.0001:
                if obj_2[i] < 0.0001:
                    pass
                else:
                    meaningful_changes.append(names[i])
            else:
                if obj_2[i] < 0.0001:
                    meaningful_changes.append(names[i])
                else:
                    pass
            #nitty gritty magnitude of change and percent change

            if abs(obj_2[i] - obj_1[i]) > 0.0001: #is the magnitude of change greater then threshold
                if (abs(obj_2[i]/obj_1[i] - 1) > .001): #is the percent change significant
                    results.append([ names[i], 1, obj_2[i] - obj_1[i], 1, 100 * (obj_2[i]/obj_1[i] - 1) ] )
                else:
                    results.append([ names[i],1, obj_2[i] - obj_1[i], 0, 100 * (obj_2[i]/obj_1[i] - 1) ] )
            else:
                if (abs(obj_2[i]/obj_1[i] - 1) > .001):
                    results.append([ names[i],0, obj_2[i] - obj_1[i], 1, 100 * (obj_2[i]/obj_1[i] - 1) ] )
                else:
                    results.append([ names[i],0, obj_2[i] - obj_1[i], 0, 100 * (obj_2[i]/obj_1[i] - 1) ] )
    else:
        print('List were not the same length or not list')
    if len(meaningful_changes) == 0:
        print('\nNo meaningful changes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

    return [results, meaningful_changes]


# In[7]:


get_ipython().system('jupyter nbconvert --to script function_farm.ipynb')

