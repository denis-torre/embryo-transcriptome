#################################################################
#################################################################
############### Embryo Illumina - Python Support ############
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Sebra and Guccione Laboratories
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import pandas as pd
import numpy as np
import re

##### 2. Custom modules #####
# Pipeline running

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

#######################################################
#######################################################
########## S6. rMATS
#######################################################
#######################################################

#############################################
########## 3. Summary
#############################################

def createRmatsSummary(rmats_file, filter_events=False):

    # Read
    rmats_dataframe = pd.read_table(rmats_file)

    # Fix chr
    rmats_dataframe['chr'] = [x.replace('chr', '') for x in rmats_dataframe['chr']]

    # Event ID format
    if 'SE' in rmats_file:
        rmats_dataframe['exonStart_1base'] = rmats_dataframe['exonStart_0base']+1
        rmats_dataframe['downstreamES_1base'] = rmats_dataframe['downstreamES']+1
        rmats_dataframe['Event_id'] = ['{GeneID};SE:{chr}:{upstreamEE}-{exonStart_1base}:{exonEnd}-{downstreamES_1base}:{strand}'.format(**rowData) for index, rowData in rmats_dataframe.iterrows()]
        rmats_dataframe['event_type'] = 'SE'
    elif 'A5SS' in rmats_file:
        rmats_dataframe['flankingES_1base'] = rmats_dataframe['flankingES']+1
        rmats_dataframe['longExonStart_1base'] = rmats_dataframe['longExonStart_0base']+1
        rmats_dataframe['shortES_1base'] = rmats_dataframe['shortES']+1
        rmats_dataframe['Event_id'] = ['{GeneID};A5:{chr}:{longExonEnd}-{flankingES_1base}:{shortEE}-{flankingES_1base}:{strand}'.format(**rowData) if rowData['strand'] == '+' else '{GeneID};A5:{chr}:{flankingEE}-{longExonStart_1base}:{flankingEE}-{shortES_1base}:{strand}'.format(**rowData) for index, rowData in rmats_dataframe.iterrows()] 
        rmats_dataframe['event_type'] = 'A5'
    elif 'A3SS' in rmats_file:
        rmats_dataframe['flankingES_1base'] = rmats_dataframe['flankingES']+1
        rmats_dataframe['longExonStart_1base'] = rmats_dataframe['longExonStart_0base']+1
        rmats_dataframe['shortES_1base'] = rmats_dataframe['shortES']+1
        rmats_dataframe['Event_id'] = ['{GeneID};A3:{chr}:{flankingEE}-{longExonStart_1base}:{flankingEE}-{shortES_1base}:{strand}'.format(**rowData) if rowData['strand'] == '+' else '{GeneID};A3:{chr}:{longExonEnd}-{flankingES_1base}:{shortEE}-{flankingES_1base}:{strand}'.format(**rowData) for index, rowData in rmats_dataframe.iterrows()] 
        rmats_dataframe['event_type'] = 'A3'
    elif 'RI' in rmats_file:
        rmats_dataframe['riExonStart_1base'] = rmats_dataframe['riExonStart_0base']+1
        rmats_dataframe['downstreamES_1base'] = rmats_dataframe['downstreamES']+1
        rmats_dataframe['Event_id'] = ['{GeneID};RI:{chr}:{riExonStart_1base}-{upstreamEE}:{downstreamES_1base}-{downstreamES_1base}:{riExonEnd}:{strand}'.format(**rowData) for index, rowData in rmats_dataframe.iterrows()]
        rmats_dataframe['event_type'] = 'RI'
    elif 'MXE' in rmats_file:
        rmats_dataframe['1stExonStart_1base'] = rmats_dataframe['1stExonStart_0base']+1
        rmats_dataframe['downstreamES_1base'] = rmats_dataframe['downstreamES']+1
        rmats_dataframe['2ndExonStart_1base'] = rmats_dataframe['2ndExonStart_0base']+1
        rmats_dataframe['Event_id'] = ['{GeneID};MX:{chr}:{upstreamEE}-{1stExonStart_1base}:{1stExonEnd}-{downstreamES_1base}:{upstreamEE}-{2ndExonStart_1base}:{2ndExonEnd}-{downstreamES_1base}:{strand}'.format(**rowData) for index, rowData in rmats_dataframe.iterrows()]
        rmats_dataframe['event_type'] = 'MX'
        
    # Get average values
    for x in ['IJC', 'SJC']:
        for y in ['1', '2']:
            col = '{x}_SAMPLE_{y}'.format(**locals())
            rmats_dataframe['mean_{x}_{y}'.format(**locals())] = [np.mean([int(y) for y in x.split(',')]) for x in rmats_dataframe[col]]
        rmats_dataframe['mean_{x}'.format(**locals())] = rmats_dataframe[['mean_{x}_1'.format(**locals()), 'mean_{x}_2'.format(**locals())]].mean(axis=1)

    # select columns
    rmats_dataframe = rmats_dataframe[['Event_id', 'event_type', 'mean_IJC', 'mean_SJC', 'IncLevel1', 'IncLevel2', 'IncLevelDifference', 'PValue', 'FDR']]

    # Add coordinates
    rmats_dataframe['igv_coordinates'] = ['chr{0}:{1}-{2}'.format(*[re.match('.*?-(.*?)-(.*?)-.*-(.*?)-.', x.replace(':', '-')).group(i) for i in [1,2,3]]) for x in rmats_dataframe['Event_id']]

    # Significance
    if filter_events:
        rmats_dataframe['significant'] = [rowData['FDR'] < 0.05 and abs(rowData['IncLevelDifference']) > 0.2 and rowData['mean_IJC'] > 5 and rowData['mean_SJC'] > 5 for index, rowData in rmats_dataframe.iterrows()]

    # Return
    return rmats_dataframe

#######################################################
#######################################################
########## S8. SUPPA
#######################################################
#######################################################

#############################################
########## 5. Event Summary
#############################################

def createSuppaSummary(basename, filter_events=False, event_type='AS'):

    # Read p-value
    pvalue_file = basename+'.dpsi'
    pvalue_dataframe = pd.read_table(pvalue_file).rename_axis('Event_id')
    pvalue_dataframe = pvalue_dataframe.rename(columns={x: x.split('_')[-1].replace('-', '') for x in pvalue_dataframe.columns}).reset_index()

    # Read PSI
    psi_file = basename+'.psivec'
    psi_dataframe = pd.read_table(psi_file).rename_axis('Event_id').reset_index()

    # Get mean PSI
    mean_psi_dataframe = pd.melt(psi_dataframe, id_vars='Event_id')
    mean_psi_dataframe['group'] = ['PSI_'+x.rsplit('_', 1)[0] for x in mean_psi_dataframe['variable']]
    mean_psi_dataframe = mean_psi_dataframe.pivot_table(index='Event_id', columns='group', values='value', aggfunc=np.mean).reset_index()

    # Get logTPM
    tpm_file = basename+'_avglogtpm.tab'
    tpm_dataframe = pd.read_table(tpm_file, header=None, names=['Event_id', 'avg_logTPM'])

    # Merge
    result_dataframe = pvalue_dataframe.merge(tpm_dataframe, on='Event_id', how='left')#.merge(mean_psi_dataframe, on='Event_id', how='left')

    # Add AS-specific events, if not isoform results
    if event_type == 'AS':

        # Add event type
        result_dataframe['event_type'] = basename.split('-')[-1]

        # Add coordinates
        result_dataframe['igv_coordinates'] = ['chr{0}:{1}-{2}'.format(*[re.match('.*?-(.*?)-(.*?)-.*-(.*?)-.', x.replace(':', '-')).group(i) for i in [1,2,3]]) for x in result_dataframe['Event_id']]

    # Significance
    if filter_events:
        result_dataframe['significant'] = [rowData['FDR'] < 0.05 and abs(rowData['IncLevelDifference']) > 0.2 and rowData['mean_IJC'] > 5 and rowData['mean_SJC'] > 5 for index, rowData in result_dataframe.iterrows()]

    # Return
    return result_dataframe



#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################