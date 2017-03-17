#!/usr/bin/env python
# coding: utf-8
# Accessibility.py
# 
"""This module calculates accessibility scores given sets of travel time data and land use data.

    Make sure the names of input data files are correct. 
    Also specify the land use activities wanted. 

    There is the option use a cumulative opportunities measure or a gravity-based measure for the accessibility score. 

"""


import pandas as pd
import numpy as np
import os

root = ''
data_dir = 'sample_data'
results_dir = 'results'
"""Directories to find data and to store results"""

taxi_file_2011 = 'ttimes_taxi_2011.csv'
taxi_file_2015 = 'ttimes_taxi_2015.csv'
myciti_file = 'myciti_ttimes.csv'
lu_file = 'lu_data.csv'
walk_file = 'walk_ttimes.csv'
"""Names of data files"""


# set parameters for cumulative opportunities measure
measure = 'cumul_opps' #  'cumul_opps' or 'gravity'
"""Accessibility measure to use"""

tt_max_list = [10,20,30]
"""List of max ttimes (in minutes)"""

modes =['taxi2015','taxi2011','myciti']
"""List of modes for which to get accessiiblity scores"""

beta = 0.1
"""Parameter for gravity measure; needs to be calibrated"""

standardize = False
"""Whether or not to standardize the resulting scores. (Don't standardize if comparing different modes)"""

def load_cost_matrix(path, fname):
    """Load cost matrix from file, calculated from travel time model. 
    Format should be long, with columns named ['start','end','agg_cost']. 
    The columns 'start' and 'end' contain id's of TAZ's, as string. 
    Returns:
        DataFrame: cost matrix with TAZ id's as columns
    """
    df = pd.read_csv(os.path.join(root,path,fname), dtype = {'start':str, 'end':str})
    df = df.sort_values(by=['start','end'])
    return df

def correct_ttime_values(df, update_val=999):
    """Update travel times between TAZs so that when travel is impossible, ttime is very large. 
    This is necessary because in some of the original files, ttime=0 when travel between TAZs is impossible
    Use for when the original matrix is not missing any TAZ pairs. 
    Args: 
        df: matrix with ttimes between TAZs, where 0 indicate
        update_val (num): val to replace 0
    Returns: 
        DataFrame: matrix with updated ttimes between TAZs
    """
    return df.replace(to_replace=0, value=999)

def fill_missing_costs(m,taz_index,missing_val=999):
    """Create a cost matrix with all pairs of TAZs. Use when the given cost matrix is missing pairs. 
    Args: 
        m (DataFrame): cost matrix in long format with TAZ ids as multiindex. Cost col should be named 'agg_cost'
        taz_index: full index of TAZs
        missing_val: value to use for missing values. Default 999.
    Returns: 
        DataFrame: cost matrix in long format with all TAZ pairs as columns, with missing values filled in.
    """

    new_ind=pd.MultiIndex.from_product([taz_index,taz_index], names=['start', 'end'])
    blank = pd.DataFrame(pd.Series(np.nan, index=new_ind))
    blank = blank.reset_index()

    merged=pd.merge(blank, m,on=['start','end'], how='left')

    merged['agg_cost'] = merged['agg_cost'].replace(to_replace=np.nan, value=999)
    return merged[['start','end','agg_cost']]


def add_walk_time(row):
    """Add walking time to start and end of trip"""
    return row['agg_cost']+walk_dict[row['start']]+walk_dict[row['end']]
    
def set_initial_zero(df):
    """set diagonal of ttime matrix to zero"""
    new_values = df.values
    np.fill_diagonal(new_values,0)  # set diagonal to zero. 
    return pd.DataFrame(new_values,index=df.index, columns=df.columns)

def load_activities(path,fname):
    """Load activities data"""
    df = pd.read_csv(os.path.join(root,path,fname), dtype={'taz':str})
    return df.set_index('taz')


def A_i_gravity(i, D, lu, J, df_tt , beta=0.1, imped_fun='neg_exp'):
    """Calculate the accessibility of a single zone i, using a gravity-based measure. 
    
    Args: 
        i (str): id of origin zone
        D (DataFrame): dataframe with values of destination activities, and taz id's as index.
        lu (str): name of land use activities column
        J (list, Series, or Index): index or list of taz id's
        df_tt (DataFrame): matrix of travel times for the given mode and year. 
        beta (float): parameter for impedance factor. Default 0.1. This should be calibrated for the particular context.
        imped_fun (str): functional form of the impedance factor. 'neg_exp' (default), 'power'
        
    Returns: 
        float: Accessibility of zone i. 
    
    Note: I'm putting the impedance factor calculation inside the function, but could also move it outside. 
    """
    to_sum = pd.Series(index=J)

    for j in J:
        c_ij = df_tt.ix[i,j] # cost of travel from i to j. 
        

        # calculate impedance factor. A higher impedance -> costlier travel. 
        if i==j: 
            impedance=1  # impedance of same zone=1. (A location in a zone has full access to all activities in that zone)
        elif imped_fun == 'neg_exp':
            impedance = np.exp(-1*c_ij* beta)   # impedance factor between i and ij. 
        elif imped_fun == 'power':
            impedance = c_ij**(-1*beta)
            
        D_j = D.ix[j,lu] # destination activities (number or weight) in zone j

        to_sum.ix[j] = D_j*impedance   # calculate accessibility for each zone j
    return to_sum.sum() 

def A_i_cumul_opps(i, D, lu, tt_max, df_tt):
    """Calculate the accessibility of a single zone i, using a cumulative opportunities measure. 
    Args: 
        i (str): id of origin zone
        D (DataFrame): dataframe with values of destination activities, and taz id's as index. 
        lu (str): name of land use activities column
        tt_max (float or int): maximum travel time, i.e., travel time cutoff
        df_tt (DataFrame): matrix of travel times for the given mode and year. 
    Returns: 
        float: Accessibility of zone i. 
    """
    J = df_tt.ix[i][df_tt.ix[i]<=tt_max].index  # find set of zones J reachable within max ttime. 
    return D.ix[J].sum()[lu]   # sum of activities available in J


def calculate_all_A(taz_ind, measure='gravity', **kwargs):
    """Calculate A for all taz's. 
    
    Args: 
        taz_ind(list, Series, or Index): index or list of taz id's
        measure (str): 'gravity' (default) or 'cumul_opps'
        **kwargs: keyword arguments needed to calculate A_i
    Returns: 
        Series: Accessibility values for all taz, with taz id's as index. 
    """

    result = pd.Series(index = taz_ind, name='A')
    for i in taz_ind: 
        if measure == 'gravity':
            A_i = A_i_gravity(i=i, D=kwargs['D'], lu=kwargs['lu'], J=taz_ind, beta=kwargs['beta'], df_tt=kwargs['df_tt'])
        elif measure=='cumul_opps':
            A_i = A_i_cumul_opps(i=i, D=kwargs['D'], lu=kwargs['lu'], tt_max=kwargs['tt_max'], df_tt=kwargs['df_tt'])
        result.ix[i] = A_i
    return result

def standardize_A(A_series):
    """Compute standardized accessibility scores. This is sometimes useful, but don't use it when comparing different modes.
    
    Args: 
        A_series(Series): Accessibility scores, with index as taz ids
    Returns: 
        Series: Standardized A scores. 
    """
    return A_series.apply(lambda x: (x-A_series.mean())/A_series.std()) 



if __name__ == '__main__':

    # create a directory for results, if doesn't already exist.
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # load walk times. 
    walktimes = pd.read_csv(os.path.join(root,data_dir, walk_file), dtype={'taz_id':str,'walk_time':float})
    walktimes = walktimes.fillna(0)
    walk_dict = dict(zip(walktimes['taz_id'],walktimes['walk_time']))

   
    # load land use activities data. 
    activ = load_activities(data_dir, lu_file)

    # make an index of taz labels
    taz_ind = activ[(activ.index!='0')&(pd.notnull(activ.index)&(activ.index!='None'))].index

    for lu_name in activ.columns:

        # load travel times by mode 
        mode_tables = {}

        for mode in ['myciti','taxi2011','taxi2015']:
            if mode=='taxi2015':
                infile = taxi_file_2015
            elif mode=='taxi2011':
                infile = taxi_file_2011
            elif mode=='myciti':
                infile = myciti_file
            
            costs_long = load_cost_matrix(data_dir,infile)
            costs_long = fill_missing_costs(costs_long, taz_ind)
            
            # add walk times for taxi (already added for MyCiTi)
            if (mode=='taxi2011')|(mode=='taxi2015'):
                costs_long['tot_cost'] = costs_long.apply(add_walk_time, axis=1)
                costs_long['agg_cost'] = costs_long['tot_cost']
            costs=costs_long.pivot(index='start',columns='end', values='agg_cost')
            costs = set_initial_zero(costs)
            mode_tables[mode]=costs

        ## Calculate A
        print('Calculating A for',len(taz_ind),'TAZs | land use:',lu_name)

        D = activ # dataframe with activities


        for mode in modes:
            df_tt = mode_tables[mode]
            for tt_max in tt_max_list:
                result = calculate_all_A(taz_ind = taz_ind, measure=measure, D=D, lu=lu_name, tt_max=tt_max, df_tt=df_tt, beta=beta)
                
                if standardize:
                    result = standardize_A(result)  # don't standardize if comparing different modes. 

                outfile = 'acc_{m}_{p}_{lu}.csv'.format(meas=measure, m=mode, p=tt_max,lu=lu_name)
                print('saving file as',outfile)
                result.to_csv(os.path.join(root,results_dir,outfile), index=True)



