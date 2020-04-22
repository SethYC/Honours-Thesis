'''
The (true) final version of part 1 of analysis pipeline. Reads in fMRI csv files and computes
mMI or mTE measures on them, producing a pickle file for each particpant with the
results data object stored.

For the original 20 ROIs i choice, mMI/mTE took about 35 minutes or so to run.

April 21, 2020
Seth Campbell
'''

# Import classes
from idtxl.multivariate_mi import MultivariateMI
from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
# from idtxl.visualise_graph import plot_network
# import matplotlib.pyplot as plt
import numpy as np
# import networkx as nx
import time
import datetime
import pickle
import sys
import re
import os

#Init global variables
mMI_pickle_path = "C:/Users/sethc/Documents/Summer Research 2019/misc/April Analysis/last mMI pickles/"
mTE_pickle_path = "C:/Users/sethc/Documents/Summer Research 2019/misc/April Analysis/last max lag 1 mTE pickles/"
data_path = "C:/Users/sethc/Documents/Summer Research 2019/dtd_mist_data/"
all_adj_matrix = []
# type = 1 #1 for MI, 2 for mTE


'''
reads TSV file, parses into a list. Based on Ford's preprocessed files, the resulting format
is a 2D array, outer array is brain region, inner array is time series.
   input:   -ROI's as a list to analyze (based on MIST atlas),
            -file name to read
   return: numpy array of parsed data from TSV data file
'''
def set_up_data(ROI_nodes, file_to_read):
    n_nodes = len(ROI_nodes) #used to init blank 2D array and for looping through the data
    file = open(data_path + file_to_read,"r")
    columns = [[] for _ in range(n_nodes)]  #create a blank 2D array
    rows = file.readlines()
    headers = rows.pop(0)    #remove the headers of the tsv file (i.e. the brain region names)

    print(headers)

    for row in rows:  #read the rows to create new rows from the columns of each process (i.e. reshaping the data)
        j = 0   #acts as column counter
        for i in ROI_nodes:
            # print("i: "+str(i))
            columns[j].append(float(row.split('\t')[i-1]))
            j += 1

    data = np.array(columns) #convert to numpy array
    # print(data)
    return(data)


'''
performs one iteration of multivariate mutual information analysis on the data
    input:  -ROI's as a list to analyze (based on MIST atlas)
            -file_to_read as path of file to analyze
            -counter ?
    output: time to complete analysis
'''
def analyze_data(ROI_nodes, file_to_read, counter, analysis):
    start_time = time.time() #start timer
    n_nodes = len(ROI_nodes)

    raw_data = set_up_data(ROI_nodes, file_to_read) #parse tsv file and convert to formatted numpy array

    data = Data(raw_data, dim_order='ps')   #convert to Data type for use with IDTXL

    # b) Initialise analysis object and define settings
    if analysis == 1: network_analysis = MultivariateMI()
    else: network_analysis = MultivariateTE()

    settings = {'cmi_estimator': 'JidtGaussianCMI',
                'max_lag_sources': 5,
                'min_lag_sources': 1}

    results = network_analysis.analyse_network(settings=settings, data=data)
    time_elapsed = (time.time())-start_time


    adj_matrix = results.get_adjacency_matrix(weights='max_te_lag', fdr=False)
    # adj_matrix.print_matrix()
    # print(adj_matrix._weight_matrix)
    all_adj_matrix.append(adj_matrix._weight_matrix)
    # print(type(all_adj_matrix))
    # print(type(all_adj_matrix[0]))


    if type == 1:
        file_name = mMI_pickle_path + file_to_read[:-4] + " mMI - " + str(n_nodes) + \
            " ROIs results (min lag 1, max lag 5) pickled.p" #create pickle file name based on analyzed file's name
    else:
        file_name = mTE_pickle_path + file_to_read[:-4] + " mTE - " + str(n_nodes) + \
            " ROIs results (min lag 1, max lag 1) pickled.p"
    pickle.dump(results, open(file_name,'wb')) #pickle results to a file
    #results = pickle.load(open('results.p', 'rb')) #open the pickled file

    return(time_elapsed)


####################################################################
input = input("\nEnter 'm' for mMI analysis, and 't' for mTE analysis: ")
if input == "mMI" or input == "m":
    print("\nmMI analysis selected")
    type = 1
else:
    print("\nmTE analysis selected")
    type = 2
time.sleep(3) #delay to confirm correct type of analysis when executing on command line

times = []
file_list = []
ROI_nodes = [3,4,16,17,41,42,47,48,59,60,73,74,79,80,96,97,115,116,135,136]  #MIST atlas numbers for which ROIs to include in the analysis

with os.scandir(data_path) as files:
    for file in files:
        if re.search(r"filtcleants*", file.name): #only use the files that start with "filtcleants..."
            #print(file.name)
            file_list.append(file.name)

# print(file_list)

j = 0
for i,file in enumerate(file_list): #iterate through all files
    # j += 1
    # if j <= 15:
    #     continue
    # # if j >= 2:
    # #     break
    try:
        times.append(analyze_data(ROI_nodes,file,i,type))
    except:
        if type == 1:
            file = open(mMI_pickle_path + str(i) + file[:-4] + " - couldn't be analyzed","w")
            print("error: couldn't analyze it")
        else:
            file = open(mTE_pickle_path + str(i) + file[:-4] + " - couldn't be analyzed","w")
            print("error: couldn't analyze it")


avg_adj_matrix = all_adj_matrix[0]

# print(times)
