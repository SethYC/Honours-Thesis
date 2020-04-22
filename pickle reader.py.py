"""
Reads pickled results from multivaraite TE or MI outputs and creates
a file with the adjaceny matrix of the information bits. This is a light
modification of the original pickle reader file script.
    Created: July 25, 2019
    Updated: November 12, 2019
    Seth Campbell
"""

#Import classes
from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
from idtxl.visualise_graph import plot_network
import numpy as np
import re, os, pickle, datetime, time

#Init global variables
# pickle_path = "C:/Users/sethc/Documents/Summer Research 2019/misc/April Analysis/last mMI pickles/" #path for where to extract pickles
pickle_path = "C:/Users/sethc/Documents/Summer Research 2019/misc/April Analysis/last max lag 1 mTE pickles/" #path for where to extract pickles
output_path = "C:/Users/sethc/Documents/Summer Research 2019/misc/April Analysis/" #path for where to dump result file
all_adj_matrix = []
results = []

#takes a list of coordinate, value triplets, and the side length of the matrix, and creates an array out of it by falttening a matrix of the information values
def make_array(coordinates_triplets, matrix_length):
    matrix = np.zeros((matrix_length,matrix_length)) #init an appropriately sized arrays with zeroes

    for (x,y,value) in coordinates_triplets:
        matrix[x][y] = value #simply assign the value to the appropriate position in the matrix
    print(matrix)
    array = matrix.flatten() #flatten matrix to 1 dimension for easier handling in analysis (i.e. each subject is 1 row!)
    return array


##############################################################################################################################
#iterate over files in pickle_path & unpickle appropriate files to a list
with os.scandir(pickle_path) as files:
    for file in files:
        if re.search(r"filtcleants.*\.p", file.name): #only use the files that start with "filtcleants..." & ends in ".p"
        # if re.search(r"filtcleants_2011_dtd_c01 mTE.*\.p", file.name): #TESTING
            #note: results is a 2d array eith each entry containing a 2 item list with the unpickled result, and file name
            results.append([pickle.load(open(pickle_path + file.name, 'rb')),file.name]) #unpickle the file and add to a list

#extract raw adjacenty matrices from unpickled Results objects (from IDTXL package)
ROI_range = len((results[0][0].get_single_target(0,False))["sources_tested"])+1 #number of ROIS to look at based on the read data
subjects_matrix = np.zeros((len(results),(ROI_range**2)+1)) #init matrix of appropriate size with zeros, len(results+1) is the number of columns, plus columns for participant type
for number,result in enumerate(results):

    coordinates = []
    adj_matrix = result[0].get_adjacency_matrix(weights='max_te_lag', fdr=False) #get adj matrix instance (note that element 0 is accessed because each result in results is a list of two items: the result object and its file name)
    adj_matrix = adj_matrix._weight_matrix #get actual adj matrix of numbers
    print(adj_matrix)
    for target_num in range(0 , ROI_range): #iterate through all target nodes (i.e. all ROIs)
        # print(target_num)
        target_info = result[0].get_single_target(target=target_num, fdr=False) #get the info about sources for the current target, this returns a dictionary of information
        sig_sources = [i[0] for i in target_info["selected_vars_sources"]] #extract the y values of the significant source because they come in tuples paired with the time lag also (which is not needed)
        sig_sources = [(target_num,i) for i in sig_sources] #append the x axis component of the coordinate
        # print(target_info)
        te_or_mi_values   = target_info["selected_sources_te"] #get TE values (which are in the same order as the "selected_var_sources")
        # te_or_mi_values   = target_info["selected_sources_mi"] #get MI values (which are in the same order as the "selected_var_sources")

        for i,(j,k) in enumerate(sig_sources): #combine the coordinates (j,k) with the matching TE/MI value into a list of triplets
            coordinates += [(j,k,te_or_mi_values[i])] #triplet format: (x position, y position, TE/MI value)

    # setting group labels for SPSS or to fill in later, i.e. if a control patient then 1, if DTD then 2
    if re.search(r"dtd_d",result[1]):
        group_type = 2
    else: #controls either have "dtd_c" in the name or something else for the few exceptions from the new data from 5 participants or so
        group_type = 1

    subjects_matrix[number] = np.append(group_type,(make_array(coordinates, ROI_range))) #use coordinates to make array of information values for analysis
    # subjects_matrix[number][0] = results[1]

np.savetxt(output_path + "foo FINAL all subjects mTE info bits.csv", subjects_matrix, delimiter=",", fmt='%.8f') #save each matrix to a csv file
# np.savetxt(output_path + "FINAL all subjects mMI info bits.csv", subjects_matrix, delimiter=",", fmt='%.8f') #save each matrix to a csv file

# print(subjects_matrix)
print("\n~Completed combination of " + str(len(results)) + " adjacency matrix files\n")
