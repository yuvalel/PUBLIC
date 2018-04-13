#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 16:37:28 2018

@author: administrator
"""

def get_sharing_histogram(usage_filename, sharing_filename): 
    #builds a histogram of sharing numbers from usage file
    import numpy as np
    import pickle
    
    with open(usage_filename,'r') as f:

        N = 0 #N is how many individuals in usage file
        for line in f: #run over usage file once to get N
            splt = line.split()
            for entry in splt[1:]:
                N = max(N, int(entry.split(':')[0]))
        f.seek(0) #go back to start of usage file
            
        sharing_hist = np.zeros(N+1, dtype=int)
        for line in f:
            splt = line.split()
            share = 0
            for entry in splt[1:]: #count 1 for each entry in usage file for a CDR3 to get its sharing number
                share += 1 
            sharing_hist[share] += 1 #add one to sharing historgram for current sharing number 
    if sharing_filename!='NO_SAVE': #don't save sharing as file if file name given is 'NO_SAVE'
        with open(sharing_filename, 'wb') as f : #save sharing to pickle file
            pickle.dump(sharing_hist, f, protocol=pickle.HIGHEST_PROTOCOL)
    return sharing_hist
#%%
def get_sharing_nums(usage_filename): 
    #construct a vector of sharing numbers per sequence from usage file
    import numpy as np
    
    with open(usage_filename,'r') as f:

        N = 0 #N is how many individuals in usage file
        for k,line in enumerate(f): #run over usage file once to get N
            splt = line.split()
            for entry in splt[1:]:
                N = max(N, int(entry.split(':')[0]))
        f.seek(0) #go back to start of usage file
            
        sharing = np.zeros(k, dtype=int)
        for i,line in enumerate(f):
            splt = line.split()
            share = 0
            for entry in splt[1:]: #count 1 for each entry in usage file for a CDR3 to get its sharing number
                share += 1 
            sharing[i] = share #add one to sharing historgram for current sharing number 
    return sharing
#%%
def get_usage(samplesizes_filename, usage_filename, rawdata_folder):
    # creates a usage file of CDR3s from many files of different samples
    # input files are located in rawdata_folder
    # expected are column oriented tet files with the CDR3 in one colume and the status of the CDR3 in a different one (will take only seqeunces classified as 'In', for inframe)
    # first output is usage file, with one line for each CDR3 and multiple entries for every sample that contained this CDR3. sorted by CDR3
    # second output is sample size file, with list of samples and id numbers (used in the usage file to indetify the sample) and number of unique reads and CDR3s
    import os
    
    delmiter = '\t'
    CDR3aa_column = 1 #CDR3aa column
    status_column = 38 #seqeunces status column. set to -1 to not check status
    
    inframe_filenames = [filename for filename in os.listdir(rawdata_folder)]
    with open(samplesizes_filename, 'w') as output_file: #open new empty file for writing samples sizes
        pass
    
    for i,filename in enumerate(inframe_filenames):
        with open(rawdata_folder + filename, 'r') as data_file:
            CDR3s_dict = {}
            print "processing file: #" + str(i) + ' - ' + filename
            k = 0
            
            for line in data_file:
                splt = line.split(delmiter)
                if (status_column==-1) or (splt[status_column]=='In'): 
                    CDR3s_dict[splt[CDR3aa_column]] = CDR3s_dict.get(splt[CDR3aa_column],0) + 1 
                    k += 1
    
        with open(samplesizes_filename, 'a') as output_file:
            output_file.write(str(i) + ' ' + filename + ' ' + str(len(CDR3s_dict)) + ' ' + str(k) + '\n')
    
    
        if i==0:
            with open(usage_filename, 'w') as agg_file:
                for CDR3 in sorted(CDR3s_dict):
                    agg_file.write(CDR3 + ' 0:' + str(CDR3s_dict[CDR3]) + '\n')
        else:
            k =  0
            CDR3s_soreted_keys = sorted(CDR3s_dict)
            with open(usage_filename, 'r') as agg_file:
                with open(usage_filename + '.tmp', 'w') as temp_file:
                    for line in agg_file:
                        while (k<len(CDR3s_dict)) and (CDR3s_soreted_keys[k]<line.split()[0]): #if current CDR3 doesn't exist in file (if current line is after current CDR3)
                            temp_file.write(CDR3s_soreted_keys[k] + ' ' + str(i) + ':' + str(CDR3s_dict[CDR3s_soreted_keys[k]]) + '\n')
                            k += 1
                        if (k<len(CDR3s_dict)) and (line.split()[0]==CDR3s_soreted_keys[k]): #if found current CDR3 in file
                            temp_file.write(line.rstrip() + ' ' + str(i) + ':' + str(CDR3s_dict[CDR3s_soreted_keys[k]]) + '\n')
                            k += 1
                        else: #current line CDR3 doesn't exist in the new ind
                            temp_file.write(line)
                    while k<len(CDR3s_dict): #keep writing into file all the rest of the CDR3s in the current file
                        temp_file.write(CDR3s_soreted_keys[k] + ' ' + str(i) + ':' + str(CDR3s_dict[CDR3s_soreted_keys[k]]) + '\n')
                        k += 1
                            
            os.rename(usage_filename + '.tmp', usage_filename)
            
#%%
            
def ROC_curve(sharing, pgens, min_shr_pub):
    # returns ROC curve for PUBLIC over a sample, for a specific minimal sharing for publicness
    import numpy as np
    from sklearn import metrics
    ROC = np.array(metrics.roc_curve([x>min_shr_pub for x in sharing], pgens))
    return ROC

#%%
def AUROC_curve(sharing, pgens):
    # return AUC values for different values for the minimal sharing for being public as list of tuples
    import numpy as np
    from sklearn import metrics    
    max_share = max(sharing)
    values_AUC = np.sort(list(set([x.astype(int) for x in np.geomspace(1,max_share,200)]))) #values of minimal sharing in which to calcualte AUC
    AUC = []
    for min_shr_pub in values_AUC:
        sharing = [x>min_shr_pub for x in sharing]
        if not(all(sharing) or not any(sharing)):
            AUC.append((min_shr_pub, metrics.roc_auc_score(sharing, pgens)))
    return AUC


#%%
    
        
def main(num_samples, survival_prob, sequence_filename, marginals_file_name, parms_file_name, V_anchor_pos_file, J_anchor_pos_file):
    #generate seqeunces from a IGoR model and saves both full nt read and CDR3aa, with flat selection
    import sys
    import sequence_generation as seq_gen
    import load_IGoR_model as load_model
        
    genomic_data = load_model.load_genomic_data_igor_model(parms_file_name, V_anchor_pos_file, J_anchor_pos_file)
    genomic_data = seq_gen.construct_full_genomic_data(genomic_data)
    genomic_model = load_model.load_and_process_igor_model(marginals_file_name)
    CP_gm = seq_gen.compute_CP_gen_model_from_gen_model(genomic_model)

    k=0
    i=0
    
    output_file = open(sequence_filename, 'w')
    while k < num_samples:
        i += 1
        current_CDR3_seq, current_fullnt = seq_gen.gen_rnd_prod_CDR3_and_full_ntseq(CP_gm, genomic_data) #generate seq
        if hash(current_CDR3_seq) > ((1 - 2 * survival_prob)*sys.maxint): #seqeunce passes selection if hash value in the right range
            output_file.write(current_CDR3_seq + ' ' + current_fullnt + '\n')
            k += 1
        
    print 'generated ' + str(i) + ' recombinations, keeping and saving ' + str(k)
    output_file.close()
