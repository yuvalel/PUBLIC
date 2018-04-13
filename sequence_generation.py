#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:06:38 2017

@author: zacharysethna
"""

import numpy as np
#import sys
#sys.path.insert(0, '/mnt/biotheory/sequence_generation/')
#import load_IGoR_model as load_model

def construct_full_genomic_data(genomic_data):
    genomic_data['full_cutV_genomic_segs'] = generate_full_cutV_genomic_segs(genomic_data)
    genomic_data['full_cutJ_genomic_segs'] = generate_full_cutJ_genomic_segs(genomic_data)
    return genomic_data

def generate_full_cutV_genomic_segs(genomic_data):
    #max_palindrome = 6 #hard coded in for now
    try:
        max_palindrome = genomic_data['max_delV_palindrome']
    except KeyError:
        max_palindrome = 6 #default for matlab code
    #print 'V palindrome: ' + str(max_palindrome)
    full_cutV_genomic_segs = []
    for V_seg in [x[2] for x in genomic_data['genV']]:
        if len(V_seg) < max_palindrome:
            full_cutV_genomic_segs += [cutR_seq(V_seg, 0, len(V_seg))]
        else:
            full_cutV_genomic_segs += [cutR_seq(V_seg, 0, max_palindrome)]
    
    return full_cutV_genomic_segs


def generate_full_cutJ_genomic_segs(genomic_data):
    #max_palindrome = 6 #hard coded in for now
    try:
        max_palindrome = genomic_data['max_delJ_palindrome']
    except KeyError:
        max_palindrome = 6 #default for matlab code
    #print 'J palindrome: ' + str(max_palindrome)
    full_cutJ_genomic_segs = []
    for J_seg in [x[2] for x in genomic_data['genJ']]:
        if len(J_seg) < max_palindrome:
            full_cutJ_genomic_segs += [cutL_seq(J_seg, 0, len(J_seg))]
        else:
            full_cutJ_genomic_segs += [cutL_seq(J_seg, 0, max_palindrome)]
    
    return full_cutJ_genomic_segs


def compute_CP_gen_model_from_gen_model(genomic_model):  
    #Note that the genomic model is not perfectly normalized coming from IGoR. We need to renormalize certain arrays.
    
    
    CP_gen_model = {}
    
    CP_gen_model['CPV'] = (genomic_model['PV']/float(np.sum(genomic_model['PV']))).cumsum()
    CP_gen_model['CPDJ'] = (genomic_model['PDJ']/float(np.sum(genomic_model['PDJ']))).flatten().cumsum()
    CP_gen_model['CPinsVDinsDJ'] = (genomic_model['PinsVDinsDJ']/float(np.sum(genomic_model['PinsVDinsDJ']))).flatten().cumsum()
    
    for V in range(genomic_model['PcutV_given_V'].shape[1]):
        if np.sum(genomic_model['PcutV_given_V'][:, V])> 0:
            genomic_model['PcutV_given_V'][:, V] = genomic_model['PcutV_given_V'][:, V]/float(np.sum(genomic_model['PcutV_given_V'][:, V]))
    
    CP_gen_model['given_V_CPcutV'] = genomic_model['PcutV_given_V'].T.cumsum(axis = 1)
    
    for J in range(genomic_model['PcutJ_given_J'].shape[1]):
        if np.sum(genomic_model['PcutJ_given_J'][:, J])> 0:
            genomic_model['PcutJ_given_J'][:, J] = genomic_model['PcutJ_given_J'][:, J]/float(np.sum(genomic_model['PcutJ_given_J'][:, J]))

    CP_gen_model['given_J_CPcutJ'] = genomic_model['PcutJ_given_J'].T.cumsum(axis = 1)
    
    for D in range(genomic_model['PcutDlcutDr_given_D'].shape[2]):
        if np.sum(genomic_model['PcutDlcutDr_given_D'][:,  :, D]) > 0:
            genomic_model['PcutDlcutDr_given_D'][:, :, D] = genomic_model['PcutDlcutDr_given_D'][:, :, D]/float(np.sum(genomic_model['PcutDlcutDr_given_D'][:, :, D]))

    CP_gen_model['given_D_CPcutDlcutDr'] = np.array([ genomic_model['PcutDlcutDr_given_D'][:, :, i].flatten().cumsum() for i in range(genomic_model['PcutDlcutDr_given_D'].shape[2])])
    
    
    CP_gen_model['VD_Cntbias'] = genomic_model['RnucleotideVD_per_nucleotideVD_5prime'].T.cumsum(axis = 1)
    CP_gen_model['DJ_Cntbias'] = genomic_model['RnucleotideDJ_per_nucleotideDJ_3prime'].T.cumsum(axis = 1)
    CP_gen_model['CP_first_nt_VD'] = calc_steady_state_dist(genomic_model['RnucleotideVD_per_nucleotideVD_5prime']).cumsum()
    CP_gen_model['CP_first_nt_DJ'] = calc_steady_state_dist(genomic_model['RnucleotideDJ_per_nucleotideDJ_3prime']).cumsum()
 
    CP_gen_model['num_J_genes'] = genomic_model['PDJ'].shape[1]
    CP_gen_model['num_delDr_poss'] = genomic_model['PcutDlcutDr_given_D'].shape[1]
    CP_gen_model['num_DJ_ins_poss'] = genomic_model['PinsVDinsDJ'].shape[1]
    
    return CP_gen_model

def calc_steady_state_dist(R):
    #Calc steady state distribution for a dinucleotide bias matrix
    
    w, v = np.linalg.eig(R)
    
    for i in range(4):
        if np.abs(w[i] - 1) < 1e-8:
            return np.real(v[:, i] / np.sum(v[:, i]))
    return -1


#Scripts for sequence generation
def gen_rnd_prod_CDR3(CP_gm, genomic_data):
    #CP_gm = cumulative probability generative model
    #gd = genomic data
    coding_pass = False
    
    while ~coding_pass:
        ac = choose_random_assignment(CP_gm)
        #ac = assignment choice
        V_seq = genomic_data['cutV_genomic_CDR3_segs'][ac['V']]

#This both checks that the position of the conserved C is identified and that the V isn't fully deleted out of the CDR3 region
        if len(V_seq) <= max(ac['cutV'], 0): 
            continue
        D_seq = genomic_data['cutD_genomic_CDR3_segs'][ac['D']]
        J_seq = genomic_data['cutJ_genomic_CDR3_segs'][ac['J']]
 #We check that the D and J aren't deleted more than allowed. Note the generative model really should reflect this structure already
        if len(D_seq) < (ac['cutDl'] + ac['cutDr']) or len(J_seq) < ac['cutJ']:
            continue
        
        V_seq = V_seq[:len(V_seq) - ac['cutV']]
        D_seq = D_seq[ac['cutDl']:len(D_seq)-ac['cutDr']]
        J_seq = J_seq[ac['cutJ']:]
        
        if (len(V_seq)+ len(D_seq) + len(J_seq) + ac['insVD'] + ac['insDJ']) % 3 != 0:
            continue
        
        
        insVD_seq = rnd_ins_seq(ac['insVD'], CP_gm['VD_Cntbias'], CP_gm['CP_first_nt_VD'])
        insDJ_seq = rnd_ins_seq(ac['insDJ'], CP_gm['DJ_Cntbias'], CP_gm['CP_first_nt_DJ'])[::-1] #have to reverse the DJ seq
        
        #Translate to amino acid sequence, see if productive
        aaseq = nt2aa(V_seq + insVD_seq + D_seq + insDJ_seq + J_seq)
        
        if '*' not in aaseq and aaseq[0]=='C' and aaseq[-1] =='F':
            return aaseq

def gen_rnd_prod_CDR3_ntseq(CP_gm, genomic_data):
    #CP_gm = cumulative probability generative model
    #gd = genomic data
    coding_pass = False
    
    while ~coding_pass:
        ac = choose_random_assignment(CP_gm)
        #ac = assignment choice
        V_seq = genomic_data['cutV_genomic_CDR3_segs'][ac['V']]

#This both checks that the position of the conserved C is identified and that the V isn't fully deleted out of the CDR3 region
        if len(V_seq) <= max(ac['cutV'], 0): 
            continue
        D_seq = genomic_data['cutD_genomic_CDR3_segs'][ac['D']]
        J_seq = genomic_data['cutJ_genomic_CDR3_segs'][ac['J']]
 #We check that the D and J aren't deleted more than allowed. Note the generative model really should reflect this structure already
        if len(D_seq) < (ac['cutDl'] + ac['cutDr']) or len(J_seq) < ac['cutJ']:
            continue
        
        V_seq = V_seq[:len(V_seq) - ac['cutV']]
        D_seq = D_seq[ac['cutDl']:len(D_seq)-ac['cutDr']]
        J_seq = J_seq[ac['cutJ']:]
        
        if (len(V_seq)+ len(D_seq) + len(J_seq) + ac['insVD'] + ac['insDJ']) % 3 != 0:
            continue
        
        
        insVD_seq = rnd_ins_seq(ac['insVD'], CP_gm['VD_Cntbias'], CP_gm['CP_first_nt_VD'])
        insDJ_seq = rnd_ins_seq(ac['insDJ'], CP_gm['DJ_Cntbias'], CP_gm['CP_first_nt_DJ'])[::-1] #have to reverse the DJ seq
        
        #Translate to amino acid sequence, see if productive
        ntseq = V_seq + insVD_seq + D_seq + insDJ_seq + J_seq
        aaseq = nt2aa(V_seq + insVD_seq + D_seq + insDJ_seq + J_seq)
        
        if '*' not in aaseq and aaseq[0]=='C' and aaseq[-1] =='F':
            return ntseq
        
def gen_rnd_prod_CDR3_and_full_ntseq(CP_gm, genomic_data):
    #CP_gm = cumulative probability generative model
    #gd = genomic data
    coding_pass = False
    
    while ~coding_pass:
        ac = choose_random_assignment(CP_gm)
        #ac = assignment choice
        V_seq = genomic_data['cutV_genomic_CDR3_segs'][ac['V']]
        
#This both checks that the position of the conserved C is identified and that the V isn't fully deleted out of the CDR3 region
        if len(V_seq) <= max(ac['cutV'], 0): 
            continue
        D_seq = genomic_data['cutD_genomic_CDR3_segs'][ac['D']]
        J_seq = genomic_data['cutJ_genomic_CDR3_segs'][ac['J']]
 #We check that the D and J aren't deleted more than allowed. Note the generative model really should reflect this structure already
        if len(D_seq) < (ac['cutDl'] + ac['cutDr']) or len(J_seq) < ac['cutJ']:
            continue
        
        V_seq = V_seq[:len(V_seq) - ac['cutV']]
        D_seq = D_seq[ac['cutDl']:len(D_seq)-ac['cutDr']]
        J_seq = J_seq[ac['cutJ']:]
        
        if (len(V_seq)+ len(D_seq) + len(J_seq) + ac['insVD'] + ac['insDJ']) % 3 != 0:
            continue
        
        
        insVD_seq = rnd_ins_seq(ac['insVD'], CP_gm['VD_Cntbias'], CP_gm['CP_first_nt_VD'])
        insDJ_seq = rnd_ins_seq(ac['insDJ'], CP_gm['DJ_Cntbias'], CP_gm['CP_first_nt_DJ'])[::-1] #have to reverse the DJ seq
        
        #Translate to amino acid sequence, see if productive
        aaseq = nt2aa(V_seq + insVD_seq + D_seq + insDJ_seq + J_seq)
        
        full_V_seq = genomic_data['full_cutV_genomic_segs'][ac['V']]
        full_V_seq = full_V_seq[:len(full_V_seq) - ac['cutV']]
        
        full_J_seq = genomic_data['full_cutJ_genomic_segs'][ac['J']][ac['cutJ']:]
        
        full_ntseq = full_V_seq + insVD_seq + D_seq + insDJ_seq + full_J_seq
        
        
        if '*' not in aaseq and aaseq[0]=='C':
            return aaseq, full_ntseq


def choose_random_assignment(CP_gm):
    #CP_gm = cumulative probability generative model
    #Choose the assignment variables
    
    #Need to specify these parameters as CP_gm feeds in flattened arrays to be called
    #num_J_genes = 14 #gm['PDJ'].shape[1]
    #num_delDr_poss = 23 #gm['PcutDlcutDr_given_D'].shape[1]
    #num_DJ_ins_poss = 31 #gm['PinsVDinsDJ'].shape[1]
    
                            
    assignment_choice = {}
    assignment_choice['V'] = CP_gm['CPV'].searchsorted(np.random.random())
    
    #For 2D arrays make sure to take advantage of a mod expansion to find indicies
    DJ_choice = CP_gm['CPDJ'].searchsorted(np.random.random())
    assignment_choice['D'] = DJ_choice/CP_gm['num_J_genes']
    assignment_choice['J'] = DJ_choice % CP_gm['num_J_genes']    
    #assignment_choice['D'] = DJ_choice/gm['PDJ'].shape[1]
    #assignment_choice['J'] = DJ_choice % gm['PDJ'].shape[1]
    
    
    #Refer to the correct slices for the dependent distributions
    assignment_choice['cutV'] = CP_gm['given_V_CPcutV'][assignment_choice['V'], :].searchsorted(np.random.random())
    
    assignment_choice['cutJ'] = CP_gm['given_J_CPcutJ'][assignment_choice['J'], :].searchsorted(np.random.random())
    
    cutDlcutDr_choice = CP_gm['given_D_CPcutDlcutDr'][assignment_choice['D'], :].searchsorted(np.random.random()) 
    
    assignment_choice['cutDl'] = cutDlcutDr_choice/CP_gm['num_delDr_poss']
    assignment_choice['cutDr'] = cutDlcutDr_choice % CP_gm['num_delDr_poss']
    
    insVDinsDJ_choice = CP_gm['CPinsVDinsDJ'].searchsorted(np.random.random())
    
    assignment_choice['insVD'] = insVDinsDJ_choice/CP_gm['num_DJ_ins_poss']
    assignment_choice['insDJ'] = insVDinsDJ_choice % CP_gm['num_DJ_ins_poss']
    
    return assignment_choice


def cutR_seq(seq, cutR, max_palindrome):
    #max_palindrome = 6 #for convenience max_palindrome has been hard coded in, if max_palindrome is not 6 change here.
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #can include lower case if wanted
    if cutR < max_palindrome:
        seq = seq + ''.join([complement_dict[nt] for nt in seq[cutR - max_palindrome:]][::-1]) #reverse complement palindrome insertions
    else:
        seq = seq[:len(seq) - cutR + max_palindrome] #deletions
    
    return seq
    
def cutL_seq(seq, cutL, max_palindrome):
    #max_palindrome = 6 #for convenience max_palindrome has been hard coded in, if max_palindrome is not 6 change here.
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #can include lower case if wanted
    if cutL < max_palindrome:
        seq = ''.join([complement_dict[nt] for nt in seq[:max_palindrome - cutL]][::-1]) + seq #reverse complement palindrome insertions
    else:
        seq = seq[cutL-max_palindrome:] #deletions
    
    return seq


def rnd_ins_seq(ins_len, nt_bias, CP_first_nt):
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num2nt = 'ACGT'
    
    if ins_len == 0:
        return ''
    
    seq = num2nt[CP_first_nt.searchsorted(np.random.random())]
    ins_len += -1
    
    while ins_len > 0:
        seq += num2nt[nt_bias[nt2num[seq[-1]], :].searchsorted(np.random.random())]
        ins_len += -1
    
    return seq

    
def nt2aa(ntseq, frame = 0):
    ntseq = ntseq[frame:]
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    aa_dict ='KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF'
    
    return ''.join([aa_dict[nt2num[ntseq[i]] + 4*nt2num[ntseq[i+1]] + 16*nt2num[ntseq[i+2]]] for i in range(0, len(ntseq), 3) if i+2 < len(ntseq)])