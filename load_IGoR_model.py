#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:53:14 2017

@author: zacharysethna
"""

import numpy as np

#IGOR LOAD
def read_palindrome_parms(parms_file_name):
    parms_file = open(parms_file_name, 'r')
    
    
    in_delV = False
    in_delDl = False
    in_delDr = False
    in_delJ = False
    
    
    for line in parms_file:
        if line.startswith('#Deletion;V_gene;'):
            in_delV = True
            in_delDl = False
            in_delDr = False
            in_delJ = False
        elif line.startswith('#Deletion;D_gene;Three_prime;'):
            in_delV = False
            in_delDl = False
            in_delDr = True
            in_delJ = False
        elif line.startswith('#Deletion;D_gene;Five_prime;'):
            in_delV = False
            in_delDl = True
            in_delDr = False
            in_delJ = False
        elif line.startswith('#Deletion;J_gene;'):
            in_delV = False
            in_delDl = False
            in_delDr = False
            in_delJ = True
        elif any([in_delV, in_delDl, in_delDr, in_delJ]) and line.startswith('%'):
            if int(line.split(';')[-1]) == 0:
                if in_delV:
                    max_delV_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                elif in_delDl:
                    max_delDl_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                elif in_delDr:
                    max_delDr_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
                elif in_delJ:
                    max_delJ_palindrome = np.abs(int(line.lstrip('%').split(';')[0]))
        else:
            in_delV = False
            in_delDl = False
            in_delDr = False
            in_delJ = False
            
    return max_delV_palindrome, max_delDl_palindrome, max_delDr_palindrome, max_delJ_palindrome
        

def read_V_gene_parms(parms_file_name):
    parms_file = open(parms_file_name, 'r')
    
    V_gene_info = {}

    in_V_gene_sec = False
    for line in parms_file:
        if line.startswith('#GeneChoice;V_gene;'):
            in_V_gene_sec = True
        elif in_V_gene_sec:
            if line[0] == '%':
                split_line = line[1:].split(';')
                V_gene_info[split_line[0]] = [split_line[1] , int(split_line[2])]
            else:
                break
    parms_file.close()
    
    genV = [[]]*len(V_gene_info.keys())
    
    for V_gene in V_gene_info.keys():
        genV[V_gene_info[V_gene][1]] = [V_gene, '', V_gene_info[V_gene][0]]

    return genV
            
def read_D_gene_parms(parms_file_name):
    parms_file = open(parms_file_name, 'r')
    
    D_gene_info = {}

    in_D_gene_sec = False
    for line in parms_file:
        if line.startswith('#GeneChoice;D_gene;'):
            in_D_gene_sec = True
        elif in_D_gene_sec:
            if line[0] == '%':
                split_line = line[1:].split(';')
                D_gene_info[split_line[0]] = [split_line[1] , int(split_line[2])]
            else:
                break
    parms_file.close()
    
    genD = [[]]*len(D_gene_info.keys())
    
    for D_gene in D_gene_info.keys():
        genD[D_gene_info[D_gene][1]] = [D_gene, D_gene_info[D_gene][0]]

    return genD
    
def read_J_gene_parms(parms_file_name):
    parms_file = open(parms_file_name, 'r')
    
    J_gene_info = {}

    in_J_gene_sec = False
    for line in parms_file:
        if line.startswith('#GeneChoice;J_gene;'):
            in_J_gene_sec = True
        elif in_J_gene_sec:
            if line[0] == '%':
                split_line = line[1:].split(';')
                J_gene_info[split_line[0]] = [split_line[1] , int(split_line[2])]
            else:
                break
    parms_file.close()
    
    genJ = [[]]*len(J_gene_info.keys())
    
    for J_gene in J_gene_info.keys():
        genJ[J_gene_info[J_gene][1]] = [J_gene, '', J_gene_info[J_gene][0]]

    return genJ

def load_genomic_CDR3_anchor_pos_and_functionality(anchor_pos_file_name):
    
    anchor_pos_and_functionality = {}
    anchor_pos_file = open(anchor_pos_file_name, 'r')
    
    first_line = True
    for line in anchor_pos_file:
        if first_line:
            first_line = False
            continue
        
        split_line = line.split(',')
        split_line = [x.strip() for x in split_line]
        anchor_pos_and_functionality[split_line[0]] = [int(split_line[1]), split_line[2].strip().strip('()')]

    return anchor_pos_and_functionality


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
    
def generate_cutV_genomic_CDR3_segs(genomic_data):
    #max_palindrome = 6 #hard coded in for now
    try:
        max_palindrome = genomic_data['max_delV_palindrome']
    except KeyError:
        max_palindrome = 6 #default for matlab code
    #print 'V palindrome: ' + str(max_palindrome)
    cutV_genomic_CDR3_segs = []
    for CDR3_V_seg in [x[1] for x in genomic_data['genV']]:
        if len(CDR3_V_seg) < max_palindrome:
            cutV_genomic_CDR3_segs += [cutR_seq(CDR3_V_seg, 0, len(CDR3_V_seg))]
        else:
            cutV_genomic_CDR3_segs += [cutR_seq(CDR3_V_seg, 0, max_palindrome)]
    
    return cutV_genomic_CDR3_segs
    
def generate_cutJ_genomic_CDR3_segs(genomic_data):
    #max_palindrome = 6 #hard coded in for now
    try:
        max_palindrome = genomic_data['max_delJ_palindrome']
    except KeyError:
        max_palindrome = 6 #default for matlab code
    #print 'J palindrome: ' + str(max_palindrome)
    cutJ_genomic_CDR3_segs = []
    for CDR3_J_seg in [x[1] for x in genomic_data['genJ']]:
        if len(CDR3_J_seg) < max_palindrome:
            cutJ_genomic_CDR3_segs += [cutL_seq(CDR3_J_seg, 0, len(CDR3_J_seg))]
        else:
            cutJ_genomic_CDR3_segs += [cutL_seq(CDR3_J_seg, 0, max_palindrome)]
    
    return cutJ_genomic_CDR3_segs
    
def generate_cutD_genomic_CDR3_segs(genomic_data):
    #max_palindrome = 6 #hard coded in for now
    try:
        max_palindrome_L = genomic_data['max_delDl_palindrome']
    except KeyError:
        max_palindrome_L = 6 #default for matlab code
    try:
        max_palindrome_R = genomic_data['max_delDr_palindrome']
    except KeyError:
        max_palindrome_R = 6 #default for matlab code
    #print 'D palindrome: ' + str(max_palindrome_L) + ' ' + str(max_palindrome_R)
    cutD_genomic_CDR3_segs = []
    for CDR3_D_seg in [x[1] for x in genomic_data['genD']]:
        if len(CDR3_D_seg) < min(max_palindrome_L, max_palindrome_R):
            cutD_genomic_CDR3_segs += [cutR_seq(cutL_seq(CDR3_D_seg, 0, len(CDR3_D_seg)), 0, len(CDR3_D_seg))]
        else:
            cutD_genomic_CDR3_segs += [cutR_seq(cutL_seq(CDR3_D_seg, 0, max_palindrome_L), 0, max_palindrome_R)]
    
    return cutD_genomic_CDR3_segs
            
def load_genomic_data_igor_model(parms_file_name, V_anchor_pos_file, J_anchor_pos_file):
    genV = read_V_gene_parms(parms_file_name)
    genD = read_D_gene_parms(parms_file_name)
    genJ = read_J_gene_parms(parms_file_name)
    
    V_anchor_pos = load_genomic_CDR3_anchor_pos_and_functionality(V_anchor_pos_file)
    J_anchor_pos = load_genomic_CDR3_anchor_pos_and_functionality(J_anchor_pos_file)
    
    for V in genV:
        try:
            if V_anchor_pos[V[0]][0] > 0 and V_anchor_pos[V[0]][1] == 'F': #Check for functionality
                V[1] = V[2][V_anchor_pos[V[0]][0]:]
            else:
                V[1] = ''
        except KeyError:
            V[1] = ''

    for J in genJ:
        try:
            if J_anchor_pos[J[0]][0] > 0 and J_anchor_pos[J[0]][1] == 'F': #Check for functionality
                J[1] = J[2][:J_anchor_pos[J[0]][0]+3]
            else:
                J[1] = ''
        except KeyError:
            J[1] = ''
        
    max_delV_palindrome, max_delDl_palindrome, max_delDr_palindrome, max_delJ_palindrome = read_palindrome_parms(parms_file_name)
    genomic_data = {'genV': genV, 'genD': genD, 'genJ': genJ}
    genomic_data['max_delV_palindrome'] = max_delV_palindrome
    genomic_data['max_delDl_palindrome'] = max_delDl_palindrome
    genomic_data['max_delDr_palindrome'] = max_delDr_palindrome
    genomic_data['max_delJ_palindrome'] = max_delJ_palindrome
    genomic_data['cutV_genomic_CDR3_segs'] = generate_cutV_genomic_CDR3_segs(genomic_data)
    genomic_data['cutJ_genomic_CDR3_segs'] = generate_cutJ_genomic_CDR3_segs(genomic_data)
    genomic_data['cutD_genomic_CDR3_segs'] = generate_cutD_genomic_CDR3_segs(genomic_data)

    return genomic_data


def read_marginals_txt( filename , dim_names=False):
	with open(filename,'r') as file:
		#Model parameters are stored inside a dictionnary of ndarrays
		model_dict = {}
		dimension_names_dict = {}
		element_name=""
		first = True
		first_dim_line = False
		element_marginal_array = []
		indices_array = []

		for line in file:
			strip_line = line.rstrip('\n') #Remove end of line character
			if strip_line[0]=='@':
				first_dim_line = True
				if not(first):
					#Add the previous to the dictionnary
					model_dict[element_name] = element_marginal_array
				else:
					first = False
				
				element_name = strip_line[1:]
				#print element_name

			if strip_line[0]=='$':
				#define array dimensions
				coma_index = strip_line.find(',')
				dimensions = []

				#Get rid of $Dim[
				previous_coma_index = 4
				while coma_index != -1:
					dimensions.append(int(strip_line[previous_coma_index+1:coma_index]))
					previous_coma_index = coma_index
					coma_index = strip_line.find(',',coma_index+1)
			
				#Add last dimension and get rid of the closing bracket 
				dimensions.append(int(strip_line[previous_coma_index+1:-1]))

				element_marginal_array = np.ndarray(shape=dimensions)

			if strip_line[0]=='#':
				if first_dim_line:
					dimensions_names = []
					if len(dimensions) > 1:
						comma_index = strip_line.find(',')
						opening_bracket_index = strip_line.find('[')
						while opening_bracket_index != -1:
							dimensions_names.append(strip_line[opening_bracket_index+1:comma_index])
							opening_bracket_index = strip_line.find('[',comma_index) 
							comma_index = strip_line.find(',',opening_bracket_index)
					first_dim_line = False
					dimensions_names.append(element_name)
					dimension_names_dict[element_name] = dimensions_names
                    
                
				#update indices
				indices_array = []
				if len(dimensions) > 1:
					comma_index = strip_line.find(',')
					closing_brack_index = strip_line.find(']')					
					while closing_brack_index != -1:
						indices_array.append(int(strip_line[comma_index+1:closing_brack_index]))
						opening_bracket_index = strip_line.find('[',closing_brack_index) 
						comma_index = strip_line.find(',',opening_bracket_index)
						closing_brack_index = strip_line.find(']',closing_brack_index+1)
				

			if strip_line[0]=='%':
				#read doubles
				coma_index = strip_line.find(',')
				marginals_values = []

				#Get rid of the %
				previous_coma_index = 0
				while coma_index != -1:
					marginals_values.append(float(strip_line[previous_coma_index+1:coma_index]))
					previous_coma_index = coma_index
					coma_index = strip_line.find(',',coma_index+1)
			
				#Add last dimension and get rid of the closing bracket 
				marginals_values.append(float(strip_line[previous_coma_index+1:]))
				if len(marginals_values)!=dimensions[-1]:
					print "problem"
				element_marginal_array[tuple(indices_array)] = marginals_values
		model_dict[element_name] = element_marginal_array				
        
        
	return [model_dict,dimension_names_dict]

def load_and_process_igor_model(file_name):
    raw_model = read_marginals_txt(file_name)
    
    model = {}
    model['PV'] = raw_model[0]['v_choice']
    model['PinsVD'] = raw_model[0]['vd_ins']
    model['PinsDJ'] = raw_model[0]['dj_ins']
    model['PcutV_given_V'] = raw_model[0]['v_3_del'].T
    model['PcutJ_given_J'] = raw_model[0]['j_5_del'].T
    model['PDJ'] = np.multiply(raw_model[0]['d_gene'].T, raw_model[0]['j_choice'])
    model['PcutDlcutDr_given_D'] = np.transpose(np.multiply(np.transpose(raw_model[0]['d_3_del'], (2, 0, 1)), raw_model[0]['d_5_del']), (2, 0 , 1))
    Rvd = raw_model[0]['vd_dinucl'].reshape((4, 4)).T
    model['RnucleotideVD_per_nucleotideVD_5prime'] = np.multiply(Rvd, 1/np.sum(Rvd, axis = 0))
    Rdj = raw_model[0]['dj_dinucl'].reshape((4, 4)).T
    model['RnucleotideDJ_per_nucleotideDJ_3prime'] = np.multiply(Rdj, 1/np.sum(Rdj, axis = 0))
    
    #Independent insertion profiles
    model['PinsVDinsDJ'] = model['PinsVD'].reshape((model['PinsVD'].size, 1))*model['PinsDJ'].reshape((1, model['PinsDJ'].size))

    return model
