# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:12:23 2020

@author: sscha
"""
import numpy as np 
import copy
import scipy.integrate as spi
from sympy import symbols, Matrix, Transpose, matrix_multiply_elementwise as MME
import networkx as nx
import matplotlib.pyplot as plt
import sys
import xlrd

''' 
###############################################################################
general genelet model functions (helper functions used within GeneletNetwork)
###############################################################################
'''

# function to assemble connectivity matrices
def con_vec2mat(con_vec,ortho_nodes,ind_nodes):
    '''
    This function converts connectivity vectors into connectivity matrices
    '''
    con_mat = np.zeros([ortho_nodes,ind_nodes], dtype = int)
    
    for n in range(len(con_vec)):
        if con_vec[n] != 0:
            con_mat[con_vec[n]-1,n] = 1
 
    
    return con_mat

# function to assemble production matrices
def prod_vec2mat(prod_vec,ortho_nodes,ind_nodes):

    Cprod_mat = np.zeros([ortho_nodes,ind_nodes], dtype = int)
    Rprod_mat = np.zeros([ortho_nodes,ind_nodes], dtype = int)
    
    for n in range(len(prod_vec)):
        if prod_vec[n] < 0:
            Rprod_mat[-prod_vec[n]-1,n] = 1
        elif prod_vec[n] > 0:
            Cprod_mat[prod_vec[n]-1,n] = 1
        
    return Cprod_mat, Rprod_mat

# function to create a topology matrix for topology plotting
def to_topology_mat(act_vec,prod_vec):
    
    ortho_nodes = max(act_vec)
    ind_nodes = len(act_vec)
    
    topology_mat = [[0]*ortho_nodes for i in range(ortho_nodes)]
    
    for i in range(ind_nodes):
        
        if prod_vec[i] > 0:
            topology_mat[abs(prod_vec[i])-1][act_vec[i]-1] = 1
        elif prod_vec[i] < 0:
            topology_mat[abs(prod_vec[i])-1][act_vec[i]-1] = -1
            
    return topology_mat

# function to create act_vec and prod_vec (of indc_vec) from a topology mat
def from_topology_mat(topology_mat,ind_n):
                            
    k = 0
    act_vec =[len(topology_mat[0])]*ind_n
    prod_vec = [0]*ind_n
    for i in range(len(topology_mat[0])):
        for j in range(len(topology_mat[0])):
        
            tmv = topology_mat[j][i]
            
            if tmv == 1:
                prod_vec[k]=(j+1)
                act_vec[k]=(i+1)
                k += 1
            if tmv == -1:
                prod_vec[k]=-(j+1)
                act_vec[k]=(i+1)
                k += 1
            if (j+1) == len(topology_mat[0]) and sum(abs(np.array(topology_mat)[:,i]))==0:
                k += 1
                    
    return act_vec, prod_vec
                
# function to create genelet initial conditions
def int_genelet_states(G_int_vec,G_tot,act_vec,dA_tot):
    
    GdAin = np.zeros(len(G_int_vec), dtype = int)
    GdBin = np.zeros(len(G_int_vec), dtype = int)
    dAin = copy.deepcopy(dA_tot)
    act_vec = np.array(act_vec)
    
    for n in range(len(G_int_vec)):
        
        # where there is more activator than total genelet of that activator
        if G_int_vec[n] == 1 and dAin[act_vec[n]-1] >= sum(G_tot[act_vec==act_vec[n]]):
            GdAin[n] = G_tot[n]
            dAin[act_vec[n]-1] = dAin[act_vec[n]-1] - G_tot[n]
            
        # where there is less activator than total genelet
        elif G_int_vec[n] == 1 and dAin[act_vec[n]-1] < sum(G_tot[act_vec==act_vec[n]]):
            # currently not sure how to handle this so an error will be printed
            print('Warning: Less activator than genelets!!!')
            '''GdAin[n] = dAin[act_vec[n]-1]/sum(act_vec==act_vec[n])
               dAin[act_vec[n]-1] = 0'''
        
        # blocked genelets
        # (excess blocker concentrations are handled in the simulation script)
        elif G_int_vec[n] == -1:
            GdBin[n] = G_tot[n]
            
    return GdAin, GdBin, dAin

# function to determine reverse complement (rc), reverse (r), or complement (c) of a nucleotide sequence
def rev_comp_seq(seq,spec_out='rc'):
    
    new_seq = '3'+seq+'5'

    new_seq = list(new_seq)
    ds = ''
    
    for i in range(len(seq)):
        if seq[i].lower()=='g':
            new_seq[i+1] = 'C'
        elif seq[i].lower()=='c':
            new_seq[i+1] = 'G'
        elif seq[i].lower()=='t':
            new_seq[i+1] = 'A'
        elif seq[i].lower()=='a':
            new_seq[i+1] = 'T'
        
    if spec_out.lower() == 'rc':
        out_seq = ds.join(new_seq[::-1]) 
    elif spec_out.lower() == 'c':
        out_seq=ds.join(new_seq)   
    elif spec_out.lower() == 'r':
        out_seq = '3'+seq[::-1]+'5'
        
    return out_seq

def xlsheet_to_dict(input_file_loc,sheet_name,ncol):
    
    seq_workbook = xlrd.open_workbook(input_file_loc)
    sheet = seq_workbook.sheet_by_name(sheet_name)
    sheet_dict = {}
    dv = [[0]*ncol for i in range(0,sheet.nrows)]
    for ri in range(0,sheet.nrows):
        dk=sheet.cell_value(ri,0)
        a = 0
        for ci in range(0,ncol):
            dv[ri][a]=sheet.cell_value(ri,ci+1)
            a += 1

        sheet_dict.update({dk : dv[ri]})
        
    return sheet_dict

# function containing the ODEs for the Kim and Winfree genelets (KWG)
def general_genelet_eqs(t,x,ortho_nodes,ind_nodes,dA_tot,dB_tot,G_tot, \
                        kpr,kpc,kpi,kd_H,kd_A,kga,kgar,kar,kgb,kgab,kgbc,kbc,kir, \
                        act_mat,rep_mat,blk_mat,ca_mat,Cprod_mat,Rprod_mat, \
                        Cindc_mat,Rindc_mat,RnH,RnA,leak):
    
    k = ortho_nodes
    g = ind_nodes
    
    # defining species for convenience
    rR = x[0:k] # repressors
    dA = x[k:2*k] # activators
    rC = x[2*k:3*k] # coactivators
    dB = x[3*k:4*k] # blockers
    
    GdA = x[4*k:4*k+g] # genelet:activator complexes
    GdB = x[4*k+g:4*k+2*g] # genelet:blocker complexes
    
    rIr = x[4*k+2*g:4*k+2*g+k] # repressor inducers 
    rIc = x[4*k+2*g+k:4*k+2*g+2*k] # coactivator inducers
    
    dR = x[4*k+2*g+2*k:4*k+2*g+3*k] # DNA repressors
    dAdR = x[4*k+2*g+3*k:4*k+2*g+4*k] # DNA repressor:DNA activator complexes which cannot be degraded
    
    # mass balances to find other concentrations
    dArR = dA_tot - dA - dAdR - np.matmul(act_mat,GdA) # activator:repressor complexes
    dBrC = dB_tot - dB - np.matmul(blk_mat,GdB) # blocker:coactivator complexes
    G = G_tot - GdA - GdB # OFF genelets
    
    # Rate equations:
    
    # repressors
    dRdt = -kgar*(rep_mat@GdA)*rR - kar*dA*rR + Rprod_mat@(kpr*GdA) + Rprod_mat@(leak*kpr*GdB) - kir*rR*rIr - RnA*kd_A*rR
    
    # activators
    dAdt = RnH*kd_H*dArR+ RnA*kd_A*dArR - kga*(act_mat@G)*dA - kar*dA*rR + kgab*(blk_mat@GdA)*dB - kar*dA*dR 
    
    # coactivators
    dCdt = -kgbc*(ca_mat@GdB)*rC - kbc*dB*rC + Cprod_mat@(kpc*GdA) + Cprod_mat@(leak*kpc*GdB) - kir*rC*rIc - RnA*kd_A*rC
    
    # blockers
    dBdt = RnH*kd_H*dBrC + RnA*kd_A*dBrC - kgb*(blk_mat@G)*dB - kbc*dB*rC - kgab*(blk_mat@GdA)*dB
    
    # G (ON)
    dGondt = (act_mat.T@(kga*dA))*G - (rep_mat.T@(kgar*rR))*GdA - (blk_mat.T@(kgab*dB))*GdA - (rep_mat.T@(kgar*dR))*GdA 
    
    #G (BLK)
    dGblkdt = (blk_mat.T@(kgb*dB))*G + (blk_mat.T@(kgab*dB))*GdA - (ca_mat.T@(kgbc*rC))*GdB
    
    # repressor inducers
    dIrdt = Rindc_mat@(kpi*GdA) + Rindc_mat@(leak*kpi*GdB) - kir*rR*rIr - RnA*kd_A*rIr
    
    # coactivator inducers
    dIcdt = Cindc_mat@(kpi*GdA) + Cindc_mat@(leak*kpi*GdB) - kir*rC*rIc - RnA*kd_A*rIc
    
    # DNA repressor
    ddRdt = - kar*dA*dR  + kgar*(rep_mat@GdA)*dR
    
    # DNA repressor:DNA activator complexes
    ddAdRdt = kar*dA*dR + kgar*(rep_mat@GdA)*dR
    
    return np.concatenate([dRdt,dAdt,dCdt,dBdt,dGondt,dGblkdt,dIrdt,dIcdt,ddRdt,ddAdRdt])

# function containing the ODEs for spatiotemporal genelets (STG)
def spatial_genelet_eqs(t,x,ortho_nodes,ind_nodes,dA_tot,G_tot, \
                        kpr,kpi,kd_H,kd_A,kd_Hg,kd_Ag,kgar,kar,kir, \
                        act_mat,rep_mat,Rprod_mat,Rindc_mat,RnH,RnA,leak):
    
    k = ortho_nodes
    g = ind_nodes
    
    # defining species for convenience
    rR = x[0:k] # repressors
    dA = x[k:2*k] # activators
    
    Gon = x[2*k:2*k+g] # active genelets complexes
    
    rIr = x[2*k+g:2*k+g+k] # repressor inducers 
    
    # mass balances to find other concentrations
    dArR = dA_tot - dA # free activator:repressor complexes
    GrR = G_tot - Gon # OFF genelets
    
    # Rate equations:
    
    # repressors
    dRdt = -kgar*(rep_mat@Gon)*rR - kar*dA*rR + Rprod_mat@(kpr*Gon) + Rprod_mat@(leak*kpr*GrR) - kir*rR*rIr - RnA*kd_A*rR
    
    # activators
    dAdt = RnH*kd_H*dArR + RnA*kd_A*dArR - kar*dA*rR
    
    # G (ON)
    dGondt = RnH*kd_Hg*GrR + RnA*kd_Ag*GrR - (rep_mat.T@(kgar*rR))*Gon
    
    # repressor inducers
    dIrdt = Rindc_mat@(kpi*Gon) + Rindc_mat@(leak*kpi*GrR) - kir*rR*rIr - RnA*kd_A*rIr
    
    return np.concatenate([dRdt,dAdt,dGondt,dIrdt])

'''
###############################################################################
GeneletNetwork CLASS DEFINITION
###############################################################################
'''
class GeneletNetwork:
    def __init__(self,act_vec='req1',prod_vec='req1',indc_vec=[],blk_vec=[],top_mat='req2',Itop_mat=[],genelet_type='KWG'):
        
        self.genelet_type = genelet_type                
        self.exit = 0
        
        # if you don't use topology matrix notation you have to input at least
        # act_vec and prod_vec
        if top_mat == 'req2':
            if act_vec  == 'req1' or prod_vec == 'req1':
                print(" Need inputs: 'act_vec' and 'prod_vec' or 'top_mat' ")

            else:
                self.topology_mat = to_topology_mat(act_vec,prod_vec)
                ortho_nodes = max(act_vec)
                ind_nodes = len(act_vec)
                if indc_vec == []:
                    indc_vec=np.zeros(ind_nodes, dtype = int)
                if Itop_mat == []:
                    self.I_topology_mat = to_topology_mat(act_vec,indc_vec)
        
        # if you use topology matrix notation you have to input topology matrix
        # if the inducer topology matrix is not supplied it is assumed to be all 0s        
        elif top_mat != 'req2':
            
            self.topology_mat = top_mat
            ortho_nodes = len(top_mat[0])
            
            if Itop_mat == []:
                Itop_mat = [[0]*ortho_nodes for a in range(ortho_nodes)]
                self.I_topology_mat = Itop_mat
            else:
                self.I_topology_mat = Itop_mat  

            # counting total number of unique RNA's produced
            np_top_mat = np.array(top_mat)
            np_Itop_mat = np.array(Itop_mat)
            ind_nodes = ortho_nodes # initialization
            # finding the total ind nodes, each column that has more than a single 
            # 1 in it produces more than a single RNA
            for i in range(ortho_nodes):
                if sum(abs(np_top_mat[:,i]))+sum(abs(np_Itop_mat[:,i])) > 1:
                    ind_nodes += sum(abs(np_top_mat[:,i]))+sum(abs(np_Itop_mat[:,i]))-1
            
            # create act_vec and prod_vec from top_mat
            act_vec,prod_vec = from_topology_mat(top_mat,ind_nodes)
            # create indc_vec from Itop_mat
            _,indc_vec = from_topology_mat(Itop_mat,ind_nodes)
              
                
        # computing everything else now that we know the topology
        
        self.act_vec = act_vec
        self.ortho_nodes = ortho_nodes
        self.ind_nodes = ind_nodes
            
        if blk_vec == []:
            # if the blk vector is not included as an input then this vector will
            # be populated assuming that any nodes that have coactivators coming to them
            # have blockers and those without coactivators do not
            c = 0
            np_top_mat = np.array(self.topology_mat)
            blk_vec=np.zeros(self.ind_nodes, dtype = int)
            for i in range(self.ortho_nodes):
                np_mat = np.array(self.topology_mat[i])
                if sum(abs(np_mat))==0:
                    # node with no connections to it will no have blocker
                    if sum(abs(np_top_mat[:,i])) == 0:
                        c = c + 1
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        c = c + 1
                elif any(np_mat==1):
                    # If a coactivator is produced for a node then it will be blocked initially
                    if sum(abs(np_top_mat[:,i])) == 0:
                        blk_vec[c] = i+1
                        c = c + 1
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        blk_vec[c] = i+1
                        c = c + 1
                elif sum(np_mat)<0 and all(np_mat!=1):
                    # nodes that are only repressed will no have blockers
                    if sum(abs(np_top_mat[:,i])) == 0:
                        c = c + 1
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        c = c + 1
                             
        self.blk_vec = blk_vec
        
        self.rep_vec = copy.deepcopy(self.act_vec)
        self.ca_vec = copy.deepcopy(self.blk_vec)
        
        # Connectivity matrices
        self.act_mat = con_vec2mat(self.act_vec,self.ortho_nodes,self.ind_nodes)
        self.blk_mat = con_vec2mat(self.blk_vec,self.ortho_nodes,self.ind_nodes)
        self.rep_mat = con_vec2mat(self.rep_vec,self.ortho_nodes,self.ind_nodes)
        self.ca_mat = con_vec2mat(self.ca_vec,self.ortho_nodes,self.ind_nodes)
        
        # Production vectors
        self.prod_vec = prod_vec
        self.indc_vec = indc_vec
        
        # defining the production matrices
        self.Cprod_mat,self.Rprod_mat = prod_vec2mat(self.prod_vec,self.ortho_nodes,self.ind_nodes)
        self.Cindc_mat,self.Rindc_mat = prod_vec2mat(self.indc_vec,self.ortho_nodes,self.ind_nodes)
            
    ''' DEFINING INITIAL CONDITIONS '''
    def initial_conditions(self,dA_tot,G_tot,G_int_vec=[],dB_added=[],rRin=[],rCin=[],rIrin=[],rIcin=[],dRin=[],dAdRin=[]):
    
        # default initial genelet states are deteremined from the topology matrix
        # if not provided by the user
        np_top_mat = np.array(self.topology_mat)
        if G_int_vec == []:
            for i in range(self.ortho_nodes):
                np_mat = np.array(self.topology_mat[i])
                if sum(abs(np_mat))==0:
                    # node with no connections to it will be ON initally
                    if sum(abs(np_top_mat[:,i])) == 0:
                        G_int_vec.append(1)
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        G_int_vec.append(1)
                elif any(np_mat==1):
                    # If a coactivator is produced for a node then it will be blocked initially
                    if sum(abs(np_top_mat[:,i])) == 0:
                        G_int_vec.append(-1)
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        G_int_vec.append(-1)
                elif sum(np_mat)<0 and all(np_mat!=1):
                    # node that is only repressed will be ON initally
                    if sum(abs(np_top_mat[:,i])) == 0:
                        G_int_vec.append(1)
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        G_int_vec.append(1)
            
        # defining default zero values for all of these inputs
        if dB_added == []:
            dB_added=np.zeros(self.ortho_nodes, dtype = int)
        if rRin == []:
            rRin=np.zeros(self.ortho_nodes, dtype = int)
        if rCin == []:
            rCin=np.zeros(self.ortho_nodes, dtype = int)
        if rIrin == []:
            rIrin=np.zeros(self.ortho_nodes, dtype = int)
        if rIcin == []:
            rIcin=np.zeros(self.ortho_nodes, dtype = int)
        if dRin == []:
            dRin=np.zeros(self.ortho_nodes, dtype = int)
        if dAdRin == []:
            dAdRin=np.zeros(self.ortho_nodes, dtype = int)
        
        self.dA_tot = dA_tot # total activator added
                        
        self.G_tot = G_tot # total genelet added
        
        # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
        self.G_int_vec = G_int_vec
        
        if self.genelet_type.lower() == 'kwg':
        
            [self.GdAin,self.GdBin,self.dAin] = int_genelet_states(self.G_int_vec,self.G_tot,self.act_vec,self.dA_tot)
            
            self.dB_anneal = np.matmul(self.blk_mat,self.GdBin) # blocker attached to the annealed genelets
            self.dB_ex_anneal = 0.5*self.dB_anneal # excess blocker from the anneal
            self.dB_added = dB_added # added blocker beyond that present from the anneal
            self.dB_tot = self.dB_anneal + self.dB_ex_anneal + self.dB_added # total blocker present 
            
            self.rRin = rRin # initial repressor concentration
            self.rCin = rCin # initial coactivator concentration
            self.dBin = self.dB_ex_anneal + self.dB_added # initial free blocker concentration
            self.rIrin = rIrin # initial repressor inducer concentration
            self.rIcin = rIcin # initial coactivator inducer concentration
            self.dRin = dRin # initial DNA repressor concentration
            self.dAdRin = dAdRin # initial DNA repressor:DNA activator concentration
            
            # final initial conditions
            self.int_cond = np.concatenate([self.rRin,self.dAin,self.rCin,\
                                            self.dBin,self.GdAin,self.GdBin,\
                                            self.rIrin,self.rIcin,self.dRin,self.dAdRin])
    
        elif self.genelet_type.lower() == 'stg':
        # be careful here about excess activator and repressors when setting inital contions
            self.Gon=np.zeros(self.ind_nodes, dtype = int)
            
            for i in range(len(self.G_int_vec)):
                if self.G_int_vec[i] == 1:
                    self.Gon[i] = self.G_tot[i]
                
            self.rRin = rRin # initial repressor concentration
            self.rIrin = rIrin # initial repressor inducer concentration
            
            # final initial conditions
            self.int_cond = np.concatenate([self.rRin,self.dA_tot,self.Gon,self.rIrin])

    
    ''' DEFINING PARAMETERS FOR SIMULATION AND SIMULATING MODEL '''
    def simulate(self,t_vec,iteration,rnase='RnH',leak=0,rate_constants=[],rR=[],dA=[],rC=[],dB=[],rIr=[],rIc=[],dR=[]):
        ''' Currently there is no capability to add genelet species mid experiment
            to add more G then update G_tot and to add GdB update G_tot and dB_tot
        '''
        # defining which rnases are being used in the simulation
        if rnase.lower() == 'rnh':
            self.RnH = 1
            self.RnA = 0
        elif rnase.lower() == 'rna':
            self.RnH = 0
            self.RnA = 1
        elif rnase.lower() == 'both':
            self.RnH = 1
            self.RnA = 1
        
        if isinstance(leak,list) and leak[0] <= 1:
            self.leak = leak
        elif isinstance(leak,list) and leak[0] > 1:
            self.leak = [e/100 for e in leak] # converting an input percent into a fraction
        elif leak <= 1:
            self.leak = leak
        elif leak > 1:
            self.leak = leak/100 # converting an input percent into a fraction
        
        if len(rate_constants) == 0:
            # defining the rate constants (not sure the best place to put these?)
            self.kpr = 0.02*np.ones(self.ind_nodes) # repressor production rates
            self.kpc = 0.02*np.ones(self.ind_nodes) # coactivator production rates
            self.kpi = 0.02*np.ones(self.ind_nodes) # inducer production rates
            self.kd_H = 0.0003*np.ones(self.ortho_nodes) # RNaseH degradation rates
            self.kd_A = 0.0003*np.ones(self.ortho_nodes) # RNaseA degradation rates
            self.kga = 1e4*np.ones(self.ortho_nodes)/1e9 # activation rates
            self.kgar = 5e3*np.ones(self.ortho_nodes)/1e9 # repression rates
            self.kar = 1e4*np.ones(self.ortho_nodes)/1e9 # activator inhibition rates
            self.kgb = 1e4*np.ones(self.ortho_nodes)/1e9 # free blocking rates
            self.kgbc = 5e3*np.ones(self.ortho_nodes)/1e9 # coactivation rates
            self.kbc = 1e4*np.ones(self.ortho_nodes)/1e9 # blocker inhibition rates
            self.kgab = 5e3*np.ones(self.ortho_nodes)/1e9 # active blocking rates
            self.kir = 1e4*np.ones(self.ortho_nodes)/1e9 # inducer binding rates
            
            # RNase degradation rates for RNAs bound to spatiotemporal genelets
            self.kd_Hg = 0.0003*np.ones(self.ind_nodes)
            self.kd_Ag = 0.0003*np.ones(self.ind_nodes)      
            
        elif len(rate_constants) != 0:
            if self.genelet_type.lower() == 'kwg':
                # defining the rate constants
                if isinstance(rate_constants[0],list) and len(rate_constants[0])>1:
                    self.kpr = np.array(rate_constants[0])
                else:
                    self.kpr = rate_constants[0]*np.ones(self.ind_nodes) # repressor production rates
                if isinstance(rate_constants[1],list) and len(rate_constants[1])>1:
                    self.kpc = np.array(rate_constants[1])
                else:
                    self.kpc = rate_constants[1]*np.ones(self.ind_nodes) # coactivator production rates
                if isinstance(rate_constants[2],list) and len(rate_constants[2])>1:
                    self.kpi = np.array(rate_constants[2])
                else:
                    self.kpi = rate_constants[2]*np.ones(self.ind_nodes) # inducer production rates
                if isinstance(rate_constants[3],list) and len(rate_constants[3])>1:
                    self.kd_H = np.array(rate_constants[3])
                else:
                    self.kd_H = rate_constants[3]*np.ones(self.ortho_nodes) # RNaseH degradation rates
                if isinstance(rate_constants[4],list) and len(rate_constants[4])>1:
                    self.kd_A = np.array(rate_constants[4])
                else:
                    self.kd_A = rate_constants[4]*np.ones(self.ortho_nodes) # RNaseA degradation rates
                if isinstance(rate_constants[5],list) and len(rate_constants[5])>1:
                    self.kga = np.array(rate_constants[5])
                else:
                    self.kga = rate_constants[5]*np.ones(self.ortho_nodes)/1e9 # activation rates
                if isinstance(rate_constants[6],list) and len(rate_constants[6])>1:
                    self.kgar = np.array(rate_constants[6])
                else:
                    self.kgar = rate_constants[6]*np.ones(self.ortho_nodes)/1e9 # repression rates
                if isinstance(rate_constants[7],list) and len(rate_constants[7])>1:
                    self.kar = np.array(rate_constants[7])
                else:
                    self.kar = rate_constants[7]*np.ones(self.ortho_nodes)/1e9 # activator inhibition rates
                if isinstance(rate_constants[8],list) and len(rate_constants[8])>1:
                    self.kgb = np.array(rate_constants[8])
                else:
                    self.kgb = rate_constants[8]*np.ones(self.ortho_nodes)/1e9 # free blocking rates
                if isinstance(rate_constants[9],list) and len(rate_constants[9])>1:
                    self.kgbc = np.array(rate_constants[9])
                else:
                    self.kgbc = rate_constants[9]*np.ones(self.ortho_nodes)/1e9 # coactivation rates
                if isinstance(rate_constants[10],list) and len(rate_constants[10])>1:
                    self.kbc = np.array(rate_constants[10])
                else:
                    self.kbc = rate_constants[10]*np.ones(self.ortho_nodes)/1e9 # blocker inhibition rates
                if isinstance(rate_constants[11],list) and len(rate_constants[11])>1:
                    self.kgab = np.array(rate_constants[11])
                else:
                    self.kgab = rate_constants[11]*np.ones(self.ortho_nodes)/1e9 # active blocking rates
                if isinstance(rate_constants[12],list) and len(rate_constants[12])>1:
                    self.kir = np.array(rate_constants[12])
                else:
                    self.kir = rate_constants[12]*np.ones(self.ortho_nodes)/1e9 # inducer binding rates
                
            elif self.genelet_type.lower() == 'stg':
                # defining the rate constants
                if isinstance(rate_constants[0],list) and len(rate_constants[0])>1:
                    self.kpr = np.array(rate_constants[0])
                else:
                    self.kpr = rate_constants[0]*np.ones(self.ind_nodes) # repressor production rates
                if isinstance(rate_constants[1],list) and len(rate_constants[1])>1:
                    self.kpi = np.array(rate_constants[1])
                else:
                    self.kpi = rate_constants[1]*np.ones(self.ind_nodes) # inducer production rates
                if isinstance(rate_constants[2],list) and len(rate_constants[2])>1:
                    self.kd_H = np.array(rate_constants[2])
                else:
                    self.kd_H = rate_constants[2]*np.ones(self.ortho_nodes) # RNaseH degradation rates
                if isinstance(rate_constants[3],list) and len(rate_constants[3])>1:
                    self.kd_A = np.array(rate_constants[3])
                else:
                    self.kd_A = rate_constants[3]*np.ones(self.ortho_nodes) # RNaseA degradation rates
                if isinstance(rate_constants[6],list) and len(rate_constants[6])>1:
                    self.kgar = np.array(rate_constants[6])
                else:
                    self.kgar = rate_constants[6]*np.ones(self.ortho_nodes)/1e9 # repression rates
                if isinstance(rate_constants[7],list) and len(rate_constants[7])>1:
                    self.kar = np.array(rate_constants[7])
                else:
                    self.kar = rate_constants[7]*np.ones(self.ortho_nodes)/1e9 # activator inhibition rates
                if isinstance(rate_constants[8],list) and len(rate_constants[8])>1:
                    self.kir = np.array(rate_constants[8])
                else:
                    self.kir = rate_constants[8]*np.ones(self.ortho_nodes)/1e9 # inducer binding rates
                
                # RNase degradation rates for RNAs bound to spatiotemporal genelets
                if isinstance(rate_constants[4],list) and len(rate_constants[4])>1:
                    self.kd_Ag = np.array(rate_constants[4])
                else:
                    self.kd_Ag = rate_constants[4]*np.ones(self.ind_nodes)      
                if isinstance(rate_constants[5],list) and len(rate_constants[5])>1:
                    self.kd_Hg = np.array(rate_constants[5])
                else:
                    self.kd_Hg = rate_constants[5]*np.ones(self.ind_nodes)
                
        k = self.ortho_nodes
        g = self.ind_nodes
        
        if iteration == 1:
            
            # calling the ODE equations 
            # here this is operating only on the initial conditions
            if self.genelet_type.lower() == 'kwg':
                # using the rate equations for the Kim and Winfree genelets
                rates = lambda t, x: general_genelet_eqs(t,x, \
                                self.ortho_nodes,self.ind_nodes,self.dA_tot,self.dB_tot, \
                                self.G_tot,self.kpr,self.kpc,self.kpi,self.kd_H,self.kd_A, \
                                self.kga,self.kgar,self.kar,self.kgb,self.kgab,self.kgbc, \
                                self.kbc,self.kir,self.act_mat,self.rep_mat,self.blk_mat, \
                                self.ca_mat,self.Cprod_mat,self.Rprod_mat,self.Cindc_mat,self.Rindc_mat, \
                                self.RnH,self.RnA,self.leak)
                
                # creating a new instance of the dA_tot and dB_tot variables for any additional sim steps
                self.dA_tot_n = copy.deepcopy(self.dA_tot)
                self.dB_tot_n = copy.deepcopy(self.dB_tot)
                                
            elif self.genelet_type.lower() == 'stg':
                # using the rate equations for the spatiotemporal genelets
                rates = lambda t, x: spatial_genelet_eqs(t,x, \
                        self.ortho_nodes,self.ind_nodes,self.dA_tot,self.G_tot, \
                        self.kpr,self.kpi,self.kd_H,self.kd_A,self.kd_Hg,self.kd_Ag, \
                        self.kgar,self.kar,self.kir,self.act_mat,self.rep_mat, \
                        self.Rprod_mat,self.Rindc_mat,self.RnH,self.RnA,self.leak)
                
                # creating a new instance of the dA_tot and dB_tot variables for any additional sim steps
                self.dA_tot_n = copy.deepcopy(self.dA_tot)
                
            # time interval for simulation
            t = t_vec # seconds
            
            # solving the ODEs
            self.sol = spi.solve_ivp(rates,[t[0],t[-1]],self.int_cond,t_eval=t)
            
            
        elif iteration > 1:
            
            # time interval for simulation
            t_n = t_vec # seconds
            
            # the initial conditions are now the last values of the previous
            int_cond_n = self.sol.y[:,-1]
            
            # calling the ODE equations
            # here this is operating on any updated conditions wrt dA_tot and dB_tot
            if self.genelet_type.lower() == 'kwg':
                # using the rate equations and updating conditions for the Kim and Winfree genelets
                rates_n = lambda t, x: general_genelet_eqs(t,x, \
                                self.ortho_nodes,self.ind_nodes,self.dA_tot_n,self.dB_tot_n, \
                                self.G_tot,self.kpr,self.kpc,self.kpi,self.kd_H,self.kd_A, \
                                self.kga,self.kgar,self.kar,self.kgb,self.kgab,self.kgbc, \
                                self.kbc,self.kir,self.act_mat,self.rep_mat,self.blk_mat, \
                                self.ca_mat,self.Cprod_mat,self.Rprod_mat,self.Cindc_mat,self.Rindc_mat, \
                                self.RnH,self.RnA,self.leak)
                
                # defining individual species for convenience below
                rR_n = int_cond_n[0:k] # repressors
                dA_n = int_cond_n[k:2*k] # activators
                rC_n = int_cond_n[2*k:3*k] # coactivators
                dB_n = int_cond_n[3*k:4*k] # blockers
                
                GdA_n = int_cond_n[4*k:4*k+g] # genelet:activator complexes
                GdB_n = int_cond_n[4*k+g:4*k+2*g] # genelet:blocker complexes
                
                rIr_n = int_cond_n[4*k+2*g:4*k+2*g+k] # repressor inducers 
                rIc_n = int_cond_n[4*k+2*g+k:4*k+2*g+2*k] # coactivator inducers
                
                dR_n = int_cond_n[4*k+2*g+2*k:4*k+2*g+3*k] # DNA repressors
                dAdR_n = int_cond_n[4*k+2*g+3*k:4*k+2*g+4*k] # DNA repressor:DNA activator complexes
                                
            elif self.genelet_type.lower() == 'stg':
                # using the rate equations for the spatiotemporal genelets
                rates_n = lambda t, x: spatial_genelet_eqs(t,x, \
                        self.ortho_nodes,self.ind_nodes,self.dA_tot_n,self.G_tot, \
                        self.kpr,self.kpi,self.kd_H,self.kd_A,self.kd_Hg,self.kd_Ag, \
                        self.kgar,self.kar,self.kir,self.act_mat,self.rep_mat, \
                        self.Rprod_mat,self.Rindc_mat,self.RnH,self.RnA,self.leak)
                
                # defining individual species for convenience below
                rR_n = int_cond_n[0:k] # repressors
                dA_n = int_cond_n[k:2*k] # activators
                
                Gon_n = int_cond_n[2*k:2*k+g] # genelet:activator complexes
                
                rIr_n = int_cond_n[2*k+g:2*k+g+k] # repressor inducers 
            
            # Update any new initial conditions based on non-default inputs
            # use multiple ifs here so more than one condition can be updated
            if rR != []:
                for i in range(len(rR_n)):
                    if isinstance(rR[i],str) == False:
                        rR_n[i] = rR_n[i] + rR[i]
                    
            if dA != []:
               for i in range(len(dA_n)):
                    if isinstance(dA[i],str) == False:
                        dA_n[i] = dA_n[i] + dA[i]
                        # need to update the total dA for mass balance
                        self.dA_tot_n[i] = self.dA_tot_n[i] + dA[i]
                        
            if rC != []:
                for i in range(len(rC_n)):
                    if isinstance(rC[i],str) == False:
                        rC_n[i] = rC_n[i] + rC[i]
                
            if dB != []:
                for i in range(len(dB_n)):
                    if isinstance(dB[i],str) == False:
                        dB_n[i] = dB_n[i] + dB[i]
                        # need to update the total dA for mass balance
                        self.dB_tot_n[i] = self.dB_tot_n[i] + dB[i]

            if rIr != []:
                for i in range(len(rIr_n)):
                    if isinstance(rIr[i],str) == False:
                        rIr_n[i] = rIr_n[i] + rIr[i]
            
            if rIc != []:
                for i in range(len(rIc_n)):
                    if isinstance(rIc[i],str) == False:
                        rIc_n[i] = rIc_n[i] + rIc[i]
                        
            if dR != []:
                for i in range(len(dR_n)):
                    if isinstance(dR[i],str) == False:
                        dR_n[i] = dR_n[i] + dR[i]
                
            if self.genelet_type.lower() == 'kwg':
                int_cond_n = np.concatenate([rR_n,dA_n,rC_n,dB_n,GdA_n,GdB_n,rIr_n,rIc_n,dR_n,dAdR_n])
            elif self.genelet_type.lower() == 'stg':
                int_cond_n = np.concatenate([rR_n,dA_n,Gon_n,rIr_n])
            
            # solving the odes
            sol_n = spi.solve_ivp(rates_n,[t_n[0],t_n[-1]],int_cond_n,t_eval=t_n)
            
            # appending each nth iteration to the previous solutions
            self.sol.t = np.concatenate([self.sol.t,sol_n.t])
            self.sol.y = np.concatenate([self.sol.y,sol_n.y],axis=1)
    
        # Saving a dictionary of all solutions for easy parsing
        # saving all concentrations for convience in assembling dictionary
        out_cons = self.sol.y
        
        self.output_concentration = {}
        
        if self.genelet_type.lower() == 'kwg':
            for i in range(k):
                # repressors
                self.output_concentration['rR'+str(i+1)] = out_cons[i]
                # activators
                self.output_concentration['dA'+str(i+1)] = out_cons[k+i]
                # coactivators
                self.output_concentration['rC'+str(i+1)] = out_cons[2*k+i]
                # blockers
                self.output_concentration['dB'+str(i+1)] = out_cons[3*k+i]
                # inducers for repressors
                self.output_concentration['rIr'+str(i+1)] = out_cons[4*k+2*g+i]
                # inducers for coactivators
                self.output_concentration['rIc'+str(i+1)] = out_cons[4*k+2*g+k+i]
                # DNA repressors
                self.output_concentration['dR'+str(i+1)] = out_cons[4*k+2*g+2*k+i]
                 # DNA repressor:DNA activator complexes
                self.output_concentration['dAdR'+str(i+1)] = out_cons[4*k+2*g+3*k+i]
                
            for i in range(g):
                # genelet:activator complexes
                self.output_concentration['GdA'+str(i+1)] = out_cons[4*k+i]
                # genelet:blocker complexes
                self.output_concentration['GdB'+str(i+1)] = out_cons[4*k+g+i]
                
        elif self.genelet_type.lower() == 'stg':
            for i in range(k):
                # repressors
                self.output_concentration['rR'+str(i+1)] = out_cons[i]
                # activators
                self.output_concentration['dA'+str(i+1)] = out_cons[k+i]
                # inducers for repressors
                self.output_concentration['rIr'+str(i+1)] = out_cons[2*k+g+i]
                
            for i in range(g):
                # genelet:activator complexes
                self.output_concentration['Gon'+str(i+1)] = out_cons[2*k+i]

    ''' PRINTING THE SYMBOLIC EQUATIONS FOR THE MODEL '''
    def export_sym_eqs(self,file_name = 'sym_eqs',path_name = ''):
        
        o_n = self.ortho_nodes + 1
        i_n = self.ind_nodes + 1
        
        RnH = self.RnH
        RnA = self.RnA
        leak = self.leak
        
        if isinstance(leak,list) and len(leak) == 1:
            leak = leak*self.ind_nodes
            leak = np.array([[leak[i]] for i in range(self.ind_nodes)])
            #sleak = [str(e) for e in leak]
            #leak = Matrix(symbols(sleak))
        elif isinstance(leak,list)!=True:
            leak = [leak]*self.ind_nodes
            leak = np.array([[leak[i]] for i in range(self.ind_nodes)])
            #sleak = [str(e) for e in leak]
            #leak = Matrix(symbols(sleak))
        else:
            leak = np.array([[leak[i]] for i in range(self.ind_nodes)])
            #sleak = [str(e) for e in leak]
            #leak = Matrix(symbols(sleak))
        
        if self.genelet_type.lower() == 'stg':
      
            # defining specie symbols
            rR = Matrix(symbols('rR1:%d'%o_n)) # repressors
            dA = Matrix(symbols('dA1:%d'%o_n)) # activators
            Gon = Matrix(symbols('Gon1:%d'%i_n)) # active genelets complexes
            rIr = Matrix(symbols('rIr1:%d'%o_n)) # repressor inducers 
            dA_tot = Matrix(symbols('dA_tot1:%d'%o_n)) # dA total
            G_tot = Matrix(symbols('G_tot1:%d'%i_n)) # dA total
            
            # rate constant symbols
            kp, kd_H, kd_A, kss, kgcs = symbols('kp kd_H kd_A kss kgcs')
            
            rep_mat = Matrix(self.rep_mat)
            Rprod_mat = Matrix(self.Rprod_mat)
            Rindc_mat = Matrix(self.Rindc_mat)
            
            # deleting repressor or inducer species that are not produced from equations
            for i in range(self.ortho_nodes):
                if sum(self.Rprod_mat[i,:])==0:
                    rR[i] = 0;
                if sum(self.Rindc_mat[i,:])==0:
                    rIr[i] = 0;                    
            
            # mass balances to find other concentrations
            dArR = dA_tot - dA # free activator:repressor complexes
            GrR = G_tot - Gon # OFF genelets
            
            # Rate equations:
            
            # repressors
            dRdt = -kgcs*MME((Matrix(rep_mat*Gon)),rR) - kss*MME(dA,rR) + kp*(Matrix(Rprod_mat*Gon)) + kp*Matrix(Rprod_mat*np.multiply(leak,GrR)) - kss*MME(rR,rIr) - RnA*kd_A*rR
            
            # activators
            dAdt = RnH*kd_H*dArR + RnA*kd_A*dArR - kss*MME(dA,rR)
            
            # G (ON)
            dGondt = RnH*kd_H*GrR + RnA*kd_A*GrR - kgcs*MME(Matrix(Transpose(rep_mat)*rR),Gon)
            
            # repressor inducers
            dIrdt = kp*(Matrix(Rindc_mat*Gon)) + kp*(Matrix(Rindc_mat*np.multiply(leak,GrR))) - kss*MME(rR,rIr) - RnA*kd_A*rIr
            
            f = open(path_name+file_name+'.txt','w')
            
            f.write("Rate Equations\n")
            for i in range(self.ind_nodes):
                f.write(" Gon%d' = " %(i+1) + "%s\n" %str(dGondt[i]))
            for i in range(self.ortho_nodes):
                f.write(" rR%d' = " %(i+1) + "%s\n" %str(dRdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" dA%d' = " %(i+1) + "%s\n" %str(dAdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" rIr%d' = " %(i+1) + "%s\n" %str(dIrdt[i]))
                
            f.write("Rate Constants\n")    
            f.write(" kp = " + "%s" %str(self.kpr[0]) + " | 1/s\n")
            f.write(" kd_H = " + "%s" %str(self.kd_H[0]) + " | 1/s\n")
            f.write(" kd_A = " + "%s" %str(self.kd_A[0]) + " | 1/s\n")
            f.write(" kss = " + "%s" %str(self.kga[0]*1e9) + " | 1/M-s\n")
            f.write(" kgcs = " + "%s" %str(self.kgar[0]*1e9) + " | 1/M-s\n")
            
            f.close()
            
        elif self.genelet_type.lower() == 'kwg':
            
            # defining specie symbols
            rR = Matrix(symbols('rR1:%d'%o_n)) # repressors
            rC = Matrix(symbols('rC1:%d'%o_n)) # repressors
            dA = Matrix(symbols('dA1:%d'%o_n)) # activators
            dB = Matrix(symbols('dB1:%d'%o_n)) # blockers
            GdA = Matrix(symbols('GdA1:%d'%i_n)) # active genelets complexes
            GdB = Matrix(symbols('GdB1:%d'%i_n)) # blocked genelets complexes
            rIr = Matrix(symbols('rIr1:%d'%o_n)) # repressor inducers 
            rIc = Matrix(symbols('rIc1:%d'%o_n)) # coactivator inducers 
            dA_tot = Matrix(symbols('dA_tot1:%d'%o_n)) # dA total
            dB_tot = Matrix(symbols('dB_tot1:%d'%o_n)) # dB total
            G_tot = Matrix(symbols('G_tot1:%d'%i_n)) # G total
            
            # rate constant symbols
            kp, kd_H, kd_A, kss, kgcs  = symbols('kp kd_H kd_A kss kgcs')
            
            act_mat = Matrix(self.act_mat)
            blk_mat = Matrix(self.blk_mat)
            ca_mat = Matrix(self.ca_mat)
            rep_mat = Matrix(self.rep_mat)
            Rprod_mat = Matrix(self.Rprod_mat)
            Cprod_mat = Matrix(self.Cprod_mat)
            Rindc_mat = Matrix(self.Rindc_mat)
            Cindc_mat = Matrix(self.Cindc_mat)
            
            # deleting repressor/inducer/blocker species that are not present equations
            for i in range(self.ortho_nodes):
                if sum(self.Rprod_mat[i,:])==0:
                    rR[i] = 0
                if sum(self.Rindc_mat[i,:])==0:
                    rIr[i] = 0
                if sum(self.Cindc_mat[i,:])==0:
                    rIc[i] = 0
                if sum(self.blk_mat[i,:])==0:
                    dB[i] = 0 
                    dB_tot[i] = 0
            # deleting GdB species that do not exist in equations
            for i in range(self.ind_nodes):
                if sum(self.blk_mat[:,i])==0:
                    GdB[i] = 0;
            
            # mass balances to find other concentrations
            dArR = dA_tot - dA - Matrix(act_mat*GdA) # activator:repressor complexes
            dBrC = dB_tot - dB - Matrix(blk_mat*GdB) # blocker:coactivator complexes
            G = G_tot - GdA - GdB # OFF genelets
            
            # Rate equations:
            
            # repressors
            dRdt = -kgcs*MME(Matrix(rep_mat*GdA),rR) - kss*MME(dA,rR) + kp*Matrix(Rprod_mat*GdA) + kp*Matrix(Rprod_mat*np.multiply(leak,GdB)) - kss*MME(rR,rIr) - RnA*kd_A*rR
            
            # activators
            dAdt = RnH*kd_H*dArR+ RnA*kd_A*dArR - kss*MME(Matrix(act_mat*G),dA) - kss*MME(dA,rR) + kgcs*MME(Matrix(blk_mat*GdA),dB)
            
            # coactivators
            dCdt = -kgcs*MME(Matrix(ca_mat*GdB),rC) - kss*MME(dB,rC) + kp*Matrix(Cprod_mat*GdA) + kp*Matrix(Cprod_mat*np.multiply(leak,GdB)) - kss*MME(rC,rIc) - RnA*kd_A*rC
            
            # blockers
            dBdt = RnH*kd_H*dBrC + RnA*kd_A*dBrC - kss*MME(Matrix(blk_mat*G),dB) - kss*MME(dB,rC) - kgcs*MME(Matrix(blk_mat*GdA),dB)
            
            # G (ON)
            dGondt = kss*MME(Matrix(Transpose(act_mat)*dA),G) - kgcs*MME(Matrix(Transpose(rep_mat)*rR),GdA) - kgcs*MME(Matrix(Transpose(blk_mat)*dB),GdA)
            
            #G (BLK)
            dGblkdt = kss*MME(Matrix(Transpose(blk_mat)*(dB)),G) + kgcs*MME(Matrix(Transpose(blk_mat)*dB),GdA) - kgcs*MME(Matrix(Transpose(ca_mat)*rC),GdB)
            
            # repressor inducers
            dIrdt = kp*Matrix(Rindc_mat*GdA) + kp*Matrix(Rindc_mat*np.multiply(leak,GdB)) - kss*MME(rR,rIr) - RnA*kd_A*rIr
            
            # coactivator inducers
            dIcdt = kp*Matrix(Cindc_mat*GdA) + kp*Matrix(Cindc_mat*np.multiply(leak,GdB)) - kss*MME(rC,rIc) - RnA*kd_A*rIc
    
            f = open(path_name+file_name+'.txt','w')
            
            f.write("Rate Equations\n")  
            for i in range(self.ind_nodes):
                f.write(" GdA%d' = " %(i+1) + "%s\n" %str(dGondt[i]))
            for i in range(self.ind_nodes):
                f.write(" GdB%d' = " %(i+1) + "%s\n" %str(dGblkdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" rR%d' = " %(i+1) + "%s\n" %str(dRdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" rC%d' = " %(i+1) + "%s\n" %str(dCdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" dA%d' = " %(i+1) + "%s\n" %str(dAdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" dB%d' = " %(i+1) + "%s\n" %str(dBdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" rIr%d' = " %(i+1) + "%s\n" %str(dIrdt[i]))
            for i in range(self.ortho_nodes):
                f.write(" rIc%d' = " %(i+1) + "%s\n" %str(dIcdt[i]))
                
            f.write("Rate Constants\n")    
            f.write(" kp = " + "%s" %str(self.kpr[0]) + " | 1/s\n")
            f.write(" kd_H = " + "%s" %str(self.kd_H[0]) + " | 1/s\n")
            f.write(" kd_A = " + "%s" %str(self.kd_A[0]) + " | 1/s\n")
            f.write(" kss = " + "%s" %str(self.kga[0]*1e9) + " | 1/M-s\n")
            f.write(" kgcs = " + "%s" %str(self.kgar[0]*1e9) + " | 1/M-s\n")
            
            f.close()
            
           
    ''' GENERATING THE NETWORK TOPOLOGY PLOT '''
    def plot_topology(self,layout='spring',pos=[],plot_title='',show_rnas=1): 
    
        top_mat = self.topology_mat
        Itop_mat = self.I_topology_mat
        ortho_n = self.ortho_nodes
        try:
            G_int_vec = self.G_int_vec
            default_int_states = 0
        except AttributeError:
            default_int_states = 1
        
        # counting total number of unique repressors and coactivators
        np_top_mat = np.array(top_mat)
        np_Itop_mat = np.array(Itop_mat)
        rna_ct = 0
        for i in range(ortho_n):
            if any(np_top_mat[i]==1):
                rna_ct += 1
            if any(np_top_mat[i]==-1):
                rna_ct += 1
            if any(np_Itop_mat[i]!=0) and all(np_top_mat[i]==0):
                rna_ct += 1
        
        NL = 'TN' # label for the nodes
        net_edges = []
        to_C_edges = []
        to_R_edges = []
        from_R_edges = []
        from_C_edges = []
        Cnode_edges = []
        aCnode_edges = []
        Rnode_edges = []
        aRnode_edges = []
        rna_nodes = []
        ind_edges = []
        ind_nodes = []
        
        if show_rnas == 1:
            n_color_val = [[1,1,1]]*(ortho_n+rna_ct)
            nsz = [1000]*(ortho_n+rna_ct)
        elif show_rnas == 0:
            n_color_val = [[1,1,1]]*(ortho_n)
            nsz = [1000]*(ortho_n)
        
        node_pos = []
        
        for i in range(ortho_n):
            for j in range(ortho_n):
            
                tmv = top_mat[i][j]
                Itmv = Itop_mat[i][j]
                
                # For topology matrix
                if tmv == 1:
                   
                    if i == j and show_rnas == 0:
                        # dummy nodes for autoactivation
                        aCnode_edges.append(['rC'+str(j+1),NL+str(i+1)])
                        nsz.append(1000)
                        n_color_val.append([1,1,1])
                    else:
                        Cnode_edges.append([NL+str(j+1),NL+str(i+1)])
                    
                    net_edges.append([NL+str(j+1),'rC'+str(i+1)])
                    to_C_edges.append([NL+str(j+1),'rC'+str(i+1)])
                  
                    net_edges.append(['rC'+str(i+1),NL+str(i+1)])
                    from_C_edges.append(['rC'+str(i+1),NL+str(i+1)])
                  
                if tmv == -1:
                    
                    if i == j and show_rnas == 0:
                        # dummy nodes for autorepression
                        aRnode_edges.append(['rR'+str(j+1),NL+str(i+1)])
                        nsz.append(1000)
                        n_color_val.append([1,1,1])
                    else:
                        Rnode_edges.append([NL+str(j+1),NL+str(i+1)])
                        
                    net_edges.append([NL+str(j+1),'rR'+str(i+1)])
                    to_R_edges.append([NL+str(j+1),'rR'+str(i+1)])
                    
                    net_edges.append(['rR'+str(i+1),NL+str(i+1)])
                    from_R_edges.append(['rR'+str(i+1),NL+str(i+1)])
                    
                # For the inducer topology matrix    
                if Itmv == 1:
                        
                    net_edges.append([NL+str(j+1),'rC'+str(i+1)])
                    ind_edges.append([NL+str(j+1),'rC'+str(i+1)])
                    ind_nodes.append([NL+str(j+1), NL+str(j+1)])
                    
                if Itmv == -1:
                    
                    net_edges.append([NL+str(j+1),'rR'+str(i+1)])
                    ind_edges.append([NL+str(j+1),'rR'+str(i+1)])
                    ind_nodes.append([NL+str(j+1), NL+str(j+1)])
                
        all_nodes = list(set(sum(net_edges,[])))
        self.net_edges = net_edges # setting attribute in case user wants to see
        
        G2 = nx.MultiDiGraph()

        if show_rnas == 1:
            G0 = nx.MultiDiGraph() # hold all the nodes without RNAs for default plot
    
            G2.add_edges_from(to_C_edges)
            G2.add_edges_from(from_C_edges)
            G2.add_edges_from(to_R_edges)
            G2.add_edges_from(from_R_edges)
            G2.add_edges_from(ind_edges)
            
            G0.add_edges_from(Cnode_edges)
            G0.add_edges_from(Rnode_edges)
            G0.add_edges_from(aCnode_edges)
            G0.add_edges_from(aRnode_edges)
            
        elif show_rnas == 0:
            G1 = nx.MultiDiGraph() # dummy for RNA nodes
            G3 = nx.MultiDiGraph() # dummy for TN nodes
            G2.add_edges_from(Cnode_edges)
            G2.add_edges_from(Rnode_edges)
            G2.add_edges_from(aCnode_edges)
            G2.add_edges_from(aRnode_edges)
            G3.add_edges_from(Cnode_edges)
            G3.add_edges_from(Rnode_edges)
            G1.add_edges_from(aCnode_edges)
            G1.add_edges_from(aRnode_edges)
            for i in range(len(ind_nodes)):
                G2.add_nodes_from(ind_nodes[i])
                G3.add_nodes_from(ind_nodes[i])
        
        # finding positions of TN nodes
        tni = 0        
        TN_nodes = []
        for d in G2.nodes():
            if d[0:2] == 'TN':
                node_pos.append(tni)
                TN_nodes.append(d)
            tni += 1
        
        # changing TN nodes to be bigger than RNA nodes                
        for x in node_pos:
            nsz[x] = 3000
        
        # reordering the nodes based on their position in the network layout
        node_pos_temp = copy.deepcopy(node_pos)
        node_order = [int(n[-1]) for n in TN_nodes] # actual order 
        for h in range(len(node_pos)):
            node_pos[node_order[h]-1] = node_pos_temp[h]
        
        # changing the colors of the TN nodes for the initial node states 
        if default_int_states == 1:
            for i in range(ortho_n):
                np_mat = np.array(top_mat[i])
                if sum(abs(np_mat))==0:
                    # node with no connections to it will be ON initally
                    n_color_val[node_pos[i]]=[0,0.8,0]
                elif any(np_mat==1):
                    # If a coactivator is produced for a node then it will be blocked initially
                    n_color_val[node_pos[i]]=[0.75,0.75,0.75]
                elif sum(np_mat)<0 and all(np_mat!=1):
                    # node that is only repressed will be ON initally
                    n_color_val[node_pos[i]]=[0,0.8,0]
        else:
            
            k = 0
            for i in range(ortho_n):
                np_mat = np.array(self.topology_mat[i])
                if G_int_vec[k] == 1:
                    # ON initially
                    n_color_val[node_pos[i]]=[0,0.8,0]
                    if sum(abs(np_top_mat[:,i])) == 0:
                        k = k + 1
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        k = k + 1
                elif G_int_vec[k] == -1:
                    # BLK initially
                    n_color_val[node_pos[i]]=[0.75,0.75,0.75]
                    if sum(abs(np_top_mat[:,i])) == 0:
                        k = k + 1
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        k = k + 1
                elif G_int_vec[k] == 0:
                    # OFF initially
                    n_color_val[node_pos[i]]=[1,0,0]
                    if sum(abs(np_top_mat[:,i])) == 0:
                        k = k + 1
                    for j in range(sum(abs(np_top_mat[:,i]))):
                        k = k + 1
                        
        # making the node sizes a bit bigger so there is room for the arrows
        pnsz = [e+1000 for e in nsz]
        
        if pos == []:
            
            # creating the layout based on only the TN nodes
            if layout.lower() == 'spring':
                if show_rnas == 1:
                    pos = nx.spring_layout(G0)
                elif show_rnas == 0:
                    pos = nx.spring_layout(TN_nodes)
            elif layout.lower() == 'spectral':
                if show_rnas == 1:
                    pos = nx.spectral_layout(G0)
                elif show_rnas == 0:
                    pos = nx.spectral_layout(TN_nodes)
            elif layout.lower() == 'circular':
                if show_rnas == 1:
                    pos = nx.circular_layout(G0)
                elif show_rnas == 0:
                    pos = nx.circular_layout(TN_nodes)
            elif layout.lower() == 'shell':
                if show_rnas == 1:
                    pos = nx.shell_layout(G0)
                elif show_rnas == 0:
                    pos = nx.shell_layout(TN_nodes)
            else:
                if show_rnas == 1:
                    pos = nx.spring_layout(G0)
                elif show_rnas == 0:
                    pos = nx.spring_layout(TN_nodes)
            
        if show_rnas == 1:
            for j in range(len(to_C_edges)):
                do_not_update = 0
                rC = to_C_edges[j][1]
                TNi = to_C_edges[j][0]
                TNi_1 = from_C_edges[j][1]
                # if the position has already been used don't append
                for v in pos.values():
                    # if a node self activates its RNA cannot be positioned between two nodes
                    if TNi == TNi_1:
                        do_not_update = -2
                    # otherwise check if there is a position available to place it between 2 nodes
                    elif v[0] == (pos[TNi][0]+pos[TNi_1][0])/2 and v[1] == (pos[TNi][1]+pos[TNi_1][1])/2:
                        if rC not in sum(to_C_edges[j+1:len(to_C_edges)],[]) or j+1 == len(to_C_edges):
                            # if there is not another place to put this node, offset its position a bit
                            do_not_update = -1
                            break
                        else: 
                            # don't append if there is another pair of nodes to put this between
                            do_not_update = 1
                            break
                        
                if do_not_update == 0 and rC not in pos.keys():        
                    pos.update({rC:[(pos[TNi][0]+pos[TNi_1][0])/2,(pos[TNi][1]+pos[TNi_1][1])/2]})
                elif do_not_update == -1 and rC not in pos.keys(): 
                    # if there is not another place to put this node, offset its position a bit
                    pos.update({rC:[(pos[TNi][0]+pos[TNi_1][0])/2,(pos[TNi][1]+pos[TNi_1][1])/2+0.25]})
                elif do_not_update == -2 and rC not in pos.keys(): 
                    # if a node self activates its RNA cannot be positioned between two nodes
                    pos.update({rC:[(pos[TNi][0]+pos[TNi_1][0])/2-0.5,(pos[TNi][1]+pos[TNi_1][1])/2+0.5]})
            
            for j in range(len(to_R_edges)):
                do_not_update = 0
                rR = to_R_edges[j][1]
                TNi = to_R_edges[j][0]
                TNi_1 = from_R_edges[j][1]
                # checking if the position has already been used
                for v in pos.values():
                     # if a node self activates its RNA cannot be positioned between two nodes
                    if TNi == TNi_1:
                        do_not_update = -2
                        
                    elif v[0] == (pos[TNi][0]+pos[TNi_1][0])/2 and v[1] == (pos[TNi][1]+pos[TNi_1][1])/2:
                        if rR not in sum(to_R_edges[j+1:len(to_R_edges)],[]) or j+1 == len(to_R_edges):
                            # if there is not another place to put this node, offset its position a bit
                            do_not_update = -1
                            break
                        else: 
                            # don't append if there is another pair of nodes to put this between
                            do_not_update = 1
                            break
                
                if do_not_update == 0 and rR not in pos.keys():
                    pos.update({rR:[(pos[TNi][0]+pos[TNi_1][0])/2,(pos[TNi][1]+pos[TNi_1][1])/2]})
                elif do_not_update == -1 and rR not in pos.keys(): 
                    # if there is not another place to put this node, offset its position a bit
                    pos.update({rR:[(pos[TNi][0]+pos[TNi_1][0])/2,(pos[TNi][1]+pos[TNi_1][1])/2+0.25]})
                elif do_not_update == -2 and rR not in pos.keys(): 
                    # if a node self activates its RNA cannot be positioned between two nodes
                    pos.update({rR:[pos[TNi][0]-0.5,pos[TNi][1]+0.5]})
                    
            for j in range(len(ind_edges)):
                rI = ind_edges[j][1]
                TNi = ind_edges[j][0]
                
                if TNi not in pos.keys():
                    if -1 < pos[rI][0] < 0.1 and -1 < pos[rI][1] < 0.1:
                            pos.update({TNi:[pos[rI][0]-0.25,pos[rI][1]-0.25]})
                    elif 0.1 < pos[rI][0] < 1.1 and 0.1 < pos[rI][1] < 1.1:
                            pos.update({TNi:[pos[rI][0]+0.25,pos[rI][1]+0.25]}) 
                    elif -1 < pos[rI][0] < 0.1 and 0.1 < pos[rI][1] < 1.1:
                            pos.update({TNi:[pos[rI][0]-0.25,pos[rI][1]+0.25]}) 
                    elif 0.1 < pos[rI][0] < 1.1 and -1 < pos[rI][1] < 0.1:
                            pos.update({TNi:[pos[rI][0]+0.25,pos[rI][1]-0.25]}) 
                                
            plt.figure()
            nx.draw_networkx_nodes(G2,pos=pos,node_color=n_color_val,node_size=nsz)
            nx.draw_networkx_labels(G2,pos=pos)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=to_C_edges,width=2,edge_color='k',arrowstyle='-|>',arrowsize=10,node_size=nsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=to_R_edges,width=2,edge_color='r',arrowstyle='-|>',node_size=nsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=from_C_edges,width=2,edge_color='k',arrowstyle='->',arrowsize=25,node_size=nsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=from_R_edges,width=2,edge_color='r',arrowstyle='|-|',arrowsize=10,node_size=pnsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=ind_edges,width=2,edge_color='r',arrowstyle='|-|',arrowsize=10,node_size=pnsz)
            left, right = plt.xlim()
            plt.xlim(left-0.25,right+0.25)
            down, up = plt.ylim()
            plt.ylim(down-0.25,up+0.25)
            if default_int_states == 1:
                plt.title(plot_title+' (Default Initial States)',fontsize=11,weight='bold')
            elif default_int_states == 0:
                plt.title(plot_title+' (User Initial States)',fontsize=11,weight='bold')
            
        elif show_rnas == 0:
            
            for j in range(len(aCnode_edges)):
                rC = aCnode_edges[j][0]
                TNi = aCnode_edges[j][1]
                pos.update({rC:[pos[TNi][0]-0.5,pos[TNi][1]+0.5]})
                
            for j in range(len(aRnode_edges)):
                rR = aRnode_edges[j][0]
                TNi = aRnode_edges[j][1]
                pos.update({rR:[pos[TNi][0]-0.5,pos[TNi][1]+0.5]})
            
            plt.figure()
            nx.draw_networkx_nodes(G2,pos=pos,node_color=n_color_val,node_size=nsz)
            #nx.draw_networkx_labels(G2,pos=pos)
            nx.draw_networkx_labels(G1,pos=pos,font_color='w')
            nx.draw_networkx_labels(G3,pos=pos,font_color='k')
            nx.draw_networkx_edges(G2,pos=pos,edgelist=Cnode_edges,width=2,edge_color='k',arrowstyle='->',arrowsize=25,node_size=pnsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=Rnode_edges,width=2,edge_color='r',arrowstyle='|-|',arrowsize=10,node_size=pnsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=aCnode_edges,width=2,edge_color='k',arrowstyle='-|>',arrowsize=15,node_size=pnsz)
            nx.draw_networkx_edges(G2,pos=pos,edgelist=aRnode_edges,width=2,edge_color='r',arrowstyle='-[',arrowsize=15,node_size=pnsz)
            left, right = plt.xlim()
            plt.xlim(left-0.25,right+0.25)
            down, up = plt.ylim()
            plt.ylim(down-0.25,up+0.25)
            if default_int_states == 1:
                plt.title(plot_title+' (Default Initial States)',fontsize=11,weight='bold')
            elif default_int_states == 0:
                plt.title(plot_title+' (User Initial States)',fontsize=11,weight='bold')
            
    ''' GENERATING THE GENELET SEQUENCES '''
    def compile_network_sequences(self,input_file_loc,desired_nodes=[],min_design=0,bth=1,save_file_name='genelet_seqs',save_path_name=''): 
        
        fun_key = xlsheet_to_dict(input_file_loc,'function_key',2)
        score_key = xlsheet_to_dict(input_file_loc,'scores',6)
        G_domains = xlsheet_to_dict(input_file_loc,'G_in_domains',1)
        activators = xlsheet_to_dict(input_file_loc,'activators',1)
        blockers = xlsheet_to_dict(input_file_loc,'blockers',1)
        R_nt = xlsheet_to_dict(input_file_loc,'R_nt_domains',1)
        C_nt = xlsheet_to_dict(input_file_loc,'C_nt_domains',1)
        
        ocw = {}
        orw = {}
        ocrw = {}
        icw = {}
        irw = {}
        icrw = {}
        # splitting each score into its own dictionary
        for k in score_key:
            ocw.update({k : score_key[k][0]})
            orw.update({k : score_key[k][1]})
            ocrw.update({k : score_key[k][2]})
            icw.update({k : score_key[k][3]})
            irw.update({k : score_key[k][4]})
            icrw.update({k : score_key[k][5]})
            
        # sorting all the domains by different scores
        sort_ocw = sorted(ocw.items(), key=lambda x: x[1], reverse=True)
        sort_orw = sorted(orw.items(), key=lambda x: x[1], reverse=True)
        sort_ocrw = sorted(ocrw.items(), key=lambda x: x[1], reverse=True)
        sort_icw = sorted(icw.items(), key=lambda x: x[1], reverse=True)
        sort_irw = sorted(irw.items(), key=lambda x: x[1], reverse=True)
        sort_icrw = sorted(icrw.items(), key=lambda x: x[1], reverse=True)
        
        topology_mat = self.topology_mat
        Itop_mat = self.I_topology_mat
        
        np_top_mat = np.array(topology_mat)
        np_Itop_mat = np.array(Itop_mat)
        ortho_nodes = self.ortho_nodes
        C_R_nodes = []
        R_nodes = []
        C_nodes = []
        any_nodes = []
        selected_G_domains = {}
        
        # determining what each topological node has to be able to do in the network
        for i in range(ortho_nodes):
            if sum(abs(np_top_mat[i,:]))==0:
                # node with no connections to it doesnt need blocker and can be any node
                any_nodes.append('TN'+str(i+1))
                
            elif any(np_top_mat[i,:]==1) and any(np_top_mat[i,:]==-1):
                # If a coactivator and repressor is produced for a node
                # then it will need to work in both directions
                # Needs blocker
                C_R_nodes.append('TN'+str(i+1))
                
            elif any(np_top_mat[i,:]==-1):
                # node that is only repressed will need to be repressible
                # Doesn't need a blocker
                R_nodes.append('TN'+str(i+1))
                
            elif any(np_top_mat[i,:]==1):
                # node that is onlycoactivated will need to be activatable
                # Needs blocker
                C_nodes.append('TN'+str(i+1))
           
        # picking nodes based on rankings because none were specified
        if desired_nodes == []:
            # if the user doesnt input the nodes to use then select some for them
            
            for i in range(len(C_R_nodes)):
                for m in range(len(sort_ocrw)):
                    k = sort_ocrw[m][0] # the key is now from the sorted list based on ocrw
                    if k not in selected_G_domains.values():
                        # find keys that can both be coact / rep (first one found is used)
                        if fun_key[k][0] == 'C1' and fun_key[k][1] == 'R1':
                            selected_G_domains.update({C_R_nodes[i] : k})
                            break
                    
            for i in range(len(any_nodes)):
                for m in range(len(sort_ocrw)):
                    k = sort_ocrw[m][0] # the key is now from the sorted list based on ocrw
                    if k not in selected_G_domains.values():
                        # find keys that can both be coact / rep (first one found is used)
                        if fun_key[k][0] == 'C1' and fun_key[k][1] == 'R1':
                            selected_G_domains.update({any_nodes[i] : k})
                            break
                        # then pick nodes that can only be coactivated
                        elif fun_key[k][0] == 'C1' and fun_key[k][1] == 'R0':
                            selected_G_domains.update({any_nodes[i] : k})
                            break  
                        # then last pick nodes that can only be repressed
                        elif fun_key[k][0] == 'C0' and fun_key[k][1] == 'R1':
                            selected_G_domains.update({any_nodes[i] : k})
                            break  
            
            c_check = 0
            for i in range(len(C_nodes)):
                # if doing minimal design then want to pick nodes ranked by just coactivation
                if min_design == 1 and c_check == 0:
                    for m in range(len(sort_ocw)):
                        k = sort_ocw[m][0] # the key is now from the sorted list based on ocw
                        # can't reuse input domains
                        if k not in selected_G_domains.values():
                            # find the first node that can only be coact
                            if fun_key[k][0] == 'C1' and fun_key[k][1] == 'R0':
                                selected_G_domains.update({C_nodes[i] : k})
                                break  
                            elif m == len(sort_ocw)-1: # if you run out of min design nodes need to pick from other nodes
                                c_check = 1 # exit condition to drop to the below else
                            
                # select based on combined coact/rep ranks if not doing min design
                # or if you run out of minimal domains
                if min_design == 0 or c_check == 1: 
                    
                    for m in range(len(sort_ocw)):
                        k = sort_ocw[m][0] # the key is now from the sorted list based on ocrw
                        # can't reuse input domains
                        if k not in selected_G_domains.values():
                            # find the first node that can at least be coact
                            if fun_key[k][0] == 'C1':
                                selected_G_domains.update({C_nodes[i] : k})
                                break  
            r_check = 0     
            for i in range(len(R_nodes)):
                # if doing minimal design then want to pick nodes ranked by just repression
                if min_design == 1 and r_check == 0:
                    for m in range(len(sort_orw)):
                        k = sort_orw[m][0] # the key is now from the sorted list based on orw
                        # can't reuse input domains
                        if k not in selected_G_domains.values():
                            # find the first node that can only be rep
                            if fun_key[k][0] == 'C0' and fun_key[k][1] == 'R1':
                                selected_G_domains.update({R_nodes[i] : k})
                                break  
                            elif m == len(sort_orw)-1: # if you run out of min design nodes need to pick from other nodes
                                r_check = 1 # exit condition to drop to the below condition
                
                # select based on combined coact/rep ranks if not doing min design
                # or if you run out of minimal domains
                if min_design == 0 or r_check == 1: 
                    
                    for m in range(len(sort_orw)):
                        k = sort_orw[m][0] # the key is now from the sorted list based on ocrw
                        # can't reuse input domains
                        if k not in selected_G_domains.values():
                            # find the first node that can at least be coact
                            if fun_key[k][1] == 'R1':
                                selected_G_domains.update({R_nodes[i] : k})
                                break 
        
       # using the nodes specified by the user                
        else:
            
            i = 0
                
            for k in desired_nodes:
                
                if k[0] == 'G':
                    if k not in selected_G_domains.values():
                        # check to make sure the selected node meets regulation criteria
                        if C_R_nodes.count('TN'+str(i+1)) > 0 and fun_key[k][0] == 'C1' and fun_key[k][1] == 'R1':
                            selected_G_domains.update({'TN'+str(i+1) : k})
                        elif C_nodes.count('TN'+str(i+1)) > 0 and fun_key[k][0] == 'C1':
                            selected_G_domains.update({'TN'+str(i+1) : k})
                        elif R_nodes.count('TN'+str(i+1)) > 0 and fun_key[k][1] == 'R1':
                            selected_G_domains.update({'TN'+str(i+1) : k})
                        elif any_nodes.count('TN'+str(i+1)) > 0:
                            selected_G_domains.update({'TN'+str(i+1) : k})
                        else: # none of the criteria were met so your selection is no good
                            sys.exit('The position for domain '+k+' is not valid')
                        i += 1
                   
                else: # if the specified node does not start with G (i.e. left blank), this node should be selected 
                    
                    # First figure out what the unspecified node needs to do
                    # then execute the algorithm for selecting domains
                    if C_R_nodes.count('TN'+str(i+1)) > 0:
                        for m in range(len(sort_ocrw)):
                            k = sort_ocrw[m][0] # the key is now from the sorted list based on ocrw
                            if k not in selected_G_domains.values():
                                # find keys that can both be coact / rep (first one found is used)
                                if fun_key[k][0] == 'C1' and fun_key[k][1] == 'R1':
                                    selected_G_domains.update({'TN'+str(i+1) : k})
                                    break
                                                    
                    elif any_nodes.count('TN'+str(i+1)) > 0:
                        for m in range(len(sort_ocrw)):
                            k = sort_ocrw[m][0] # the key is now from the sorted list based on ocrw
                            if k not in selected_G_domains.values():
                                # find keys that can both be coact / rep (first one found is used)
                                if fun_key[k][0] == 'C1' and fun_key[k][1] == 'R1':
                                    selected_G_domains.update({'TN'+str(i+1) : k})
                                    break
                                # then pick nodes that can only be coactivated
                                elif fun_key[k][0] == 'C1' and fun_key[k][1] == 'R0':
                                    selected_G_domains.update({'TN'+str(i+1) : k})
                                    break  
                                # then last pick nodes that can only be repressed
                                elif fun_key[k][0] == 'C0' and fun_key[k][1] == 'R1':
                                    selected_G_domains.update({'TN'+str(i+1) : k})
                                    break  
                        
                    elif C_nodes.count('TN'+str(i+1)) > 0:
                        c_check = 0
                        # if doing minimal design then want to pick nodes ranked by just coactivation
                        if min_design == 1 and c_check == 0:
                            for m in range(len(sort_ocw)):
                                k = sort_ocw[m][0] # the key is now from the sorted list based on ocw
                                # can't reuse input domains
                                if k not in selected_G_domains.values():
                                    # find the first node that can only be coact
                                    if fun_key[k][0] == 'C1' and fun_key[k][1] == 'R0':
                                        selected_G_domains.update({'TN'+str(i+1) : k})
                                        break  
                                    elif m == len(sort_ocw)-1: # if you run out of min design nodes need to pick from other nodes
                                        c_check = 1 # exit condition to drop to the below else
                                    
                        # select based on combined coact/rep ranks if not doing min design
                        # or if you run out of minimal domains
                        if min_design == 0 or c_check == 1: 
                            
                            for m in range(len(sort_ocw)):
                                k = sort_ocw[m][0] # the key is now from the sorted list based on ocrw
                                # can't reuse input domains
                                if k not in selected_G_domains.values():
                                    # find the first node that can at least be coact
                                    if fun_key[k][0] == 'C1':
                                        selected_G_domains.update({'TN'+str(i+1) : k})
                                        break  
                        
                    elif R_nodes.count('TN'+str(i+1)) > 0:
                        r_check = 0
                        # if doing minimal design then want to pick nodes ranked by just repression
                        if min_design == 1 and r_check == 0:
                            for m in range(len(sort_orw)):
                                k = sort_orw[m][0] # the key is now from the sorted list based on orw
                                # can't reuse input domains
                                if k not in selected_G_domains.values():
                                    # find the first node that can only be rep
                                    if fun_key[k][0] == 'C0' and fun_key[k][1] == 'R1':
                                        selected_G_domains.update({'TN'+str(i+1) : k})
                                        break  
                                    elif m == len(sort_orw)-1: # if you run out of min design nodes need to pick from other nodes
                                        r_check = 1 # exit condition to drop to the below condition
                        
                        # select based on combined coact/rep ranks if not doing min design
                        # or if you run out of minimal domains
                        if min_design == 0 or r_check == 1: 
                            
                            for m in range(len(sort_orw)):
                                k = sort_orw[m][0] # the key is now from the sorted list based on ocrw
                                # can't reuse input domains
                                if k not in selected_G_domains.values():
                                    # find the first node that can at least be coact
                                    if fun_key[k][1] == 'R1':
                                        selected_G_domains.update({'TN'+str(i+1) : k})
                                        break 
                        
                    else: # none of the criteria were met so your selection is no good
                        sys.exit('No valid domains could be identified for '+'TN'+str(i+1))
                    i += 1
        
        if len(selected_G_domains)<ortho_nodes:
            sys.exit('There are not enough nodes to create the topology with the given inputs')
            
        # now that we know what nodes go where lets stitch this shit together
        c = 0
        genelet_sequences = []
        genelet_seq_dict = {}
        rc_append = 'CGACTCACTATA'
        PSHP = 'GGGAGATTCGTCTCCCAA'
        for k in range(len(selected_G_domains)):
            c += 1
            gdk = selected_G_domains['TN'+str(k+1)] # key to the G_domain
            genelet_sequences = []
            for i in range(len(np_top_mat[:,1])):
                
                if np_top_mat[i,c-1] == 1:
                # this node coactivates another
                    
                    Cpos = selected_G_domains['TN'+str(i+1)]
                    if 'TN'+str(k+1) in R_nodes and bth == 0:
                    # removing the BTH sequence if desired by user for only repressed nodes
                        genelet_sequences.append(['*'+gdk+'C'+Cpos[1:]+'-nt: 5'+G_domains[gdk][0][8:]+C_nt['C'+Cpos[1:]+'-nt'][0]])
                    else:
                        genelet_sequences.append([gdk+'C'+Cpos[-1]+'-nt: 5'+G_domains[gdk][0]+C_nt['C'+Cpos[1:]+'-nt'][0]])
                        
                    genelet_sequences.append(['C'+Cpos[1:]+'-t: '+rev_comp_seq(rc_append+C_nt['C'+Cpos[1:]+'-nt'][0])[:-1]])
                    
                elif np_top_mat[i,c-1] == -1:
                # this node represses another
                    
                    Rpos = selected_G_domains['TN'+str(i+1)]
                    if 'TN'+str(k+1) in R_nodes and bth == 0:
                    # removing the BTH sequence if desired by user for only repressed nodes
                        genelet_sequences.append(['*'+gdk+'R'+Rpos[1:]+'-nt: 5'+G_domains[gdk][0][8:]+R_nt['R'+Rpos[1:]+'-nt'][0]])
                    else:
                        genelet_sequences.append([gdk+'R'+Rpos[1:]+'-nt: 5'+G_domains[gdk][0]+R_nt['R'+Rpos[1:]+'-nt'][0]])
                        
                    genelet_sequences.append(['R'+Rpos[1:]+'-t: '+rev_comp_seq(rc_append+R_nt['R'+Rpos[1:]+'-nt'][0])[:-1]])
                
                elif np_Itop_mat[i,c-1] == 1:
                # this node produces an inducer for another nodes repressor
                    
                    Cpos = selected_G_domains['TN'+str(i+1)]
                    Ck = 'C'+Cpos[1:]+'-nt'
                    if 'TN'+str(k+1) in R_nodes and bth == 0:
                    # removing the BTH sequence if desired by user for only repressed nodes
                        genelet_sequences.append(['*'+gdk+'Ic'+Cpos[1:]+'-nt: 5'+G_domains[gdk][0][8:]+PSHP+rev_comp_seq(C_nt[Ck][0][-16:])[1:-1]])
                    else:
                        genelet_sequences.append([gdk+'Ic'+Cpos[1:]+'-nt: 5'+G_domains[gdk][0]+PSHP+rev_comp_seq(C_nt['C'+Cpos[1:]+'-nt'][0][-16:])[1:-1]])
                        
                    genelet_sequences.append(['Ic'+Cpos[-1]+'-t: 5'+C_nt[Ck][0][-16:]+rev_comp_seq(rc_append+PSHP)[1:-1]])
                    
                elif np_Itop_mat[i,c-1] == -1:
                # this node produces an inducer for another nodes repressor
                    
                    Rpos = selected_G_domains['TN'+str(i+1)]
                    Rk = 'R'+Rpos[1:]+'-nt'
                    if 'TN'+str(k+1) in R_nodes and bth == 0:
                    # removing the BTH sequence if desired by user for only repressed nodes
                        genelet_sequences.append(['*'+gdk+'Ir'+Rpos[1:]+'-nt: 5'+G_domains[gdk][0][8:]+PSHP+rev_comp_seq(R_nt['R'+Rpos[1:]+'-nt'][0][-16:])[1:-1]])
                    else:
                        genelet_sequences.append([gdk+'Ir'+Rpos[1:]+'-nt: 5'+G_domains[gdk][0]+PSHP+rev_comp_seq(R_nt['R'+Rpos[1:]+'-nt'][0][-16:])[1:-1]])
                        
                    genelet_sequences.append(['Ir'+Rpos[1:]+'-t: 5'+R_nt[Rk][0][-16:]+rev_comp_seq(rc_append+PSHP)[1:-1]])
                    
                elif sum(abs(np_top_mat[:,c-1])) == 0 and sum(abs(np_Itop_mat[:,c-1])) == 0 and i == len(np_top_mat[:,1])-1:
                # this node doesn't produce an RNA that connects to anything 
                                
                    if 'TN'+str(k+1) in R_nodes and bth == 0:
                    # removing the BTH if desired by the user
                        genelet_sequences.append(['*'+gdk+'D-nt'+': 5'+G_domains[gdk][0][8:len(G_domains[gdk][0])]+'GGGAGA'])
                    else:
                        genelet_sequences.append([gdk+'D-nt'+': 5'+G_domains[gdk][0]+'GGGAGA'])
                   
                    genelet_sequences.append([gdk+'D-t: '+rev_comp_seq(rc_append+'GGGAGA')[:-1]])
                    
                # add the activators and blockers at the end    
                if i == len(np_top_mat[:,1])-1:
                    genelet_sequences.append(['dA'+gdk[1:]+': 5'+activators['dA'+gdk[1:]][0]])
                    if any(np_top_mat[c-1,:] == 1):
                        # only add blockers if this node has a coactivator coming to it
                        genelet_sequences.append(['dB'+gdk[1:]+': 5'+blockers['dB'+gdk[1:]][0]])
                
                genelet_seq_dict.update({'TN'+str(k+1) : genelet_sequences})
        
        # write to a text file
        f = open(save_path_name+save_file_name+'.txt','w')
                
        for k in genelet_seq_dict.keys():
            f.write('------'+k+'-----'+'\n')
            for v in range(len(genelet_seq_dict[k])):
                f.write(genelet_seq_dict[k][v][0]+'\n')
            
        f.close()
        
        
        
        