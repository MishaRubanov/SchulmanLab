# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 20:37:10 2020

@author: sscha
"""
import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM

'''
Examples you want to run
'''
ex1 = 1 # KWG: IFFL with 150 nM added blockers and how to use with topology mat input
ex2 = 0 # KWG: IFFL dA3 titration
ex3 = 0 # KWG: extended IFFL with 5% leak with and without added blockers
ex4 = 0 # KWG: TSN state 3 -> state 1 and topology plotting with default/user int conditions
ex5 = 0 # KWG: TSN state 2 -> state 1 -> state 3
ex6 = 0 # KWG: IFFL with 150 nM added blockers and user defined rate consants
ex7 = 0 # STG: TSN with RNase A
ex8 = 0 # STG: G1-|G2-|G3 cascade with 5% leak with and without added dA3
ex9 = 0 # STG: G1-|G3 and G2-|rR3 (subtraction circuit) topology/Itopology mat inputs as well 
ex10 = 0 # KWG: TSN with upstream inducer nodes state 1 -> 2 -> 3 
ex11 = 0 # KWG: BS_IFFL1|2 network
ex12 = 0 # KWG: I_BS_IFFL1|2 network
ex13 = 0 # KWG: FBI_BS_IFFL1|2 network (IFFL1 -> State 2)
ex14 = 0 # KWG: FBI_BS_IFFL1|2 network (induce to State2 and IFFL2 -> State 1)
ex15 = 0 # KWG: FBI_BS_IFFL1|2 network (fully connected oscillator)
ex16 = 0 # KWG: simpler induction / bistable / pulse feedback potential oscillator
ex17 = 0 # KWG: Misha's simple pulse network

'''
###############################################################################
###############################################################################
KIM AND WINFREE GENELET EXAMPLES
###############################################################################
###############################################################################
'''

''' 
EXAMPLE 1
###############################################################################
This simulates a IFFL with 150 nM excess blocker for G3 and G4
-------------------------------------------------------------------------------
   -This example primarily illustrates how to to simply use the GGM
###############################################################################
'''
if ex1 == 1:
    
    # OG Sam method of network definition
    #          G1 G2 G3 G4      
    act_vec =  [1, 1, 2, 3]
#    blk_vec =  [0, 0, 2, 3]
#    prod_vec = [3, 2,-3, 0]
#    indc_vec = [0, 0, 0, 0]
    
    # Misha's topology matrix method of network definition
    topology_mat = [[0,0,0],[1,0,0],[1,-1,0]]
    
    ''' initializing topology '''
    #IFFL = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
    IFFL = GGM.GeneletNetwork(top_mat=topology_mat)
    IFFL.plot_topology()
    IFFL.plot_topology(show_rnas=0)
    
    # Define initial conditions
    #                  dA1  dA2  dA3
    dA_tot = np.array([250, 250, 250]) # total activator added
    #                 G1  G2  G3  G4
    G_tot = np.array([5, 25, 25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4
    #G_int_vec = [0, 0,-1,-1]
    
    ''' initializing initial conditions '''
    IFFL.initial_conditions(dA_tot,G_tot,G_int_vec=[],dB_added = [0,150,150]) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds
    
    ''' simulating the IFFL'''
    IFFL.simulate(t_vec1,1)
    
    # pulling out the desired concentrations for plotting
    G1 = IFFL.output_concentration['GdA1']
    G3 = IFFL.output_concentration['GdA3']
    G4 = IFFL.output_concentration['GdA4']
    
    sim_t = IFFL.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('IFFL with added blockers',fontsize=fs+2,weight='bold')
    plt.subplot(3,4,1)
    plt.plot(sim_t/3600,G1/G_tot[0],color='g',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,5)
    plt.plot(sim_t/3600,G3/G_tot[2],color='r',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,9)
    plt.plot(sim_t/3600,G4/G_tot[3],color='b',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 2
###############################################################################
This simulates a IFFL with a titration of dA3 concentrations
-------------------------------------------------------------------------------
   -This example primarily illustrates how to to simply use the GGM
###############################################################################
'''
if ex2 == 1:
    #          G1 G2 G3 G4      
    act_vec =  [1, 1, 2, 3]
    blk_vec =  [0, 0, 2, 3]
    prod_vec = [3, 2,-3, 0]
    indc_vec = [0, 0, 0, 0]
    
    ''' initializing topology '''
    IFFLt = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
    IFFLt.plot_topology()
    IFFLt.plot_topology(show_rnas=0)
    
    # Define initial conditions
    dA_tot =[0]*3 # initializing dA list
    #                     dA1  dA2  dA3
    dA_tot[0] = np.array([250, 250, 125]) # total activator added
    dA_tot[1] = np.array([250, 250, 250]) # total activator added
    dA_tot[2] = np.array([250, 250, 500]) # total activator added
    
    #                 G1  G2  G3  G4
    G_tot = np.array([25, 5, 25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4
    G_int_vec = [0, 0,-1,-1]
    
    G4 = [0]*3 # initializing G4 for storage
    cv = [[0,0,0.5],[0,0,1],[0,0.5,1]] # colors for plotting
    fs = 13 # font size for plotting
    plt.figure()
    
    for i in range(len(dA_tot)):
    
        ''' initializing initial conditions '''
        IFFLt.initial_conditions(dA_tot[i],G_tot,G_int_vec) # default of 0 for all other initial conditions
        
        t_vec1 = np.linspace(0,3,1001)*3600 # seconds
        
        ''' simulating the IFFL'''
        IFFLt.simulate(t_vec1,1)
        
        # pulling out the desired concentrations for plotting
        G4[i] = IFFLt.output_concentration['GdA4']
        sim_t = IFFLt.sol.t
        
        plt.subplot(3,4,1)
        plt.plot(sim_t/3600,G4[i]/G_tot[3],color=cv[i],linewidth=4)
        
    plt.suptitle('IFFL with dA3 titration',fontsize=fs+2,weight='bold')    
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 3
###############################################################################
This simulates an extended IFFL where an upstream genelet produces the coactivator
of the first two nodes of the IFFL with and without 5% leak
This simulates the network both with and without 150 nM dB2,dB3,dB4 which suppress the leak
-------------------------------------------------------------------------------
   -This example illustrates how to include the leak reactions
   -This example illustrates how to export the KWG symbolic equations
###############################################################################
'''
if ex3 == 1:
    #          G1 G2 G3 G4 G5     
    act_vec =  [1, 2, 2, 3, 4]
    blk_vec =  [0, 2, 2, 3, 4]
    prod_vec = [2, 3, 4,-4, 0]
    indc_vec = [0, 0, 0, 0, 0]
    
#    top_mat = [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 1, -1, 0]]
#    eIFFL = GGM.GeneletNetwork(top_mat=top_mat)
    
    ''' initializing topology '''
    eIFFL = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
    eIFFL.plot_topology(layout='circular')
    eIFFL.plot_topology(layout='circular',show_rnas=0)
    
    # Define initial conditions
    dA_tot = [0]*2
    #                     dA1  dA2  dA3  dA4
    dA_tot[0] = np.array([250, 250, 250, 250]) # total activator added
    dA_tot[1] = np.array([  0, 250, 250, 250]) # total activator added
    #                 G1  G2  G3 G4  G5
    G_tot = np.array([25, 5, 25, 25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4 G5
    G_int_vec = [0, -1,-1,-1,-1]
    
    G2 = [0]*2
    G4 = [0]*2
    G5 = [0]*2
    ls = ['-',':']
    '''
    0 nM ADDED BLOCKERS
    '''
    plt.figure()
    for i in range(len(dA_tot)):
 
        ''' initializing initial conditions '''
        eIFFL.initial_conditions(dA_tot[i],G_tot,G_int_vec) # default of 0 for all other initial conditions
        
        t_vec1 = np.linspace(0,3,1001)*3600 # seconds
        
        ''' simulating the IFFL'''
        eIFFL.simulate(t_vec1,1,leak=0.05)
        
        # pulling out the desired concentrations for plotting
        
        G2[i] = eIFFL.output_concentration['GdA2']
        G4[i] = eIFFL.output_concentration['GdA4']
        G5[i] = eIFFL.output_concentration['GdA5']
        
        sim_t = eIFFL.sol.t
        
        fs = 13
        
        plt.suptitle('Extended IFFL with leak',fontsize=fs+2,weight='bold')
        plt.subplot(3,4,2)
        plt.title('No addded blockers')
        plt.plot(sim_t/3600,G2[i]/G_tot[1],color='g',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,6)
        plt.plot(sim_t/3600,G4[i]/G_tot[3],color='r',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,10)
        plt.plot(sim_t/3600,G5[i]/G_tot[4],color='b',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xlabel('time (hr)',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
    G2 = [0]*2
    G4 = [0]*2
    G5 = [0]*2
    '''
    150 nM ADDED BLOCKERS
    '''
    for i in range(len(dA_tot)):
        ''' initializing initial conditions '''
        eIFFL.initial_conditions(dA_tot[i],G_tot,G_int_vec,dB_added=[0,150,150,150]) # default of 0 for all other initial conditions
        
        t_vec1 = np.linspace(0,3,1001)*3600 # seconds
        
        ''' simulating the IFFL'''
        eIFFL.simulate(t_vec1,1,leak=0.05)
        
        # pulling out the desired concentrations for plotting
        G2[i] = eIFFL.output_concentration['GdA2']
        G4[i] = eIFFL.output_concentration['GdA4']
        G5[i] = eIFFL.output_concentration['GdA5']
        
        sim_t = eIFFL.sol.t
        
        fs = 13
        plt.suptitle('Extended IFFL with leak',fontsize=fs+2,weight='bold')
        plt.subplot(3,4,4)
        plt.title('150 nM blockers')
        plt.plot(sim_t/3600,G2[i]/G_tot[1],color='g',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,8)
        plt.plot(sim_t/3600,G4[i]/G_tot[3],color='r',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,12)
        plt.plot(sim_t/3600,G5[i]/G_tot[4],color='b',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xlabel('time (hr)',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

    ''' Exporting symbolic equations '''
    eIFFL.export_sym_eqs(file_name='kwg_testing_sym_export', \
                         path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\')


''' 
EXAMPLE 4
###############################################################################
This simulates setting the TSN into State 3 by initally including rR1 and rR2
and then switching to State 1 by addding rIr1 after 2.5 hours
-------------------------------------------------------------------------------
   -This example illustrates the use of multiple iterations of model.simulate
   -This example also illustrates how to add non-default initial conditions for rR
###############################################################################
'''
if ex4 == 1:
    
     # Sam's method for network definition
#    #          G1 G2 G3 G4 G5 G6      
#    act_vec =  [1, 1, 2, 2, 3, 3]
#    blk_vec =  [0, 0, 0, 0, 0, 0]
#    prod_vec = [-3,-2,-1,-3,-1,-2]
#    indc_vec = [0, 0, 0, 0, 0, 0]
    
    # Misha's topology matrix method of network definition
    topology_mat = [[0,-1,-1],[-1,0,-1],[-1,-1,0]]
    
    ''' initializing TSN topology '''
    #TSN = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
    TSN = GGM.GeneletNetwork(top_mat=topology_mat)
    TSN.plot_topology(plot_title='TSN State 3')
    
    # Define initial conditions
    #                  dA1  dA2  dA3
    dA_tot = np.array([250, 250, 250]) # total activator added
    #                 G1  G2  G3  G4  G5  G6
    G_tot = np.array([50, 50, 50, 50, 50, 50]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4 G5 G6
    G_int_vec = [0, 0, 0, 0, 1, 1]
    
    ''' initializing initial conditions '''
    # specifying an initial amount of rR1 and rR2 to set the system into State 3
    TSN.initial_conditions(dA_tot,G_tot,G_int_vec,rRin = [1000,1000,0]) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,2.5,1001)*3600 # seconds
    
    ''' simulating the TSN model for initial RNA production to set State 2'''
    TSN.simulate(t_vec1,1)
    
    t_vec2 = np.linspace(t_vec1[-1]/3600,6,1001)*3600 # seconds
    
    ''' simulating the TSN model after adding rIr1 to change to State 1 '''
    TSN.simulate(t_vec2,2,rIr=[15000,'NA','NA'])
    TSN.plot_topology(plot_title='TSN State 3')
    TSN.plot_topology(plot_title='TSN State 3',show_rnas=0)
    
    ''' Plotting solution '''
    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = TSN.output_concentration['GdA1']
    G3 = TSN.output_concentration['GdA3']
    G5 = TSN.output_concentration['GdA5']
    
    sim_t = TSN.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('TSN: State 3 to State 1',fontsize=fs+2,weight='bold')
    plt.subplot(3,3,1)
    plt.plot(sim_t/3600,G1/G_tot[0],color='b',linewidth=4)
    plt.plot(sim_t/3600,G3/G_tot[2],color='r',linewidth=4)
    plt.plot(sim_t/3600,G5/G_tot[4],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 5
###############################################################################
This simulates a double flip from State 2 to State 1 to State 3
-------------------------------------------------------------------------------
   -This example sets up the initial condtions like I did experimentally where
    initially only one activator is present and the others are added after 30 m
   -This example illustrates the use of multiple iterations of model.simulate
   -This example also illustrates pulling out more concentrations for plotting
###############################################################################
'''
if ex5 == 1:
    #          G1 G2 G3 G4 G5 G6      
    act_vec =  [1, 1, 2, 2, 3, 3]
    blk_vec =  [0, 0, 0, 0, 0, 0]
    prod_vec = [-3,-2,-1,-3,-1,-2]
    indc_vec = [0, 0, 0, 0, 0, 0]
    
    ''' initializing topology '''
    TSNdf = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
    TSNdf.plot_topology()
    
    # Define initial conditions
    #                  dA1  dA2  dA3
    dA_tot = np.array([0, 250, 0]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6
    G_tot = np.array([50, 50, 50, 50, 50, 50]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4 G5 G6
    G_int_vec = [0, 0, 1, 1, 0, 0]
    
    ''' initializing initial conditions '''
    TSNdf.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions
    TSNdf.plot_topology()
    TSNdf.plot_topology(show_rnas=0)
    
    t_vec1 = np.linspace(0,0.5,101)*3600 # seconds
    
    ''' simulating the TSN model for initial RNA production to set State 2'''
    TSNdf.simulate(t_vec1,1)
    
    t_vec2 = np.linspace(t_vec1[-1]/3600,1,1001)*3600 # seconds
    
    ''' simulating the TSN model after added the other activators to start experiment '''
    TSNdf.simulate(t_vec2,2,dA=[250,'NA',250])
    
    t_vec3 = np.linspace(t_vec2[-1]/3600,3.5,1001)*3600 # seconds
    
    ''' simulating the TSN model after adding rIr1 to change to State 1 '''
    TSNdf.simulate(t_vec3,3,rIr=[10000,'NA','NA'])
    
    t_vec4 = np.linspace(t_vec3[-1]/3600,8.5,1001)*3600 # seconds
    
    ''' simulating the TSN model after adding rIr3 to change to State 3 '''
    TSNdf.simulate(t_vec4,4,rIr=['NA','NA',25000])
    
    ''' Plotting solution '''
    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = TSNdf.output_concentration['GdA1']
    G3 = TSNdf.output_concentration['GdA3']
    G5 = TSNdf.output_concentration['GdA5']
    
    R1 = TSNdf.output_concentration['rR1']
    R2 = TSNdf.output_concentration['rR2']
    R3 = TSNdf.output_concentration['rR3']
    
    rI1 = TSNdf.output_concentration['rIr1']
    rI2 = TSNdf.output_concentration['rIr2']
    rI3 = TSNdf.output_concentration['rIr3']
    
    sim_t = TSNdf.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('TSN: State 2 to State 1 to State 3',fontsize=fs+2,weight='bold')
    plt.subplot(3,3,1)
    # 0.5 hr shift as that is when the experiment actually starts here
    plt.plot(sim_t/3600-0.5,G1/G_tot[0],color='b',linewidth=4)
    plt.plot(sim_t/3600-0.5,G3/G_tot[2],color='r',linewidth=4)
    plt.plot(sim_t/3600-0.5,G5/G_tot[4],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,8)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    fs = 13
    plt.subplot(3,3,2)
    plt.plot(sim_t/3600-0.5,R1,color='b',linewidth=3,linestyle=':')
    plt.plot(sim_t/3600-0.5,R2,color='r',linewidth=3,linestyle=':')
    plt.plot(sim_t/3600-0.5,R3,color=[0.4,0,0.4],linewidth=3,linestyle=':')
    plt.ylabel('[ Repressor ] nM',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,8)
    #ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    fs = 13
    plt.subplot(3,3,3)
    plt.plot(sim_t/3600-0.5,rI1,color='b',linewidth=3,linestyle=':')
    plt.plot(sim_t/3600-0.5,rI2,color='r',linewidth=3,linestyle=':')
    plt.plot(sim_t/3600-0.5,rI3,color=[0.4,0,0.4],linewidth=3,linestyle=':')
    plt.ylabel('[ Inducer ] nM',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,8)
    #ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 6
###############################################################################
This simulates a IFFL with 150 nM excess blocker for G3 and G4
-------------------------------------------------------------------------------
   -This example shows how to import your own rate constant values
###############################################################################
'''
if ex6 == 1:
    #          G1 G2 G3 G4      
    act_vec =  [1, 1, 2, 3]
    blk_vec =  [0, 0, 2, 3]
    prod_vec = [3, 2,-3, 0]
    indc_vec = [0, 0, 0, 0]
    
    # below this demonstrates the different ways you can input reaction rate constants
    # it can be a list of length ind_nodes for production rates
    # it can be a list of ortho_nodes for other rates
    # Or it can be a single value which will be assumed to be the same for all nodes
    # defining the rate constants
         # G1   G2   G3  G4
    kpr = [0.1,0.2,0.15,0.05] # repressor production rates
    kpc = [0.1,0.2,0.15,0.05]# coactivator production rates
    kpi = [0.1,0.2,0.15,0.05] # inducer production rates
          # rC2   rC3    rR3
    kd_H = [0.003,0.002,0.0025] # RNaseH degradation rates
    kd_A = [0.003,0.001,0.004] # RNaseA degradation rates
   #    All nodes
    kga = [1e4] # activation rates
    kgar = 5e3 # repression rates
    kar = 1e4 # activator inhibition rates
    kgb = 1e4 # free blocking rates
    kgbc = 5e3# coactivation rates
    kbc = 1e4# blocker inhibition rates
    kgab = 5e3 # active blocking rates
    kir = 1e4 # inducer binding rates
    
    my_rc = [kpr,kpc,kpi,kd_H,kd_A,kga,kgar,kar,kgb,kgbc,kbc,kgab,kir]
    
    ''' initializing TSN topology '''
    IFFL = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
    IFFL.plot_topology()
    IFFL.plot_topology(show_rnas=0)
    
    # Define initial conditions
    #                  dA1  dA2  dA3
    dA_tot = np.array([250, 250, 250]) # total activator added
    #                 G1  G2  G3  G4
    G_tot = np.array([25, 5, 25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4
    G_int_vec = [0, 0,-1,-1]
    
    ''' initializing initial conditions '''
    IFFL.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added = [0,150,150]) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds
    
    ''' simulating the IFFL'''
    IFFL.simulate(t_vec1,1,rate_constants=my_rc)
    
    # pulling out the desired concentrations for plotting
    G1 = IFFL.output_concentration['GdA1']
    G3 = IFFL.output_concentration['GdA3']
    G4 = IFFL.output_concentration['GdA4']
    
    sim_t = IFFL.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('IFFL with added blockers',fontsize=fs+2,weight='bold')
    plt.subplot(3,4,1)
    plt.plot(sim_t/3600,G1/G_tot[0],color='g',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,5)
    plt.plot(sim_t/3600,G3/G_tot[2],color='r',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,9)
    plt.plot(sim_t/3600,G4/G_tot[3],color='b',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

'''
###############################################################################
###############################################################################
SPATIOTEMPORAL GENELET EXAMPLES
###############################################################################
###############################################################################
'''

''' 
EXAMPLE 7
###############################################################################
This simulates setting the TSN into State 3 by initally including rR1 and rR2
and then switching to State 1 by addding rIr1 after 2.5 hours 
-------------------------------------------------------------------------------
   -This example illustrates how to use the STG model with RNase A
   -This example illustrates the use of multiple iterations of model.simulate
   -This example also illustrates how to add non-default initial conditions for rR
   -This example shows how to export the symbolic expression to a text file
###############################################################################
'''
if ex7 == 1:
    #          G1 G2 G3 G4 G5 G6      
    act_vec =  [1, 1, 2, 2, 3, 3]
    prod_vec = [-3,-2,-1,-3,-1,-2]
    indc_vec = [0, 0, 0, 0, 0, 0]
    
    ''' initializing TSN topology '''
    # blk_vec is not needed as an input here
    TSNsg = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,genelet_type='STG')
    TSNsg.plot_topology()
    TSNsg.plot_topology(show_rnas=0)
    
    # Define initial conditions
    #                  dA1  dA2  dA3
    dA_tot = np.array([500, 500, 500]) # total activator added
    #                 G1  G2  G3  G4  G5  G6
    G_tot = np.array([25, 25, 25, 25, 25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON) OFF here is GrR for STG genelets
    #            G1 G2 G3 G4 G5 G6
    G_int_vec = [0, 0, 0, 0, 1, 1]
    
    ''' initializing initial conditions '''
    # specifying an initial amount of rR1 and rR2 to set the system into State 3
    TSNsg.initial_conditions(dA_tot,G_tot,G_int_vec,rRin=[1000,1000,0]) # default of 0 for all other initial conditions
    
    t_vec1 = np.linspace(0,2.5,1001)*3600 # seconds
    
    ''' simulating the TSN model for initial RNA production to set State 2'''
    TSNsg.simulate(t_vec1,1,rnase='RnA')
    
    t_vec2 = np.linspace(t_vec1[-1]/3600,6,1001)*3600 # seconds
    
    ''' simulating the TSN model after adding rIr1 to change to State 1 '''
    TSNsg.simulate(t_vec2,2,rnase='RnA',rIr=[15000,'NA','NA'])
    
    
    ''' Plotting solution '''
    # pulling out solutions for plotting from output_concentration attribute
    G1 = TSNsg.output_concentration['Gon1']
    G3 = TSNsg.output_concentration['Gon3']
    G6 = TSNsg.output_concentration['Gon6']
    
    R1 = TSNsg.output_concentration['rR1']
    R2 = TSNsg.output_concentration['rR2']
    R3 = TSNsg.output_concentration['rR3']
    
    A1 = TSNsg.output_concentration['dA1']
    A2 = TSNsg.output_concentration['dA2']
    A3 = TSNsg.output_concentration['dA3']
    
    sim_t = TSNsg.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('TSN: State 3 to State 1',fontsize=fs+2,weight='bold')
    plt.subplot(3,3,1)
    plt.plot(sim_t/3600,G1/G_tot[0],color='b',linewidth=4)
    plt.plot(sim_t/3600,G3/G_tot[2],color='r',linewidth=4)
    plt.plot(sim_t/3600,G6/G_tot[5],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,2)
    plt.plot(sim_t/3600,R1,color='b',linewidth=4)
    plt.plot(sim_t/3600,R2,color='r',linewidth=4)
    plt.plot(sim_t/3600,R3,color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('[Repressors]',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    #ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,3)
    plt.plot(sim_t/3600,A1,color='b',linewidth=4)
    plt.plot(sim_t/3600,A2,color='r',linewidth=4)
    plt.plot(sim_t/3600,A3,color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('[Activators]',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    #ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    ''' Exporting symbolic equations '''
    TSNsg.export_sym_eqs(file_name='stg_testing_sym_export', \
                         path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\')


''' 
EXAMPLE 8
###############################################################################
This simulates a G1-|G2-|G3 (aka G1->G3) STG model with 5% leak
This simulates the network both with and without 150 nM dA1,dA2,dA3 which suppress the leak
-------------------------------------------------------------------------------
   -This example illustrates how to include the leak reactions and use both RNase H and A
   -This example illustrates how to export the STG symbolic equations
###############################################################################
'''
if ex8 == 1:
    #          G1 G2 G3   
    act_vec =  [1, 2, 3]
    prod_vec = [-2,-3,0]
    indc_vec = [0, 0, 0]
    
    
    ''' initializing topology '''
    G123 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,genelet_type='STG')
    G123.plot_topology()
    G123.plot_topology(show_rnas=0)
    
    # Define initial conditions
    G_tot = [0]*2
    #                 dA1 dA2 dA3
    dA_tot = np.array([0,  0,  0]) # total activator added
    #                    G1  G2  G3
    G_tot[0] = np.array([25,  25, 25]) # total genelet added
    G_tot[1] = np.array([ 0,  25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3
    G_int_vec = [1, 1, 0]
    
    G1 = [0]*2
    G2 = [0]*2
    G3 = [0]*2
    ls = ['-',':']
    '''
    0 nM ADDED activators
    '''
    plt.figure()
    for i in range(len(G_tot)):
 
        ''' initializing initial conditions '''
        G123.initial_conditions(dA_tot,G_tot[i],G_int_vec) # default of 0 for all other initial conditions
        
        t_vec1 = np.linspace(0,3,1001)*3600 # seconds
        
        ''' simulating the network '''
        G123.simulate(t_vec1,1,rnase='both',leak=0.05)
        
        # pulling out the desired concentrations for plotting
        
        G1[i] = G123.output_concentration['Gon1']
        G2[i] = G123.output_concentration['Gon2']
        G3[i] = G123.output_concentration['Gon3']
        
        sim_t = G123.sol.t
        
        fs = 13
        
        plt.suptitle('STG G1-|G2-|G3 with leak',fontsize=fs+2,weight='bold')
        plt.subplot(3,4,2)
        plt.title('No addded activators')
        plt.plot(sim_t/3600,G1[i]/G_tot[0][0],color='g',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,6)
        plt.plot(sim_t/3600,G2[i]/G_tot[0][1],color='r',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,10)
        plt.plot(sim_t/3600,G3[i]/G_tot[0][2],color='b',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xlabel('time (hr)',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
    G1 = [0]*2
    G2 = [0]*2
    G3 = [0]*2
     #                 dA1  dA2  dA3
    dA_tot = np.array([ 0,  0,  500]) # total activator added
    
    '''
    500 nM ADDED dA3
    '''
    for i in range(len(G_tot)):
        ''' initializing initial conditions '''
        G123.initial_conditions(dA_tot,G_tot[i],G_int_vec) # default of 0 for all other initial conditions
        
        t_vec1 = np.linspace(0,3,1001)*3600 # seconds
        
        ''' simulating the network'''
        G123.simulate(t_vec1,1,rnase='both',leak=0.05)
        
        # pulling out concentrations for plotting 
        G1[i] = G123.output_concentration['Gon1']
        G2[i] = G123.output_concentration['Gon2']
        G3[i] = G123.output_concentration['Gon3']
        
        sim_t = G123.sol.t
        
        fs = 13
        
        plt.subplot(3,4,4)
        plt.title('500 nM dA3')
        plt.plot(sim_t/3600,G1[i]/G_tot[0][0],color='g',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,8)
        plt.plot(sim_t/3600,G2[i]/G_tot[0][1],color='r',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
        
        plt.subplot(3,4,12)
        plt.plot(sim_t/3600,G3[i]/G_tot[0][2],color='b',linewidth=4,linestyle=ls[i])
        plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
        plt.xlabel('time (hr)',fontsize=fs,weight='bold')
        plt.xticks(fontsize=fs,weight='bold')
        plt.yticks(fontsize=fs,weight='bold')
        ax1 = plt.gca()
        ax1.set_xlim(0,3)
        ax1.set_ylim(-0.1,1.1)
        ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
        ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

    ''' Exporting symbolic equations '''
    G123.export_sym_eqs(file_name='stg_leak_testing_sym_export', \
                         path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\')

''' 
EXAMPLE 9
###############################################################################
This simulates a G1-|G3 and G2-|rR3 (substraction circuit) STG model with 5% leak
This simulates the network both with [G1] > [G2] and [G2] > [G1]
-------------------------------------------------------------------------------
   -This example illustrates how to include the leak reactions and use both RNase H and A
   -This example illustrates how to export the STG symbolic equations
###############################################################################
'''
if ex9 == 1:
    
#    #          G1 G2 G3   
#    act_vec =  [1, 2, 3]
#    prod_vec = [-3,0 ,0]
#    indc_vec = [0,-3, 0]
    
    topology_mat = [[0,0,0],[0,0,0],[-1,0,0]]
    I_topology_mat = [[0,0,0],[0,0,0],[0,-1,0]]
    
    ''' initializing topology '''
    #G123 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,genelet_type='STG')
    G123 = GGM.GeneletNetwork(top_mat=topology_mat,Itop_mat=I_topology_mat,genelet_type='STG')
    G123.plot_topology()
    G123.plot_topology(show_rnas=0)
    
    # Define initial conditions
    G_tot = [0]*2
    #                 dA1 dA2 dA3
    dA_tot = np.array([0,  0,  0]) # total activator added
    #                    G1  G2  G3
    G_tot[0] = np.array([25,  5, 25]) # total genelet added
    G_tot[1] = np.array([ 5,  25, 25]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3
    G_int_vec = [1, 1, 1]
    
    G1 = [0]*2
    G2 = [0]*2
    G3 = [0]*2
    ls = ['-',':']
    '''
    0 nM ADDED activators
    '''
    plt.figure()
    for i in range(len(G_tot)):
 
        ''' initializing initial conditions '''
        G123.initial_conditions(dA_tot,G_tot[i],G_int_vec) # default of 0 for all other initial conditions
        
        t_vec1 = np.linspace(0,3,1001)*3600 # seconds
        
        ''' simulating the network '''
        G123.simulate(t_vec1,1,rnase='RnA',leak=0.05)
        
        # pulling out the desired concentrations for plotting
        
        G1[i] = G123.output_concentration['Gon1']
        G2[i] = G123.output_concentration['Gon2']
        G3[i] = G123.output_concentration['Gon3']
        
    sim_t = G123.sol.t
    
    fs = 13
    
    plt.suptitle('STG subtraction circuit with leak',fontsize=fs+2,weight='bold')
    plt.subplot(3,4,2)
    plt.title('[TN1] > [TN2]')
    plt.plot(sim_t/3600,G1[0]/G_tot[0][0],color='g',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,6)
    plt.plot(sim_t/3600,G2[0]/G_tot[0][1],color='r',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,10)
    plt.plot(sim_t/3600,G3[0]/G_tot[0][2],color='b',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,4)
    plt.title('[TN2] > [TN1]')
    plt.plot(sim_t/3600,G1[1]/G_tot[1][0],color='g',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,8)
    plt.plot(sim_t/3600,G2[1]/G_tot[1][1],color='r',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,4,12)
    plt.plot(sim_t/3600,G3[1]/G_tot[1][2],color='b',linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,3)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 10
###############################################################################
This simulates a double flip from State 2 to State 1 to State 3
-------------------------------------------------------------------------------
   -Complex example of a TSN with three inducer producing nodes for changing states
   -This shows how to change numerous concentrations per iterative simulation
###############################################################################
'''
if ex10 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9     
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6]
    prod_vec = [-3,-2,-1,-3,-1,-2,0, 0, 0]
    indc_vec = [0, 0, 0, 0, 0, 0,-1,-2,-3]
    
    ''' initializing topology '''
    TSNinodes = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6
    dA_tot = np.array([250, 250, 250,  0,   0,   0]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9
    G_tot = np.array([25, 25, 25, 25, 25, 25, 50, 50, 50]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #            G1 G2 G3 G4 G5 G6 G7 G8 G9
    G_int_vec = [1, 1, 0, 0, 0, 0, 0, 0, 0]
    
    ''' initializing initial conditions '''
    TSNinodes.initial_conditions(dA_tot,G_tot,G_int_vec,rRin=[0,1000,1000,0,0,0]) # default of 0 for all other initial conditions
    TSNinodes.plot_topology(layout='shell')
    TSNinodes.plot_topology(layout='shell',show_rnas=0)
    
    t_vec1 = np.linspace(0,1,101)*3600 # seconds
    
    ''' simulating the TSN model for initial hour in State 1 '''
    TSNinodes.simulate(t_vec1,1)
    
    t_vec2 = np.linspace(t_vec1[-1]/3600,4,1001)*3600 # seconds
    
    ''' simulating the TSN model after adding the activator for G8 to induce a switch to State 2 '''
    TSNinodes.simulate(t_vec2,2,dA=['NA','NA','NA','NA',250,'NA'])
    
    t_vec3 = np.linspace(t_vec2[-1]/3600,15,1001)*3600 # seconds
    
    ''' simulating the TSN model after adding the activator for G9 and shutting off G8 to induce a switch to State 3 '''
    TSNinodes.simulate(t_vec3,3,dA=['NA','NA','NA','NA','NA',250],dR=['NA','NA','NA','NA',500,'NA'])
    
    ''' Plotting solution '''
    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = TSNinodes.output_concentration['GdA1']
    G3 = TSNinodes.output_concentration['GdA3']
    G5 = TSNinodes.output_concentration['GdA5']
    G7 = TSNinodes.output_concentration['GdA7']
    G8 = TSNinodes.output_concentration['GdA8']
    G9 = TSNinodes.output_concentration['GdA9']
    
    sim_t = TSNinodes.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('TSN: State 1 to State 2 to State 3',fontsize=fs+2,weight='bold')
    plt.subplot(3,3,4)
    plt.plot(sim_t/3600,G1/G_tot[0],color='b',linewidth=4)
    plt.plot(sim_t/3600,G3/G_tot[2],color='r',linewidth=4)
    plt.plot(sim_t/3600,G5/G_tot[4],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,15)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,1)
    plt.plot(sim_t/3600,G7/G_tot[6],color='b',linewidth=2,linestyle=':')
    plt.plot(sim_t/3600,G8/G_tot[7],color='r',linewidth=2,linestyle=':')
    plt.plot(sim_t/3600,G9/G_tot[8],color=[0.4,0,0.4],linewidth=2,linestyle=':')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,15)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    ''' Exporting symbolic equations '''
    TSNinodes.export_sym_eqs(file_name='TSN_inodes_eqs', \
                         path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\')

''' 
EXAMPLE 11
###############################################################################
BS_IFFL1|2 network with excess blocker in IFFL modules
-------------------------------------------------------------------------------
   -Complex example of a large network
###############################################################################
'''
if ex11 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12   
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6, 6,  7,  8]
    prod_vec = [-2,3,-1, 6, 4, 5,-5, 0, 7, 8, -8,  0]
    
    ''' initializing topology '''
    BS_IFFL12 = GGM.GeneletNetwork(act_vec,prod_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6  dA7  dA8
    dA_tot = np.array([250, 250, 250, 250, 250, 250, 250, 250]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10 G11 G12
    G_tot = np.array([25, 15, 25, 15,  5, 25, 25, 25,  5, 25, 25, 25]) # total genelet added
    
    G_int_vec = [0]*2
    rR_int = [0]*2
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #               G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12
    G_int_vec[0] = [1, 1, 0, 0,-1,-1, -1,-1,-1,-1, -1, -1]
    G_int_vec[1] = [0, 0, 1, 1,-1,-1, -1,-1,-1,-1, -1, -1]
    rR_int[0] = [0,1000,0,0,0,0,0,0]
    rR_int[1] = [1000,0,0,0,0,0,0,0]
    
    G1=[0]*2
    G3=[0]*2
    G8=[0]*2
    G12=[0]*2
    
    for i in range(len(rR_int)):
        ''' initializing initial conditions '''
        BS_IFFL12.initial_conditions(dA_tot,G_tot,G_int_vec[i],dB_added=[0,0,150,150,150,150,150,150],rRin=rR_int[i]) # default of 0 for all other initial conditions
        BS_IFFL12.plot_topology(show_rnas=1,layout='shell',plot_title='State '+str(i+1))
        BS_IFFL12.plot_topology(show_rnas=0,layout='shell',plot_title='State '+str(i+1))
        
        t_vec1 = np.linspace(0,5,101)*3600 # seconds
        
        ''' simulating the TSN model for initial hour in State 1 '''
        BS_IFFL12.simulate(t_vec1,1)
    
        ''' Plotting solution '''
        
        # pulling out solutions for plotting from output_concentration attribute
        G1[i] = BS_IFFL12.output_concentration['GdA1']
        G3[i] = BS_IFFL12.output_concentration['GdA3']
        G8[i] = BS_IFFL12.output_concentration['GdA8']
        G12[i] = BS_IFFL12.output_concentration['GdA12']
    
    sim_t = BS_IFFL12.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
    plt.subplot(3,3,4)

    plt.plot(sim_t/3600,G8[0]/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12[0]/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    ax1 = plt.gca()
    ax1.set_xlim(0,5)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,1)    
    plt.title('State 1')
    plt.plot(sim_t/3600,G1[0]/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3[0]/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0)
    ax1 = plt.gca()
    ax1.set_xlim(0,5)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.suptitle('BS_IFFL1|2: State 2',fontsize=fs+2,weight='bold')
    plt.subplot(3,3,6)
   
    plt.plot(sim_t/3600,G8[1]/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12[1]/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    ax1 = plt.gca()
    ax1.set_xlim(0,5)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,3)
    plt.title('State 2')
    plt.plot(sim_t/3600,G1[1]/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3[1]/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0)
    ax1 = plt.gca()
    ax1.set_xlim(0,5)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')


''' 
EXAMPLE 12
###############################################################################
I_BS_IFFL1|2 network with excess blocker in IFFL modules
-------------------------------------------------------------------------------
   -Complex example of a large network
###############################################################################
'''
if ex12 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14   
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6, 6,  7,  8,  9,  10]
    prod_vec = [-2,3,-1, 6, 4, 5,-5, 0, 7, 8, -8,  0,  0,  0]
    indc_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, -1, -2]
    
    ''' initializing topology '''
    I_BS_IFFL12 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6  dA7  dA8  dA9  dA10
    dA_tot = np.array([250, 250, 250, 250, 250, 250, 250, 250,  0,  0]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10 G11 G12 G13 G15
    G_tot = np.array([25, 15, 25, 15,  5, 25, 25, 25,  5, 25, 25, 25, 100, 100]) # total genelet added
    
    G_int_vec = [0]*2
    rR_int = [0]*2
    dA_add = [0]*2
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #               G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14
    G_int_vec[0] = [1, 1, 0, 0,-1,-1, -1,-1,-1,-1, -1, -1,  0,  0]
    G_int_vec[1] = [0, 0, 1, 1,-1,-1, -1,-1,-1,-1, -1, -1,  0,  0]
    rR_int[0] = [0,1000,0,0,0,0,0,0,0,0]
    rR_int[1] = [1000,0,0,0,0,0,0,0,0,0]
    dA_add[0] = ['NA','NA','NA','NA','NA','NA','NA','NA','NA',250]
    dA_add[1] =  ['NA','NA','NA','NA','NA','NA','NA','NA',250,'NA']
    
    G1=[0]*2
    G3=[0]*2
    G8=[0]*2
    G12=[0]*2
    G13=[0]*2
    G14=[0]*2
    GB13=[0]*2
    GB14=[0]*2
    A9=[0]*2
    A10=[0]*2
    
    for i in range(len(rR_int)):
        ''' initializing initial conditions '''
        I_BS_IFFL12.initial_conditions(dA_tot,G_tot,G_int_vec[i],dB_added=[0,0,150,150,150,150,150,150,0,0],rRin=rR_int[i]) # default of 0 for all other initial conditions
        I_BS_IFFL12.plot_topology(show_rnas=1,layout='shell',plot_title='State '+str(i+1))
        #I_BS_IFFL12.plot_topology(show_rnas=0,layout='shell',plot_title='State '+str(i+1))
        
        t_vec1 = np.linspace(0,1.5,101)*3600 # seconds
        
        ''' simulating the model before turning on induce nodes  '''
        I_BS_IFFL12.simulate(t_vec1,1)
    
        t_vec2 = np.linspace(t_vec1[-1]/3600,6,1001)*3600 # seconds
        
        ''' simulating the model  '''
        I_BS_IFFL12.simulate(t_vec2,2,dA=dA_add[i])    
        
        # pulling out solutions for plotting from output_concentration attribute
        G1[i] = I_BS_IFFL12.output_concentration['GdA1']
        G3[i] = I_BS_IFFL12.output_concentration['GdA3']
        G8[i] = I_BS_IFFL12.output_concentration['GdA8']
        G12[i] = I_BS_IFFL12.output_concentration['GdA12']
        G13[i] = I_BS_IFFL12.output_concentration['GdA13']
        G14[i] = I_BS_IFFL12.output_concentration['GdA14']
        GB13[i] = I_BS_IFFL12.output_concentration['GdB13']
        GB14[i] = I_BS_IFFL12.output_concentration['GdB14']
        A9[i] = I_BS_IFFL12.output_concentration['dA9']
        A10[i] = I_BS_IFFL12.output_concentration['dA10']
    
    ''' Exporting symbolic equations '''
    I_BS_IFFL12.export_sym_eqs(file_name='I_BS_IFFL12_eqs', \
                         path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\')

    
    ''' Plotting solution '''
    sim_t = I_BS_IFFL12.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
   
    plt.subplot(3,3,1)   
    plt.title('State 1 to State 2')
    plt.plot(sim_t/3600,G13[0]/G_tot[12],color=[0.5,0.5,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G14[0]/G_tot[13],color=[1,0.75,0],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN9','TN10'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,4)    
    plt.plot(sim_t/3600,G1[0]/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3[0]/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,7)
    plt.plot(sim_t/3600,G8[0]/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12[0]/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0,loc=9)
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,3)
    plt.title('State 2 to State 1')
    plt.plot(sim_t/3600,G13[1]/G_tot[12],color=[0.5,0.5,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G14[1]/G_tot[13],color=[1,0.75,0],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN9','TN10'],frameon=0)
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,6)
    plt.plot(sim_t/3600,G1[1]/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3[1]/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

    plt.subplot(3,3,9)
    plt.plot(sim_t/3600,G8[1]/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12[1]/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0,loc=9)
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
''' 
EXAMPLE 13
###############################################################################
FBI_BS_IFFL1|2 network with excess blocker in IFFL modules
-------------------------------------------------------------------------------
   -Complex example of a large network
###############################################################################
'''
if ex13 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14   
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6, 6,  7,  8,  9,  10]
    prod_vec = [-2,3,-1, 6, 4, 5,-5, 0, 7, 8, -8,  9,  0,  0]
    indc_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, -1, -2]
    
    ''' initializing topology '''
    FBI_BS_IFFL12 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6  dA7  dA8  dA9  dA10
    dA_tot = np.array([250, 250, 250, 250, 250, 250, 250, 250,  250,  0]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10 G11 G12 G13  G14
    G_tot = np.array([25, 15, 25, 15,  5, 25, 25, 15,  5, 25, 25, 15, 75, 75]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #           G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14
    G_int_vec = [0, 0, 1, 1, -1,-1, -1,-1,-1,-1, -1, -1, -1, 0]
    rR_int = [1000,0,0,0,0,0,0,0,0,0]
    #dA_add = ['NA','NA','NA','NA','NA','NA','NA','NA','NA',250]
 
    ''' initializing initial conditions '''
    FBI_BS_IFFL12.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,0,150,150,150,150,150,150,250,0],rRin=rR_int) # default of 0 for all other initial conditions
    FBI_BS_IFFL12.plot_topology(show_rnas=1,layout='shell',plot_title='State 1')
    #FBI_BS_IFFL12.plot_topology(show_rnas=0,layout='shell',plot_title='State 1')
    
    t_vec1 = np.linspace(0,6,1001)*3600 # seconds
    
    ''' simulating the model before turning on induce nodes  '''
    FBI_BS_IFFL12.simulate(t_vec1,1)

#    t_vec2 = np.linspace(t_vec1[-1]/3600,6,1001)*3600 # seconds
#    
#    ''' simulating the model  '''
#    I_BS_IFFL12.simulate(t_vec2,2,dA=dA_add[i])    
#    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = FBI_BS_IFFL12.output_concentration['GdA1']
    G3 = FBI_BS_IFFL12.output_concentration['GdA3']
    G8 = FBI_BS_IFFL12.output_concentration['GdA8']
    G12 = FBI_BS_IFFL12.output_concentration['GdA12']
    G13 = FBI_BS_IFFL12.output_concentration['GdA13']
    G14 = FBI_BS_IFFL12.output_concentration['GdA14']

    ''' Plotting solution '''
    sim_t = FBI_BS_IFFL12.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
       
    plt.subplot(3,3,1)   
    plt.title('State 2 -> State 1')
    plt.plot(sim_t/3600,G13/G_tot[12],color=[0.5,0.5,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G14/G_tot[13],color=[1,0.75,0],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN9','TN10'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,4)    
    plt.plot(sim_t/3600,G1/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,7)
    plt.plot(sim_t/3600,G8/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0,loc=9)
    ax1 = plt.gca()
    ax1.set_xlim(0,6)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 14
###############################################################################
FBI_BS_IFFL1|2 network with excess blocker in IFFL modules
-------------------------------------------------------------------------------
   -Here we start in State 1 which is stable and we turn on TN10 to switch to
    State 2 which executes its temporal program and then switchs back to State 1
###############################################################################
'''
if ex14 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14   
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6, 6,  7,  8,  9,  10]
    prod_vec = [-2,3,-1, 6, 4, 5,-5, 0, 7, 8, -8,  9,  0,  0]
    indc_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, -1, -2]
    
    ''' initializing topology '''
    FBI_BS_IFFL12 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6  dA7  dA8  dA9  dA10
    dA_tot = np.array([250, 250, 250, 250, 325, 250, 250, 320,  250,  0]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10 G11 G12 G13  G14
    G_tot = np.array([25, 10, 25, 10,  5, 10, 10, 15,  5, 10, 10, 15, 75, 75]) # total genelet added
    
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #           G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14
    G_int_vec = [1, 1, 0, 0, -1,-1, -1,-1,-1,-1, -1, -1,  -1,  0]
    rR_int = [0,1000,0,0,0,0,0,0,0,0]
    dA_add = ['NA','NA','NA','NA','NA','NA','NA','NA','NA',250]
 
    ''' initializing initial conditions '''
    FBI_BS_IFFL12.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,0,320,320,320,320,320,320,250,0],rRin=rR_int) # default of 0 for all other initial conditions
    FBI_BS_IFFL12.plot_topology(show_rnas=1,layout='shell',plot_title='State 1')
    #I_BS_IFFL12.plot_topology(show_rnas=0,layout='shell',plot_title='State '+str(i+1))
    
    t_vec1 = np.linspace(0,1.5,1001)*3600 # seconds
    
    ''' simulating the model before turning on induce nodes  '''
    FBI_BS_IFFL12.simulate(t_vec1,1)

    t_vec2 = np.linspace(t_vec1[-1]/3600,3,1001)*3600 # seconds
    
    ''' simulating the model  '''
    FBI_BS_IFFL12.simulate(t_vec2,2,dA=dA_add)   

    t_vec3 = np.linspace(t_vec2[-1]/3600,15,1001)*3600 # seconds
    
    ''' simulating the model  '''
    FBI_BS_IFFL12.simulate(t_vec3,3,dR=['NA','NA','NA','NA','NA','NA','NA','NA','NA',500])    
    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = FBI_BS_IFFL12.output_concentration['GdA1']
    G3 = FBI_BS_IFFL12.output_concentration['GdA3']
    G8 = FBI_BS_IFFL12.output_concentration['GdA8']
    G12 = FBI_BS_IFFL12.output_concentration['GdA12']
    G13 = FBI_BS_IFFL12.output_concentration['GdA13']
    G14 = FBI_BS_IFFL12.output_concentration['GdA14']

    ''' Plotting solution '''
    sim_t = FBI_BS_IFFL12.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
       
    plt.subplot(3,3,1)   
    plt.title('State 1 induce State 2 -> State 1')
    plt.plot(sim_t/3600,G13/G_tot[12],color=[0.5,0.5,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G14/G_tot[13],color=[1,0.75,0],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN9','TN10'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,15)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,4)    
    plt.plot(sim_t/3600,G1/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,15)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,7)
    plt.plot(sim_t/3600,G8/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0,loc=9)
    ax1 = plt.gca()
    ax1.set_xlim(0,15)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    ''' Compiling genelet sequences to build the network '''
    input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
    file_name='FBI_BS_IFFL12_sequences'
    path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\'

    FBI_BS_IFFL12.compile_network_sequences(input_file_loc,bth=0,save_file_name=file_name,save_path_name=path_name)

''' 
EXAMPLE 15
###############################################################################
FBI_BS_IFFL1|2 network with excess blocker in IFFL modules fully connected
-------------------------------------------------------------------------------
   -Complex example of a large network
###############################################################################
'''
if ex15 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14   
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 5, 6, 6,  7,  8,  9,  10]
    prod_vec = [-2,3,-1, 6, 4, 5,-5, 10, 7, 8, -8,  9,  0,  0]
    indc_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0, -1, -2]
    
    ''' initializing topology '''
    FBI_BS_IFFL12 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6  dA7  dA8  dA9  dA10
    dA_tot = np.array([250, 250, 250, 250, 435, 250, 250, 435, 250,  250]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10 G11 G12 G13  G14
    G_tot = np.array([25, 10, 25, 10,  5, 10, 10, 9,  5, 10, 10, 9, 42, 42]) # total genelet added
     
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #           G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14
    G_int_vec = [1, 1, 0, 0,-1,-1, -1,-1,-1,-1, -1, -1,  -1,  -1]
    rR_int = [0,1000,0,0,0,0,0,0,0,0]
    #dA_add = ['NA','NA','NA','NA','NA','NA','NA','NA','NA',250]
 
    ''' initializing initial conditions '''
    FBI_BS_IFFL12.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,0,355,355,355,355,355,355,357.75,360],rRin=rR_int) # default of 0 for all other initial conditions
    FBI_BS_IFFL12.plot_topology(show_rnas=1,layout='circular',plot_title='State 1')
    #I_BS_IFFL12.plot_topology(show_rnas=0,layout='shell',plot_title='State '+str(i+1))
    
    t_vec1 = np.linspace(0,45,501)*3600 # seconds
    
    ''' simulating the model before turning on induce nodes  '''
    FBI_BS_IFFL12.simulate(t_vec1,1)

#    t_vec2 = np.linspace(t_vec1[-1]/3600,6,1001)*3600 # seconds
#    
#    ''' simulating the model  '''
#    I_BS_IFFL12.simulate(t_vec2,2,dA=dA_add[i])    
#    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = FBI_BS_IFFL12.output_concentration['GdA1']
    G3 = FBI_BS_IFFL12.output_concentration['GdA3']
    G8 = FBI_BS_IFFL12.output_concentration['GdA8']
    G12 = FBI_BS_IFFL12.output_concentration['GdA12']
    G13 = FBI_BS_IFFL12.output_concentration['GdA13']
    G14 = FBI_BS_IFFL12.output_concentration['GdA14']

    ''' Plotting solution '''
    sim_t = FBI_BS_IFFL12.sol.t
    
    fs = 13
    plt.figure()
    plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
       
    plt.subplot(3,3,1)   
    plt.title('State 1 to State 2')
    plt.plot(sim_t/3600,G13/G_tot[12],color=[0.5,0.5,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G14/G_tot[13],color=[1,0.75,0],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN9','TN10'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,45)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,4)    
    plt.plot(sim_t/3600,G1/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,45)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,7)
    plt.plot(sim_t/3600,G8/G_tot[7],color='b',linewidth=4)
    plt.plot(sim_t/3600,G12/G_tot[11],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0,loc=9)
    ax1 = plt.gca()
    ax1.set_xlim(0,45)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')

''' 
EXAMPLE 16
###############################################################################
Simpler induction / bistable / pulse feedback oscillator
-------------------------------------------------------------------------------
   -Complex example of a large network
###############################################################################
'''
if ex16 == 1:
    #          G1 G2 G3 G4 G5 G6 G7 G8 G9 G10   
    act_vec =  [1, 1, 2, 2, 3, 3, 4, 4, 5, 6]
    prod_vec = [-2,3,-1, 4, 6,-3, 5,-4, 0, 0]
    indc_vec = [0, 0, 0, 0, 0, 0, 0, 0,-1,-2]
    
    ''' initializing topology '''
    S_BS_OC = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec)
    
    # Define initial conditions
    #                  dA1  dA2  dA3  dA4  dA5  dA6
    dA_tot = np.array([250, 250, 250, 250, 250, 250]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4  G5  G6  G7  G8  G9  G10
    G_tot = np.array([25, 15, 25, 15, 15, 10, 15, 10, 42, 42]) # total genelet added
     
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #           G1 G2 G3 G4 G5 G6 G7 G8 G9 G10
    G_int_vec = [1, 1, 0, 0,-1,-1, -1,-1,-1,-1]
    rR_int = [0,1000,0,0,0,0]
    #dA_add = ['NA','NA','NA','NA','NA','NA','NA','NA','NA',250]
 
    ''' initializing initial conditions '''
    S_BS_OC.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,0,150,150,150,150],rRin=rR_int) # default of 0 for all other initial conditions
    S_BS_OC.plot_topology(show_rnas=1,layout='circular')
    S_BS_OC.plot_topology(show_rnas=0,layout='shell')
    
    t_vec1 = np.linspace(0,10,501)*3600 # seconds
    
    ''' simulating the model before turning on induce nodes  '''
    S_BS_OC.simulate(t_vec1,1)

#    t_vec2 = np.linspace(t_vec1[-1]/3600,6,1001)*3600 # seconds
#    
#    ''' simulating the model  '''
#    I_BS_IFFL12.simulate(t_vec2,2,dA=dA_add[i])    
#    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = S_BS_OC.output_concentration['GdA1']
    G3 = S_BS_OC.output_concentration['GdA3']
    G5 = S_BS_OC.output_concentration['GdA5']
    G7 = S_BS_OC.output_concentration['GdA7']
    G9 = S_BS_OC.output_concentration['GdA9']
    G10 = S_BS_OC.output_concentration['GdA10']

    ''' Plotting solution '''
    sim_t = S_BS_OC.sol.t
    
    fs = 13
    plt.figure()
  #  plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
       
    plt.subplot(3,3,1)   
  #  plt.title('State 1 to State 2')
    plt.plot(sim_t/3600,G9/G_tot[8],color=[0.5,0.5,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G10/G_tot[9],color=[1,0.75,0],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN5','TN6'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,10)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,4)    
    plt.plot(sim_t/3600,G1/G_tot[0],color=[0.75,0.9,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3/G_tot[2],color=[0,0.5,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    #plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,10)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    plt.subplot(3,3,7)
    plt.plot(sim_t/3600,G5/G_tot[4],color='b',linewidth=4)
    plt.plot(sim_t/3600,G7/G_tot[6],color=[0.4,0,0.4],linewidth=4)
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    plt.legend(['TN5','TN8'],frameon=0,loc=9)
    ax1 = plt.gca()
    ax1.set_xlim(0,10)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
''' 
EXAMPLE 17
###############################################################################
Misha's pulse network
-------------------------------------------------------------------------------
   -
###############################################################################
'''
if ex17 == 1:
    
# UPSTREAM AUTOREPRESSION 
    
    #          G1 G2 G3 G4  
    act_vec =  [1, 1, 2]
    prod_vec = [2,-1,-2]
    
    
    ''' initializing topology '''
    MP = GGM.GeneletNetwork(act_vec,prod_vec)
    
    # Define initial conditions
    #                  dA1  dA2
    dA_tot = np.array([250, 250]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2  G3  G4 
    G_tot = np.array([25, 25, 25]) # total genelet added
     
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #           G1 G2 G3 G4
    G_int_vec = [0,0,-1]
 
    ''' initializing initial conditions '''
    MP.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,250]) # default of 0 for all other initial conditions
    MP.plot_topology(show_rnas=0,layout='circular',plot_title = 'Upstream autorepression')
    
    t_vec1 = np.linspace(0,4,501)*3600 # seconds
    
    ''' simulating the model before turning on induce nodes  '''
    MP.simulate(t_vec1,1)

#    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = MP.output_concentration['GdA1']
    G3 = MP.output_concentration['GdA3']

    ''' Plotting solution '''
    sim_t = MP.sol.t
    
    fs = 13
    plt.figure(2)
  #  plt.suptitle('BS_IFFL1|2',fontsize=fs+2,weight='bold')
       
    plt.subplot(3,3,1)   
    plt.title('Upstream autorepression')
    plt.plot(sim_t/3600,G1/G_tot[0],color=[0,0.75,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G3/G_tot[2],color=[0,0,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,4)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
# NO UPSTREAM AUTOREPRESSION 
    
     #          G1 G2  
    act_vec =  [1, 2]
    prod_vec = [2,-2]
    
    ''' initializing topology '''
    MP2 = GGM.GeneletNetwork(act_vec,prod_vec)
    
    # Define initial conditions
    #                  dA1  dA2
    dA_tot = np.array([250, 250]) # total activator added (only dA2 is present to set that state)
    #                 G1  G2
    G_tot = np.array([25, 25]) # total genelet added
     
    # initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
    #           G1 G2 G3 G4
    G_int_vec = [0,-1]
 
    ''' initializing initial conditions '''
    MP2.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,250]) # default of 0 for all other initial conditions
    MP2.plot_topology(show_rnas=0,layout='circular',plot_title = 'No upstream autorepression')
    
    t_vec1 = np.linspace(0,4,501)*3600 # seconds
    
    ''' simulating the model before turning on induce nodes  '''
    MP2.simulate(t_vec1,1)

#    
    # pulling out solutions for plotting from output_concentration attribute
    G1 = MP2.output_concentration['GdA1']
    G2 = MP2.output_concentration['GdA2']

    ''' Plotting solution '''
    sim_t = MP2.sol.t
    
    plt.figure(2)
    plt.subplot(3,3,2)   
    plt.title('No upstream autorepression')
    plt.plot(sim_t/3600,G1/G_tot[0],color=[0,0.75,0],linewidth=4,linestyle='-')
    plt.plot(sim_t/3600,G2/G_tot[1],color=[0,0,1],linewidth=4,linestyle='-')
    plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
    plt.xlabel('time (hr)',fontsize=fs,weight='bold')
    plt.legend(['TN1','TN2'],frameon=0)
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,4)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_tick_params(which='both', size=5, width=2, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=2, direction='in', right='on')
    
    ''' Compiling genelet sequences to build the MP network '''
    input_file_loc = 'C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\Programming Club\\Python\Python 2020\\General genelet model\\git test\\PythonVersion\\all_genelet_sequences.xlsx'
    file_name='MP_genelet_sequences'
    path_name='C:\\Users\\sscha\\OneDrive - Johns Hopkins University\\Desktop\\'

    MP.compile_network_sequences(input_file_loc,save_file_name=file_name,save_path_name=path_name)